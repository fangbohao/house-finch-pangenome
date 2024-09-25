#!/bin/bash

############################################
# SNP Calling and Filtering for House Finch
############################################

# This script performs SNP calling and filtering on RAD, WGS, and HiFi samples.
# It includes mapping reads to the VGP reference genome, SNP calling with bcftools and longshot,
# merging VCF files, filtering variants, and calculating heterozygosity.

#############################
# RAD and WGS Sample Analysis
#############################

# Index the VGP reference genome using bowtie-build
bowtie-build -f bHaeMex1.pri.cur.20220203.fasta bHaeMex1.pri.cur.20220203  # Index VGP reference genome

# Mapping reads (FASTQ files) to the VGP genome
while read ID FASTQ_1 FASTQ_2 PATH; do
    echo "
#########
# ${ID}
#########
singularity exec --cleanenv /n/home00/bfang/programs/bowtie/Bowtie2.sif \\
bowtie2 -x ${GENOME} -1 ${FASTQ_1} -2 ${FASTQ_2} --very-sensitive-local -I 149 -p 15 | \\
samtools view -b | samtools sort -o ${PATH}/${ID}.sorted.bam

# Add read group information
$picard AddOrReplaceReadGroups INPUT=${PATH}/${ID}.sorted.bam OUTPUT=${PATH}/${ID}.bam \\
    RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${ID} && \\
samtools index ${PATH}/${ID}.bam

# Mark duplicates and fix mate information
$picard MarkDuplicates \\
    INPUT=${PATH}/${ID}.bam \\
    OUTPUT=${PATH}/${ID}_marked.bam \\
    METRICS_FILE=${PATH}/${ID}_marked.metrics.txt \\
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

$picard FixMateInformation \\
    INPUT=${PATH}/${ID}_marked.bam \\
    OUTPUT=${PATH}/${ID}_marked_fixmate.bam \\
    SO=coordinate CREATE_INDEX=true
" >> Codes_mapping_RAD_WGS.sh
done < SAMPLE.txt

###########################
# SNP Calling with bcftools
###########################
# Generate SNP calling script for all chromosomes (1 to 39)
for a in `seq 1 39`; do
    echo "
bcftools mpileup -r SUPER_${a} -d 10000 -q 30 -Q 30 -f \$FASTA_VGP_prim --ignore-RG \\
    -b BAM_list_141IDs.txt --threads 20 | \\
    bcftools call --threads 20 -vm -Oz -o VCFs_Chrs_RAD2016/Hofi_Allison2016_141IDs_Chr${a}.vcf.gz
" >> Script_Call_SNPs_39Chrs_RAD2016_Nov13_2023.sh
done

# Index VCF files
bcftools index *.vcf.gz

# Concatenate VCF files of autosomes
bcftools concat -f List_VCFs_39Chrs_RAD_WGS.txt | bgzip --threads 25 -c > Hofi_RAD_WGS_39Chrs.vcf.gz

######################
# HiFi Sample Analysis
######################

# Mapping and SNP calling for each of the 16 individuals
# Example for AZ_1:

# Mapping HiFi reads to VGP genome using minimap2
minimap2 -ax map-hifi $VGP_genome_autosomes $HiFi_AZ_1 -t 30 | \\
    samtools sort -m4G -@ 30 -O BAM -o Read_to_VGP/AZ_1_VGP_aln.bam

# SNP calling using Longshot
longshot --bam Read_to_VGP/AZ_1_VGP_aln.bam --ref $VGP_genome_autosomes --out AZ_1_VGP.vcf

# Index VCF files
bcftools index *vcf.gz

############################
# Merging and Filtering SNPs
############################

# Merge VCF files from 135 individuals (RAD + WGS + HiFi)
bcftools merge -l vcf_135IDs.list --threads 20 -Oz -o Hofi_HiFi_RAD_WGS_autosomes_135IDs.vcf.gz

# SNP filtering using VCFtools
vcftools --gzvcf Hofi_HiFi_RAD_WGS_autosomes_135IDs.vcf.gz \\
    --min-alleles 2 --max-alleles 2 --remove-indels \\
    --max-missing 0.95 \\
    --mac 2 \\
    --recode --recode-INFO-all --stdout | \\
    bgzip --threads 15 > Hofi_HiFi_RAD_WGS_autosomes_135IDs_missing96_mac2.vcf.gz

# Calculate heterozygosity (genome-wide and inversion-associated SNPs)
vcftools --gzvcf Hofi_HiFi_RAD_WGS_autosomes_135IDs_missing96_mac2.vcf.gz --het --out Hofi_HiFi_RAD_WGS_autosomes_135IDs_missing96_mac2
vcftools --gzvcf Hofi_HiFi_RAD_WGS_autosomes_135IDs_missing96_mac2_INV.vcf.gz --het --out Hofi_HiFi_RAD_WGS_autosomes_135IDs_missing96_mac2_INV

###############################
# Runs of Homozygosity Analysis
###############################

# Convert VCF to PLINK format
plink --vcf $VCF_HiFi_biallelicSNPs_autosomes_095 --vcf-half-call m --make-bed \\
    --out VCF_HiFi_biallelicSNPs_autosomes_095 --allow-extra-chr

# Identify runs of homozygosity
plink --bfile VCF_HiFi_biallelicSNPs_autosomes_095 --homozyg \\
    --out VCF_HiFi_biallelicSNPs_autosomes_095 --allow-extra-chr \\
    --homozyg-window-snp 50 --homozyg-snp 50 --homozyg-kb 10 \\
    --homozyg-gap 300 --homozyg-density 200 --homozyg-window-missing 2 \\
    --homozyg-het 2 --homozyg-window-het 2

################################################################
# Individual Heterozygosity based on SNPs, Insertions, Deletions
################################################################

# Calculate heterozygosity using DipCall for each individual
for a in {AL_1,AZ_1,CA_1,MA_1,NM_1,NY_1,OH_1,WA_1,AL_2,AZ_2,CA_2,MA_2,NM_2,NY_2,OH_2,WA_2}; do
    echo "
##########
## ${a}
##########
dipcall.kit/run-dipcall ${a} \$VGP_genome \$FASTA_${a}_hap1 \$FASTA_${a}_hap2 > ${a}/${a}.mak

cd Dipcall/${a}
make -j2 -f ${a}.mak 2>&1 | tee ${a}.log

# Get phased VCF
/n/home00/bfang/programs/dipcall/dipcall.kit/k8 dipcall.kit/dipcall-aux.js vcfpair -a ${a}.pair.vcf.gz | \\
dipcall/dipcall.kit/htsbox bgzip > ${a}.dip.vcf.gz

# Exclude sex chromosomes and calculate statistics
grep -v '^SUPER_Z\|^SUPER_W' ${a}.dip.bed > ${a}.dip_nosex.bed  # Exclude sex chromosomes
dipcall/dipcall.kit/dipcall-aux.js dipsum ${a}.dip_nosex.bed ${a}.dip.vcf.gz > Statistics_${a}.txt
" >> Code_dipcall_heterozygosity.sh
done
