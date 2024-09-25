#!/bin/bash

###############################################
# PGGB | Pangenome Graph Builder
###############################################

#######################################################################
# PGGB | 16 House Finches + 1 Rosefinch + 1 VGP_BF (Bullfinch)
#######################################################################

###############################
# Step 1: Rename Sequences
###############################

# Use 'fastix' to rename sequences and compress them with 'bgzip'
# Then index the compressed FASTA files with 'samtools faidx'

$fastix -p 'sample#hap1#' sample_hap1.fa | bgzip -c -@ 25 > FASTA_sample_hap1_renamed.fa.gz && samtools faidx FASTA_sample_hap1_renamed.fa.gz
$fastix -p 'sample#hap2#' sample_hap2.fa | bgzip -c -@ 25 > FASTA_sample_hap2_renamed.fa.gz && samtools faidx FASTA_sample_hap2_renamed.fa.gz

####################################
# Step 2: Sequence Partitioning
# Convert FASTA to PAF using wfmash
####################################

# Define the reference genome variable
FASTA_VGP_prim_renamed=FASTA_VGP_prim_renamed.fa.gz

# Loop over each sample and generate PAF files for both haplotypes
for a in {VGPBF,RF,AL_1,AL_2,AZ_1,AZ_2,CA_1,CA_2,MA_1,MA_2,NM_1,NM_2,NY_1,NY_2,OH_1,OH_2,WA_1,WA_2}; do
    PATH_PAF_hap1=${a}_hap1.vs.VGP.paf
    PATH_PAF_hap2=${a}_hap2.vs.VGP.paf

    # Append commands to 'Code_wfmash.sh' for each sample
    echo "
######################
## ${a}
######################
singularity exec --cleanenv /n/home00/bfang/programs/pggb/PGGB_May25_2023.sif \\
wfmash $FASTA_VGP_prim_renamed $FASTA_${a}_hap1_renamed -s 50k -p 90 -N -m -t 8 > $PATH_PAF_hap1

singularity exec --cleanenv /n/home00/bfang/programs/pggb/PGGB_May25_2023.sif \\
wfmash $FASTA_VGP_prim_renamed $FASTA_${a}_hap2_renamed -s 50k -p 90 -N -m -t 8 > $PATH_PAF_hap2
" >> Code_wfmash.sh
done

##############################################################
# Step 3: Collect Unmapped Contigs and Remap in Split Mode
##############################################################

# Loop over all haplotype FASTA files
ls *hap*gz | while read FASTA; do
    NAME=$(basename $FASTA _renamed.fa.gz | sed 's,FASTA_,,g')
    PATH_PAF=/partitioning/$NAME.vs.VGP.paf
    PATH_NO_SPLIT_PAF=/partitioning/$NAME.vs.ref.no_split.paf

    # Find unaligned contigs
    comm -23 <(cut -f1 $FASTA.fai | sort) <(cut -f1 $PATH_PAF | sort) > $NAME.unaligned.txt
    wc -l $NAME.unaligned.txt

    # If there are unaligned contigs, remap them
    if [[ $(wc -l $NAME.unaligned.txt | cut -f1 -d' ') != 0 ]]; then
        echo "
######################
## $NAME
######################
samtools faidx $FASTA $(tr '\n' ' ' < $NAME.unaligned.txt) > $NAME.unaligned.fa

# Index unaligned FASTA
samtools faidx $NAME.unaligned.fa

# Map unaligned FASTA back to reference genome using wfmash
wfmash $FASTA_VGP_prim_renamed $NAME.unaligned.fa -s 50k -p 90 -m -t 8 > $PATH_NO_SPLIT_PAF
" >> Code_wfmash_18IDs_NO_SPLIT.sh
    fi
done

#######################################################################
# Step 4: Collect Best Mappings for Split Rescue Attempts
#######################################################################

# Contigs that did not map are partitioned using a split mapping approach,
# requiring 90% identity over 50kb to seed the mappings.

cd /partitioning
ls *.vs.ref.no_split.paf | while read PAF; do
    cat $PAF | awk -v OFS='\t' '{ print $1,$11,$0 }' | sort -k1,1 -k2,2nr | awk -v OFS='\t' '$1 != last { print($3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15); last = $1; }'
done > rescues.paf

#################################
# Step 5: Process Each Chromosome
#################################

cd /parts
module load python
source activate fasrc

# List of good chromosomes
good_chrs=Good_Chrs_2022-07-28.txt
DIR_PARTITIONING=/partitioning

# For each chromosome, collect contigs
for CHR in $(cat $good_chrs); do
    echo "$CHR"

    awk ' $6 ~ "'$CHR'$" ' $(ls $DIR_PARTITIONING/*.vs.VGP.paf | grep -v unaligned | sort; echo $DIR_PARTITIONING/rescues.paf) | cut -f1 | sort | uniq > $CHR.contigs
done

# Extract sequences and create pan-genome FASTA
for CHR in $(cat $good_chrs); do
    echo "$CHR"
    samtools faidx $FASTA_VGP_prim_renamed $CHR > $CHR.pan.fa

    ls /n/holylfs04/LABS/edwards_lab/Lab/bfang/Genomes_assembled/*hap*gz | while read FASTA; do
        NAME=$(basename $FASTA _renamed.fa.gz)
        echo "$CHR $NAME"

        samtools faidx $FASTA $(comm -12 <(cut -f1 $FASTA.fai | sort) <(sort $CHR.contigs)) >> $CHR.pan.fa
    done
    bgzip -@ 10 $CHR.pan.fa && samtools faidx $CHR.pan.fa.gz
done

#################################################
# Step 6: Pangenome Graph Building per Chromosome
#################################################

cd /graphs
DIR_PARTS=/parts

# Build pangenome graph for each chromosome using PGGB
pggb -i $DIR_PARTS/$CHR.pan.fa.gz -o $CHR.pan \
    -n 37 \  # Number of haplotypes
    -p 90 \  # Minimum average nucleotide identity for segments
    -s 50k \ # Segment length for mapping
    -k 79 \  # Filter exact matches below this length
    -T 25 -t 25

################################
# Step 7: Check size of GFA file
################################
cd /graphs
ls VGP*pan/*fix.gfa -rt1g | awk '{print $4, $8}' > size.txt


#######################################
# Step 8: Call VCF Using vg deconstruct
#######################################

# Rename 'VGPBF' to 'CABF' to prevent VG from treating it as the reference genome.

for CHR in $(cat $good_chrs); do
    PATH_GFA=$PATH_out
    PATH_SED_GFA=$CHR.pan.sed.gfa

    NAME=$(echo $CHR | sed 's,VGP#prim#,,g')
    echo "
########################
## $NAME
sed 's/#/-/' $PATH_GFA | sed 's/-hap/#/' > $PATH_SED_GFA

singularity exec --cleanenv /n/home00/bfang/programs/pggb/PGGB_Dec27_2023.sif \
vg deconstruct -P \"VGP-prim\" -H '#' -e -a -t 10 $PATH_SED_GFA > $PATH_VCF_Pangenome/VCF_${NAME}_REF_VGP.vcf
" >> Code_PGGB_CallVCF_Normalization.sh
done

###############################################
# Step 9: Normalization with VCFbub and VCFwave
###############################################

# Compress VCF files for VCFbub
for a in *.vcf; do bgzip $a -@ 10; done

# Run VCFbub and VCFwave
for CHR in $(cat $good_chrs | sed 's,VGP#prim#,,g'); do
echo "
########################
## $CHR
singularity exec --cleanenv PGGB_Dec27_2023.sif \
vcfbub -l 0 -a 100000 --input VCF_${CHR}_REF_VGP.vcf.gz > vcfbub_Jan2024/VCF_${CHR}_REF_VGP_bub.vcf

singularity exec --cleanenv PGGB_Dec27_2023.sif \
vcfwave -I 1000 -t 20 vcfbub_Jan2024/VCF_${CHR}_REF_VGP_bub.vcf > vcfwave_Jan2024/VCF_${CHR}_REF_VGP_bub_wave.vcf
" >> Code_Normalize_VCF_Jan2024.sh
done

################################
# Step 10: Merge All Chromosomes
################################

# Prepare reference FASTA for normalization
zcat FASTA_VGP_prim_renamed.fa.gz | sed 's,VGP#prim#,VGP-prim#,g' | bgzip --threads 30 -c > FASTA_VGP_prim_renamed_for_bcftools.fa.gz
samtools faidx FASTA_VGP_prim_renamed_for_bcftools.fa.gz

# Concatenate VCF files
cd /vcfwave_Jan2024
bcftools concat -f files_vcf_gz.txt | bgzip --threads 30 -c > pggb_merged.vcf.gz

# Normalize using bcftools
cd /vcfwave_Jan2024

# Fix ploidy
bcftools +fixploidy pggb_merged.vcf.gz > pggb_merged_fixploidy.vcf

# Fill tags
bcftools +fill-tags pggb_merged_fixploidy.vcf -Ob -o pggb_merged_cleanup.bcf -- -t AN,AC,AF

# Normalize variants (this may introduce some INV=1 INFO lines)
bcftools norm -m+both -f FASTA_VGP_prim_renamed_for_bcftools.fa.gz pggb_merged_cleanup.bcf > pggb_merged_normalized.vcf

# Clean and sort
bcftools sort pggb_merged_normalized.vcf | bcftools +fill-tags -- -t AN,AC,AF,F_MISSING | bcftools view -o pggb_cleaned_final.vcf.gz -Oz

# Remove "AT" in the INFO field to reduce VCF loading time
bcftools annotate --remove INFO/AT pggb_cleaned_final.vcf.gz -o pggb_cleaned_final_noAT.vcf.gz -Oz --threads 30

# Compress and delete intermediate VCFs to save space
bcftools view -o pggb_merged_fixploidy.vcf.gz -Oz pggb_merged_fixploidy.vcf
bcftools view -o pggb_merged_normalized.vcf.gz -Oz pggb_merged_normalized.vcf
rm -f pggb_merged_fixploidy.vcf
rm -f pggb_merged_normalized.vcf

###############################
# Step 11: Graph Statistics
###############################
cd /graphs_2024
for a in *pan/*final.gfa; do
    NAME=$(basename "$a" .pan.fa.gz.dac1d73.c2fac19.754753f.smooth.final.gfa | sed 's,VGP#prim#,,g')
    echo "singularity exec --cleanenv PGGB_Dec27_2023.sif odgi stats -i ${a} -t 30 -S > Stats_communities_2024Mar20/${NAME}.txt" >> Sbatch_stats_communities_2024Mar20.sh
done

#######################################################################
# Minigraph | Sequence-to-Graph Mapper and Graph Constructor
#######################################################################

# Build the minigraph
$minigraph -cxggs -t 25 \
    $FASTA_VGP_prim \
    $FASTA_VGPBF_hap1_renamed \
    $FASTA_VGPBF_hap2_renamed \
    $FASTA_NY_1_hap1_renamed \
    $FASTA_NY_1_hap2_renamed \
    $FASTA_NY_2_hap1_renamed \
    $FASTA_NY_2_hap2_renamed \
    $FASTA_OH_1_hap1_renamed \
    $FASTA_OH_1_hap2_renamed \
    $FASTA_OH_2_hap1_renamed \
    $FASTA_OH_2_hap2_renamed \
    $FASTA_WA_1_hap1_renamed \
    $FASTA_WA_1_hap2_renamed \
    $FASTA_WA_2_hap1_renamed \
    $FASTA_WA_2_hap2_renamed \
    $FASTA_AL_1_hap1_renamed \
    $FASTA_AL_1_hap2_renamed \
    $FASTA_AL_2_hap1_renamed \
    $FASTA_AL_2_hap2_renamed \
    $FASTA_MA_1_hap1_renamed \
    $FASTA_MA_1_hap2_renamed \
    $FASTA_MA_2_hap1_renamed \
    $FASTA_MA_2_hap2_renamed \
    $FASTA_NM_1_hap1_renamed \
    $FASTA_NM_1_hap2_renamed \
    $FASTA_NM_2_hap1_renamed \
    $FASTA_NM_2_hap2_renamed \
    $FASTA_AZ_1_hap1_renamed \
    $FASTA_AZ_1_hap2_renamed \
    $FASTA_AZ_2_hap1_renamed \
    $FASTA_AZ_2_hap2_renamed \
    $FASTA_CA_1_hap1_renamed \
    $FASTA_CA_1_hap2_renamed \
    $FASTA_CA_2_hap1_renamed \
    $FASTA_CA_2_hap2_renamed \
    $FASTA_RF_hap1_renamed \
    $FASTA_RF_hap2_renamed \
    > Minigraph_RefVGP_37haps_Jun2_2023.gfa

# Convert GFA to FASTA
$gfatools gfa2fa Minigraph_RefVGP_37haps_Jun2_2023.gfa > Minigraph_RefVGP_37haps_Jun2_2023.fa

# Extract localized structural variations
$gfatools bubble Minigraph_RefVGP_37haps_Jun2_2023.gfa > Minigraph_RefVGP_37haps_Jun2_2023.bed

###############################################
# Minigraph | Convert BED to VCF
###############################################

# Merge SV BED files and convert to VCF
paste \
SVs_Minigraph_RefVGP_RF_hap1_.bed \
SVs_Minigraph_RefVGP_RF_hap2_.bed \
SVs_Minigraph_RefVGP_VGPBF_hap1_.bed \
SVs_Minigraph_RefVGP_VGPBF_hap2_.bed \
SVs_Minigraph_RefVGP_CA_1_hap1_.bed \
SVs_Minigraph_RefVGP_CA_1_hap2_.bed \
SVs_Minigraph_RefVGP_CA_2_hap1_.bed \
SVs_Minigraph_RefVGP_CA_2_hap2_.bed \
SVs_Minigraph_RefVGP_WA_1_hap1_.bed \
SVs_Minigraph_RefVGP_WA_1_hap2_.bed \
SVs_Minigraph_RefVGP_WA_2_hap1_.bed \
SVs_Minigraph_RefVGP_WA_2_hap2_.bed \
SVs_Minigraph_RefVGP_NM_1_hap1_.bed \
SVs_Minigraph_RefVGP_NM_1_hap2_.bed \
SVs_Minigraph_RefVGP_NM_2_hap1_.bed \
SVs_Minigraph_RefVGP_NM_2_hap2_.bed \
SVs_Minigraph_RefVGP_AZ_1_hap1_.bed \
SVs_Minigraph_RefVGP_AZ_1_hap2_.bed \
SVs_Minigraph_RefVGP_AZ_2_hap1_.bed \
SVs_Minigraph_RefVGP_AZ_2_hap2_.bed \
SVs_Minigraph_RefVGP_MA_1_hap1_.bed \
SVs_Minigraph_RefVGP_MA_1_hap2_.bed \
SVs_Minigraph_RefVGP_MA_2_hap1_.bed \
SVs_Minigraph_RefVGP_MA_2_hap2_.bed \
SVs_Minigraph_RefVGP_AL_1_hap1_.bed \
SVs_Minigraph_RefVGP_AL_1_hap2_.bed \
SVs_Minigraph_RefVGP_AL_2_hap1_.bed \
SVs_Minigraph_RefVGP_AL_2_hap2_.bed \
SVs_Minigraph_RefVGP_NY_1_hap1_.bed \
SVs_Minigraph_RefVGP_NY_1_hap2_.bed \
SVs_Minigraph_RefVGP_NY_2_hap1_.bed \
SVs_Minigraph_RefVGP_NY_2_hap2_.bed \
SVs_Minigraph_RefVGP_OH_1_hap1_.bed \
SVs_Minigraph_RefVGP_OH_1_hap2_.bed \
SVs_Minigraph_RefVGP_OH_2_hap1_.bed \
SVs_Minigraph_RefVGP_OH_2_hap2_.bed \
| /n/home00/bfang/programs/minigraph/misc/mgutils.js merge - > Minigraph_RefVGP_37haps.bed

# Convert merged BED to VCF
./mgutils.js merge2vcf Minigraph_RefVGP_37haps.bed > Minigraph_RefVGP_37haps.vcf

#######################################################################
# Pangene | Constructing a Pangenome Gene Graph
#######################################################################

# Step 1: Rename FASTA files
for a in {VGPBF,AL_1,AL_2,AZ_1,AZ_2,MA_1,MA_2,NM_1,NM_2,NY_1,NY_2,OH_1,OH_2,WA_1,WA_2,CA_1,CA_2,RF}; do
echo "
$fastix -p '${a}#hap1#' $FASTA${a}_hap1 > FASTA_${a}_hap1_renamed.fasta
$fastix -p '${a}#hap2#' $FASTA${a}_hap2 > FASTA_${a}_hap2_renamed.fasta
" >> code_rename_Oct21_2023.sh
done

###############################################
# Step 2: Run Miniprot on Renamed Haplotypes
###############################################

miniprot=/n/home00/bfang/programs/miniprot/miniprot
protein=Hofi_NCBI_protein_Oct2023_renamed.faa

for a in {VGPBF,RF,AL_1,AL_2,AZ_1,AZ_2,MA_1,MA_2,NM_1,NM_2,NY_1,NY_2,OH_1,OH_2,WA_1,WA_2,CA_1,CA_2}; do
echo "
$miniprot -t16 -d ${a}_hap1.mpi $FASTA${a}_hap1_renamed
$miniprot -Iut16 --gff ${a}_hap1.mpi $protein > ${a}_hap1.gff
$miniprot --outs=0.97 --no-cs -Iut16 ${a}_hap1.mpi $protein > ${a}_hap1.paf

$miniprot -t16 -d ${a}_hap2.mpi $FASTA${a}_hap2_renamed
$miniprot -Iut16 --gff ${a}_hap2.mpi $protein > ${a}_hap2.gff
$miniprot --outs=0.97 --no-cs -Iut16 ${a}_hap2.mpi $protein > ${a}_hap2.paf
" >> Codes_Miniprot_36haps_Oct21_2023.sh
done

##################################################
# Step 3: Build Pangenome Gene Graph Using Pangene
##################################################

pangene=/n/home00/bfang/programs/pangene/pangene

$pangene \
VGP_prim.paf \
VGPBF_hap2.paf VGPBF_hap1.paf \
WA_2_hap2.paf WA_2_hap1.paf \
WA_1_hap2.paf WA_1_hap1.paf \
OH_2_hap2.paf OH_2_hap1.paf \
OH_1_hap2.paf OH_1_hap1.paf \
NY_2_hap2.paf NY_2_hap1.paf \
NY_1_hap2.paf NY_1_hap1.paf \
NM_2_hap2.paf NM_2_hap1.paf \
NM_1_hap2.paf NM_1_hap1.paf \
MA_2_hap2.paf MA_2_hap1.paf \
MA_1_hap2.paf MA_1_hap1.paf \
CA_2_hap2.paf CA_2_hap1.paf \
CA_1_hap2.paf CA_1_hap1.paf \
AZ_2_hap2.paf AZ_2_hap1.paf \
AZ_1_hap2.paf AZ_1_hap1.paf \
AL_2_hap2.paf AL_2_hap1.paf \
AL_1_hap2.paf AL_1_hap1.paf \
RF_hap2.paf RF_hap1.paf \
> pangene_36haps_Oct21_2023.gfa
