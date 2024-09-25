#!/bin/bash

##########################################################
# Genome Assembly of 16 House Finch and 1 Rosefinch Individuals
##########################################################

# Genomes of 16 House Finch and one Rosefinch individuals were assembled using
# the scripts below for each individual.

###############################################
# Step 1: Genome Assembly using Hifiasm
###############################################

# Assemble the assembly genome using Hifiasm with 30 threads
$hifiasm -o assembly -t 30 assembly_reads_2024Apr2.fa.gz

# Extract primary and haplotype assemblies from the GFA files
awk '/^S/{print ">"$2;print $3}' assembly.bp.p_ctg.gfa > assembly_prim.fa
awk '/^S/{print ">"$2;print $3}' assembly.bp.hap1.p_ctg.gfa > assembly_hap1.fa
awk '/^S/{print ">"$2;print $3}' assembly.bp.hap2.p_ctg.gfa > assembly_hap2.fa

###############################################
# Step 2: Assembly Quality Metrics
###############################################

# Generate assembly statistics for primary assembly
assembly-stats assembly_prim.fa assembly_prim.fa assembly_prim.fa > Stats_assembly_prim.txt

# Generate assembly statistics for haplotype 1
assembly-stats assembly_prim.fa assembly_hap1.fa assembly_hap1.fa > Stats_assembly_hap1.txt

# Generate assembly statistics for haplotype 2
assembly-stats assembly_prim.fa assembly_hap2.fa assembly_hap2.fa > Stats_assembly_hap2.txt

###############################################
# Step 3: Assess Genome Assembly Completeness using BUSCO
###############################################

# Run BUSCO analysis on the primary assembly
singularity exec --cleanenv busco.sif busco -m genome -i $assembly_prim -l aves_odb10 -c 20 -o assembly_prim -f

# Run BUSCO analysis on haplotype 1 assembly
singularity exec --cleanenv busco.sif busco -m genome -i $assembly_hap1 -l aves_odb10 -c 20 -o assembly_hap1 -f

# Run BUSCO analysis on haplotype 2 assembly
singularity exec --cleanenv busco.sif busco -m genome -i $assembly_hap2 -l aves_odb10 -c 20 -o assembly_hap2 -f

##########################################################
# Genome Gene Annotation
##########################################################

###############################################
# Step 4: De novo Transcript Discovery using IsoQuant
###############################################

# Use IsoQuant to discover transcripts from PacBio CCS Iso-Seq data
isoquant.py --reference $VGP_genome \
            --fastq $isoseq_brain $isoseq_eye $isoseq_liver $isoseq_heart \
            --data_type pacbio_ccs \
            --threads 20 \
            -o IsoQuant/isoquant_from_scrtach_4tissues

###############################################
# Step 5: Mapping RNA-seq Reads to VGP Genome using HISAT2
###############################################

# Build HISAT2 index of the VGP genome
hisat2-build $VGP_genome VGP_genome

# Map paired-end RNA-seq reads to the VGP genome and sort the alignments
hisat2 -x VGP_genome -p 10 -1 $RNA_1 -2 $RNA_2 --no-discordant | samtools sort -@ 3 -O BAM -o RNA.bam

###############################################
# Step 6: Transcript Assembly using StringTie
###############################################

# Short-read RNA was downloaded from:
# Zhang, Q., Hill, G. E., Edwards, S. V., & Backström, N. (2014).
# A house finch (Haemorhous mexicanus) spleen transcriptome reveals intra-
# and interspecific patterns of gene expression, alternative splicing, and
# genetic diversity in passerines. BMC Genomics, 15, 1-17.

# Assemble transcripts from RNA-seq alignments
stringtie -o transcriptome.gtf RNA.bam -p 30 -v

###############################################
# Step 7: Identify Coding Regions within Transcripts using TransDecoder
###############################################

# Convert GTF to FASTA format for transcripts
gtf_genome_to_cdna_fasta.pl transcriptome.gtf $VGP_genome > transcriptome.fasta

# Convert GTF to GFF3 format
gtf_to_alignment_gff3.pl transcriptome.gtf > transcriptome.gff3  # GTF to GFF

# Predict long open reading frames (ORFs) using TransDecoder
TransDecoder.LongOrfs -t transcriptome.fasta  # Generate best candidate ORF predictions

# Predict the likely coding regions
TransDecoder.Predict -t transcriptome.fasta  # Generate best candidate ORF predictions

# Map the predicted ORFs back to the genome coordinates
cdna_alignment_orf_to_genome_orf.pl \
    transcriptome.fasta.transdecoder.gff3 \
    transcriptome.gff3 \
    transcriptome.fasta > transcriptome_transdecoder.gff3

###############################################
# Step 8: Genome Annotation using BRAKER2
###############################################

## RNA-Seq Evidence Based Annotation

# Run BRAKER2 with RNA-seq data to predict genes
singularity exec --no-home \
                 --home /opt/gm_key \
                 --cleanenv \
                 --env AUGUSTUS_CONFIG_PATH=${PWD}/augustus_config \
                 braker_2.1.6_5-2022-05-04.sif \
                 braker.pl --cores 20 \
                            --bam=RNA.bam \
                            --workingdir=/n/holyscratch01/edwards_lab/bfang/Hofi/Assembly/Genome_Annotation/Braker/Braker_1/ \
                            --genome=$genome_masked_VGP \
                            --softmasking

## Protein Evidence Based Annotation

# Run BRAKER2 in epmode with protein data to predict genes
singularity exec --no-home \
                 --home /opt/gm_key \
                 --cleanenv \
                 --env AUGUSTUS_CONFIG_PATH=${PWD}/augustus_config \
                 braker_2.1.6_5-2022-05-04.sif \
                 braker.pl --cores 20 \
                            --epmode \
                            --genome=$genome_masked_VGP \
                            --workingdir=/n/holyscratch01/edwards_lab/bfang/Hofi/Assembly/Genome_Annotation/Braker/Braker_2 \
                            --prot_seq=$proteins \
                            --softmasking

###############################################
# Step 9: Genome Annotation Projection using TOGA
###############################################

# Using TOGA, we projected annotations of coding genes from Chicken (Gallus gallus, GCA_016699485.1),
# Zebra Finch (Taeniopygia guttata, GCA_003957565.4), and Common Canary (Serinus canaria, GCA_022539315.2).
# The scripts below are for Zebra Finch as an example.

## Convert FASTA to 2bit format (required by TOGA)

# For Zebra Finch genome
faToTwoBit Taeniopygia_guttata.bTaeGut1_v1.p.dna_sm.toplevel.fa.gz zebrafinch.2bit
twoBitInfo zebrafinch.2bit stdout | sort -k2rn > genome.chrom.sizes

# For VGP genome
faToTwoBit $VGP_genome_soft_masked.gz VGP.2bit
twoBitInfo VGP.2bit stdout | sort -k2rn > genome.chrom.sizes

## Convert GFF3 annotations to BED12 format
gff3ToGenePred <(gzip -d -c $zebrafinch.gff3.gz) temp.genePred
genePredToBed temp.genePred zebrafinch_gff3_anno.bed

# Remove 'transcript:' prefix from the BED file
sed -i 's,transcript:,,g' zebrafinch_gff3_anno.bed

## Run LASTZ and axtChain for each chromosome
for a in `ls Hofi_chrs_softmasked`; do
    CHR=$(basename $a .fa)
    VGP_chr_now=Hofi_chrs_softmasked/$CHR.fa

    echo "
    ### $CHR ###
    ## LASTZ Alignment
    lastz ../Taeniopygia_guttata.bTaeGut1_v1.p.dna_sm.toplevel.fa[multiple] $VGP_chr_now \\
    --step=100 --strand=both --ambiguous=n --inner=2000 --gappedthresh=6000 --gapped --identity=90 --format=axt > zebrafinch_hofi_$CHR.axt

    ## axtChain
    axtChain -faQ -faT -linearGap=loose zebrafinch_hofi_$CHR.axt ../Taeniopygia_guttata.bTaeGut1_v1.p.dna_sm.toplevel.fa $VGP_chr_now zebrafinch_hofi_$CHR.chain
    " >> LastZ_axtChain_Zebrafinch_Hofi.sh
done

# Run the generated script
# sbatch LastZ_axtChain_Zebrafinch_Hofi.sh

# Merge and sort all the chain files
chainMergeSort *.chain > zebrafinch_hofi_all.chain  # Merge and sort all the chains

## Run TOGA for gene annotation projection
toga.py \
    chains/zebrafinch_hofi_all.chain zebrafinch_gff3_anno.bed zebrafinch.2bit $VGP_2bit \
    -i zebrafinch_isoform.txt \
    --kt --pn zebrafinch_VGP_sm --nc ${path_to_nextflow_config_dir} \
    --cb 3,5 --cjn 500 --ms

###############################################
# Step 10: Gene Annotation using Miniprot
###############################################

# Index the VGP genome for Miniprot
$miniprot -t16 -d VGP.mpi $VGP_genome

# Change to the working directory
cd /n/holyscratch01/edwards_lab/bfang/Hofi/Assembly/Genome_Annotation/Miniprot

# Map proteins from Zebra Finch to the VGP genome
$miniprot -Iut16 --gff VGP.mpi zebra_finch_proteins.faa > VGP_Zebra_Finch.gff

# Map proteins from Chicken to the VGP genome
$miniprot -Iut16 --gff VGP.mpi chicken_proteins.faa > VGP_Chicken.gff

# Map proteins from Canary to the VGP genome
$miniprot -Iut16 --gff VGP.mpi canary_proteins.faa > VGP_Canary.gff

###############################################
# Step 11: Annotation Cleanup using AGAT
###############################################

## Fix overlapping genes
singularity exec $AGAT agat_sp_fix_overlaping_genes.pl \
    -gff GFF_merged.gff3 \
    -f $VGP_genome \
    -o GFF_rmdupGenes.gff3

## Remove features with duplicated locations
singularity exec $AGAT agat_sp_fix_features_locations_duplicated.pl \
    --gff $GFF_rmdupGenes \
    --out GFF_rmdupGenes_rmdup.gff3

## Keep the longest transcript for each gene
singularity exec $AGAT agat_sp_keep_longest_isoform.pl \
    --gff GFF_rmdupGenes_rmdup.gff3 \
    --out Hofi_LongestTranscript.gff3

###############################################
# Step 12: Assess Transcriptome Completeness using BUSCO
###############################################

singularity exec --cleanenv busco_v5.1.3_cv1.sif busco \
    -m transcriptome -c 20 -f \
    -i transcripts.fa \
    -l /n/holyscratch01/external_repos/INFORMATICS/BUSCO/aves_odb10 \
    -o transcriptome

##########################################################
# Genome Repeat Annotation
##########################################################

###############################################
# Step 13: Repeat Annotation using RepeatModeler and RepeatMasker
###############################################

## Step 13.1: De novo Repeat Identification using RepeatModeler

# Build a RepeatModeler database
singularity exec --cleanenv dfam_tetools.sif BuildDatabase -name VGP $VGP_genome -engine ncbi

# Run RepeatModeler with LTRStruct option
singularity exec --cleanenv dfam_tetools.sif RepeatModeler -LTRStruct -threads 25 -engine ncbi -database VGP 2>&1 | tee VGP.log

## Step 13.2: Mask Repeats using RepeatMasker

# Combine custom repeat libraries
cat VGP-families.fa repbase.fa > curated_repeats_VGP.fa

# Run RepeatMasker using the curated repeat library
singularity exec --cleanenv dfam_tetools.sif RepeatMasker -pa 25 -a -e ncbi -xsmall -lib dfam_tetools.fa -dir VGP $VGP_genome

###############################################
# Step 14: Detecting Segmental Duplications using BISER
###############################################

biser -o VGP -t 15 $VGP_softmasked --keep-contigs --keep-temp --gc-heap 10G

###############################################
# Step 15: Simple Repeat Finder (SRF) Analysis
###############################################

## Count k-mers using KMC
$kmc -fq -k151 -t16 -ci100 -cs100000 @VGP.txt VGP_count tmp_dir

## Dump k-mer counts
$kmc_dump VGP_count VGP_count.txt

## Find simple repeats using SRF
$srf -p VGP VGP_count.txt > srf_VGP.fa

## Align repeats to the genome using Minimap2
minimap2 -c -t 10 -N1000000 -f1000 -r100,100 <(k8 $srfutils enlong srf_VGP.fa) VGP.fq > VGP_srf_aln.paf

## Convert PAF to BED format
k8 $srfutils paf2bed VGP_srf_aln.paf > VGP_srf_aln.bed

## Generate abundance file
k8 $srfutils bed2abun VGP_srf_aln.bed > VGP_srf_aln.len

###############################################
# Step 16: Measuring Telomeric Repeats using seqtk telo
###############################################

# Generate a script to detect telomere lengths for multiple samples
for a in {VGP,AL,AZ,MA,NM,NY,OH,WA}; do
    echo "
    ##########
    ## ${a}
    ##########
    $seqtk telo <(cat ${a}_1_VGP_aln_no_ZW_downsampled20X_new.fq ${a}_2_VGP_aln_no_ZW_downsampled20X_new.fq) > ${a}.raw.bed 2> ${a}.raw.count" >> TelomereLength_detection.sh
done
