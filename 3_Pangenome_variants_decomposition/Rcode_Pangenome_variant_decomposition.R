library(data.table)
library(ggplot2)
library(cowplot)
library(parallel)
library(bedtoolsr)

setwd("/n/holyscratch01/edwards_lab/bfang/Hofi/Assembly/Pangenome/SFS_2024")
CORES <- 30

#########################
### Load the VCF file ###
#########################
vcf_raw <- fread("pggb_cleaned_final_noAT.vcf.gz", skip = "#CHROM", nThread = CORES)
print("1 | VCF loaded!")

# Split the VCF data by chromosome
vcf_raw <- split(vcf_raw, by = "#CHROM")

vcf_info_38chrs <- rbindlist(lapply(1:length(vcf_raw), function(x) {
  # Extract data for the current chromosome
  vcf_raw_now <- vcf_raw[[x]]
  chr_now <- gsub("VGP-prim#", "", names(vcf_raw[x]))
  
  #############################
  ### Calculate Allele Info ###
  #############################
  
  # Calculate the number of alleles (N_allele)
  N_allele <- vcf_raw_now[,.(all=paste(REF,ALT,sep = ","))][,mclapply(.SD, function(x) sapply(strsplit(x,","),function(y) length(y)),mc.cores = CORES)][,setnames(.SD,"N_allele")]
  
  # Calculate the length of the reference allele (N_bp_ref)
  N_bp_ref <- vcf_raw_now[, .(N_bp_ref = nchar(REF))]
  
  # Calculate the maximum base pair difference among alleles (N_bp_diff)
  N_bp_diff <- vcf_raw_now[, .(all = paste(REF, ALT, sep = ","))][, mclapply(.SD, function(x) sapply(strsplit(x, ","), function(y) max(nchar(y) - min(nchar(y)))), mc.cores = CORES)][, setnames(.SD, "N_bp_diff")]
  
  # Compile VCF information into a data.table
  vcf_info <- vcf_raw_now[,c(1,2,4,5)][,setnames(.SD,old="#CHROM",new="CHROM")][,CHROM:=gsub("VGP-prim#","",CHROM)]
  vcf_info <- data.table(vcf_info, N_allele, N_bp_ref, N_bp_diff)
  
  ######################################
  ### Determine General Variant Type ###
  ######################################
  # Classify variants based on allele length differences
  vcf_info[N_bp_diff == 0 & N_bp_ref == 1, type := "SNP"]
  vcf_info[N_bp_diff == 0 & N_bp_ref < 50 & N_bp_ref != 1, type := "MNP"]
  vcf_info[N_bp_diff == 0 & N_bp_ref >= 50, type := "SV_complex"]
  vcf_info[N_bp_diff > 0 & N_bp_diff < 50, type := "INDEL"]
  vcf_info[N_bp_diff > 0 & N_bp_diff >= 50, type := "SV"]
  
  # Check the counts of each variant type
  vcf_info[, .N, by = type]
  
  # Assign variant length (N_bp)
  vcf_info[type == "SNP", N_bp := 1]
  vcf_info[type != "SNP", N_bp := N_bp_diff]
  
  #############################
  ### Process Genotype Data ###
  #############################
  
  ID_haps <- unlist(lapply(colnames(vcf_raw_now[, 10:27]), function(x) list(paste0(x, "_hap1"), paste0(x, "_hap2"))))
  
  # Extract genotype data
  geno <- vcf_raw_now[, 10:27]
  
  # Convert phased genotypes from "/" to "|"
  geno <- geno[, mclapply(.SD, function(x) gsub("/", "|", x), mc.cores = CORES)]
  
  # Replace missing genotypes with ".|."
  geno[geno == "."] <- ".|."
  
  # Ensure all genotypes have two alleles
  geno <- geno[, mclapply(.SD, function(x) ifelse(!x %like% "\\|", paste0(x, "|."), x), mc.cores = CORES)]
  
  # Split genotypes into haplotype 1 and haplotype 2
  geno_hap1 <- geno[, mclapply(.SD, function(x) sapply(strsplit(x, "\\|"), function(y) y[1]), mc.cores = CORES)][
    , setnames(.SD, new = paste0(names(.SD), "_hap1"))
  ]
  geno_hap2 <- geno[, mclapply(.SD, function(x) sapply(strsplit(x, "\\|"), function(y) y[2]), mc.cores = CORES)][
    , setnames(.SD, new = paste0(names(.SD), "_hap2"))
  ]
  
  # Combine haplotype data and order columns
  geno <- cbind(geno_hap1, geno_hap2)[, ..ID_haps]
  
  print(paste0("2 | Genotype converted!", chr_now))

  # Determine invariant sites among House Finch haplotypes
  invariants_tag <- geno[
    , .(invariants_HF = pmax(
      AL_1_hap1, AL_1_hap2, AL_2_hap1, AL_2_hap2,
      AZ_1_hap1, AZ_1_hap2, AZ_2_hap1, AZ_2_hap2,
      CA_1_hap1, CA_1_hap2, CA_2_hap1, CA_2_hap2,
      MA_1_hap1, MA_1_hap2, MA_2_hap1, MA_2_hap2,
      NM_1_hap1, NM_1_hap2, NM_2_hap1, NM_2_hap2,
      NY_1_hap1, NY_1_hap2, NY_2_hap1, NY_2_hap2,
      OH_1_hap1, OH_1_hap2, OH_2_hap1, OH_2_hap2,
      WA_1_hap1, WA_1_hap2, WA_2_hap1, WA_2_hap2
    ) == 0)
  ]
  
  # Add invariants tag to VCF info
  vcf_info <- data.table(vcf_info, invariants_tag)
  
  ################################
  ### Calculate Derived Allele ###
  ################################
  
  # Create a backup of the genotype data
  geno_1 <- copy(geno)
  
  # Determine polarization using Rosefinch as the outgroup
  ancestry <- geno_1[, .(RF_hap1, RF_hap2, RF = paste0(RF_hap1, "|", RF_hap2))]
  ancestry[, polarized := ifelse(RF %in% c("0|0", "1|1", "1|.", "0|."), TRUE, FALSE)]
  vcf_info[, ancestry_polarized := ancestry$polarized]
  
  # Calculate derived alleles using RF_hap1 as the ancestral allele
  tag_outgroup <- which(names(geno_1) == "RF_hap1")
  geno_derived <- data.table(do.call(
    cbind,
    mclapply(1:ncol(geno_1), function(x) {
      ifelse(geno_1[, x, with = FALSE] == geno_1[, tag_outgroup, with = FALSE], 0, 1)
    }, mc.cores = CORES)
  ))
  geno_derived[geno_1 == "."] <- "."
  
  print(paste0("3 | Derived allele count calculated!", chr_now))
  
  ##############################################
  ### Calculate Derived Allele Counts (DAC) ###
  ##############################################
  
  # Note: Derived Allele Count is valid only for polarized bi-allelic variants
  
  # Define population groups
  ID_west <- c("CA", "NM", "AZ", "WA")
  ID_east <- c("AL", "OH", "NY", "MA")
  columns_west <- grep(paste(ID_west, collapse = "|"), colnames(geno_derived), value = TRUE)
  columns_east <- grep(paste(ID_east, collapse = "|"), colnames(geno_derived), value = TRUE)
  columns_all <- c(columns_west, columns_east)
  
  # Calculate DAC for each population
  DAC_all <- data.table(
    DAC_west = geno_derived[, ..columns_west][, rowSums(.SD == 1)],
    DAC_east = geno_derived[, ..columns_east][, rowSums(.SD == 1)],
    DAC_all = geno_derived[, ..columns_all][, rowSums(.SD == 1)],
    Geno_RF_hap1 = geno_1$RF_hap1  # Record outgroup raw genotype
  )
  
  # Calculate missing data counts (MC)
  missing_all <- data.table(
    MC_west = geno_derived[, ..columns_west][, rowSums(.SD == ".")],
    MC_east = geno_derived[, ..columns_east][, rowSums(.SD == ".")],
    MC_all = geno_derived[, ..columns_all][, rowSums(.SD == ".")]
  )[, MC_ratio := MC_all / 32]
  
  # Combine DAC and missing data with VCF info
  vcf_info <- data.table(vcf_info, DAC_all, missing_all)
  
  ###################################
  ### Allele Count per Population ###
  ###################################
  
  # Calculate the number of unique alleles per site for each population
  vcf_info$allele_count_West <- apply(geno_1[, ..columns_west],1,function(x) length(unique(x[x!="."])))
  vcf_info$allele_count_East <- apply(geno_1[, ..columns_east],1,function(x) length(unique(x[x!="."])))
  vcf_info$allele_count_All <- apply(geno_1[, ..columns_all],1,function(x) length(unique(x[x!="."])))
  
  
  ##################################
  ### Determine Variant Subtypes ###
  ##################################
  # Classify SNPs and MNPs
  vcf_info[type == "SNP", subtype := "SNP"]
  vcf_info[type == "MNP", subtype := "MNP"]
  
  # Classify INDELs (2-49 bp), valid only for polarized bi-allelic sites
  vcf_info[type == "INDEL" & Geno_RF_hap1 == 0 & nchar(REF) == 1, subtype := "INS"]
  vcf_info[type == "INDEL" & Geno_RF_hap1 == 0 & nchar(ALT) == 1, subtype := "DEL"]
  vcf_info[type == "INDEL" & Geno_RF_hap1 == 1 & nchar(REF) == 1, subtype := "DEL"]
  vcf_info[type == "INDEL" & Geno_RF_hap1 == 1 & nchar(ALT) == 1, subtype := "INS"]
  vcf_info[type == "INDEL" & is.na(subtype), subtype := "INDEL_complex"]
  
  # Classify SVs (>= 50 bp), valid only for polarized bi-allelic sites
  vcf_info[type == "SV" & Geno_RF_hap1 == 0 & nchar(REF) == 1, subtype := "SV_INS"]
  vcf_info[type == "SV" & Geno_RF_hap1 == 0 & nchar(ALT) == 1, subtype := "SV_DEL"]
  vcf_info[type == "SV" & Geno_RF_hap1 == 1 & nchar(REF) == 1, subtype := "SV_DEL"]
  vcf_info[type == "SV" & Geno_RF_hap1 == 1 & nchar(ALT) == 1, subtype := "SV_INS"]
  vcf_info[type == "SV" & is.na(subtype), subtype := "SV_complex"]
  vcf_info[type == "SV_complex", subtype := "SV_complex"]
  vcf_info[type == "SV" & is.na(subtype), subtype := "SV_complex"]
  vcf_info[type == "SV" & N_allele > 2, subtype := "Multiallelc_SV"]
  
  # Check counts by type and subtype
  vcf_info[, .N, by = type]
  vcf_info[, .N, by = subtype]
  
  # Save per-chromosome summary
  fwrite(vcf_info[order(POS)],file = paste0("Chrs_table/Sum_",gsub("VGP-prim#","",names(vcf_raw[x])),".txt"),sep = "\t")
  
  # Return processed VCF info
  vcf_info
}))

# Save the combined VCF information for all chromosomes
fwrite(vcf_info_38chrs,"VCF_INFO_allsites_2024Sep.txt",sep = "\t",quote = F,row.names = F,col.names = T)

#####################################
##### Summary of Variant Counts #####
#####################################
################
## Bi-allelic ##
################
tab_1 <- merge(vcf_info_38chrs[allele_count_West==2][,.(West=.N),by=subtype],
               vcf_info_38chrs[allele_count_East==2][,.(East=.N),by=subtype],all=T)

tab_2 <- vcf_info_38chrs[allele_count_All==2][,unique(.SD,by=c("CHROM","POS","subtype"))][,.(All=.N),by=subtype]
tab_biallelic <- merge(tab_2,tab_1,all=T)[,allelic_count:="Bi-allelic"]

###################
## Multi-allelic ##
###################
tab_3 <- merge(vcf_info_38chrs[allele_count_West>2][,.(West=.N),by=subtype],
               vcf_info_38chrs[allele_count_East>2][,.(East=.N),by=subtype],all=T)

tab_4 <- vcf_info_38chrs[allele_count_All>2][,unique(.SD,by=c("CHROM","POS","subtype"))][,.(All=.N),by=subtype]
tab_multiallelic <- merge(tab_4,tab_3,all=T)[,allelic_count:="Multi-allelic"]

#################################################
## Combine Bi-allelic and Multi-allelic Counts ##
#################################################
sum_pggb <- rbind(tab_biallelic, tab_multiallelic)[,c(5,1:4)][order(allelic_count,nchar(subtype))]
sum_pggb <- sum_pggb[subtype=="Multiallelc_SV",allelic_count:="Multi-allelic"][,.(All=sum(All),West=sum(West),East=sum(East)),by=.(allelic_count,subtype)]
sum_pggb[subtype %like% "SV",sum(All)]

sum_pggb_1 <- rbind(sum_pggb[allelic_count=="Bi-allelic"],
                    sum_pggb[allelic_count=="Multi-allelic"][subtype %like% "SNP"],
                    sum_pggb[allelic_count=="Multi-allelic"][subtype %like% "SV",.(allelic_count="Multi-allelic", subtype="SV_multiallelic",All=sum(All),West=sum(West),East=sum(East))],
                    sum_pggb[allelic_count=="Multi-allelic"][subtype %in% c("DEL","INS","INDEL_complex"),.(allelic_count="Multi-allelic", subtype="INDEL_multiallelic",All=sum(All),West=sum(West),East=sum(East))])

fwrite(sum_pggb_1,"Sum_PGGB_2024Sep19.txt",sep="\t",col.names = T)

# sum_pggb_1[subtype %like% "SV",sum(All)] # 887118 SVs recorded in PGGB VCF

