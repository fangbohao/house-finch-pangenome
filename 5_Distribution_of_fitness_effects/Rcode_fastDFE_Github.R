library(fastdfe)
library(data.table)
library(parallel)
library(ggplot2)

fastdfe <- load_fastdfe() # load the fastdfe package
BaseInference <- fastdfe$BaseInference # import classes
Spectrum <- fastdfe$Spectrum # import Spectrum function from fastDFE

setwd("Pangenome/fastDFE") # working directory

#########
## SFS ##
#########
sfs_nue_all <- SFS[type %in% c("SNP","INDEL","SV") & genomic_region=="intergenic"] ## Neutral SFS
sfs_sel_all <- SFS[type %in% c("SNP","INDEL","SV") & genomic_region!="intergenic"]

###############
### fastDFE ###
###############
out_fastDFE <- lapply(1:nrow(sfs_sel_all),function(x){
  tab <- sfs_sel_all[x]
  pop <- tab$population

  ## Using neutral sites from SNP or INDEL
  rbindlist(lapply(c("SNP","INDEL","SV"),function(y){
    sfs_neutral <- sfs_nue_all[population==pop][type==y][,SFS]
    sfs_neutral <- as.numeric(strsplit(sfs_neutral,", ")[[1]])
    sfs_selection <- as.numeric(strsplit(tab$SFS,", ")[[1]])
    
    # create inference object
    sfs_neut <- Spectrum(sfs_neutral)
    sfs_sel <- Spectrum(sfs_selection)
    
    inf <- BaseInference(
      sfs_neut = sfs_neut,
      sfs_sel = sfs_sel,
      n_runs = 10,
      fixed_params = list(all=list(eps=0, S_b=10, p_b=0)),
      do_bootstrap = T
    )
    
    # run inference
    sfs_modelled <- BaseInference$run(inf)
    p <- BaseInference$plot_discretized(inf,intervals=c(-Inf, -100, -10, -1, 0))
    DFE_now <- data.table(p$data)[,group:=NULL][,setnames(.SD,c("S","fraction","ymin","ymax"))][,S:=c("100_inf","10_100","1_10","0_1")]
    DFE_now <- data.table(tab[,.(population,type,genomic_region)], neutral=y, DFE_now)
    DFE_now

  }))

})

tab_fastDFE <- rbindlist(out_fastDFE)

fwrite(tab_fastDFE,paste0("Out_fastDFE_",Sys.Date(),".txt"),sep = "\t")





