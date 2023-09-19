#functions are written or modified by T. Oosting
#If you use these functions please cite Oosting et al., 2023
#packages
library(dplyr)
library(stringr)
library(ggplot2)
library(raster)
library(gdistance)
library(reshape2)
library(boot)
library(SNPRelate)
library(vcfR)

#Rfunctions
#####reading in vcf or gds file into R (SNPrelate)#####
snpgdsReadGDS <- function(vcf_file = NULL, gds_file = NULL){
  if(is.null(gds_file)){
    print("no gds file suplied, looking for gds file with similar name as vcf file")
    gds_file <- str_replace(vcf_file,".vcf.gz",".gds")
    if(file.exists(gds_file)){
      print("gds file found, loading now")
      gds <- SNPRelate::snpgdsOpen(gds_file,readonly = FALSE)
      return(gds)
    } else {
      print("no gds file, converting vcf to gds in same directory")
      SNPRelate::snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only")
      gds <- SNPRelate::snpgdsOpen(gds_file,readonly = FALSE)
      gdsfmt::add.gdsn(gds,"snp.id", paste0(read.gdsn(index.gdsn(gds, "snp.chromosome")),
                                            ":",
                                            read.gdsn(index.gdsn(gds, "snp.position"))) ,
                       replace = TRUE)
      return(gds)
    } 
  } else {
    print("gds file supplied")
    if(file.exists(gds_file)){
      print("gds file found, loading now")
      gds <- SNPRelate::snpgdsOpen(gds_file,readonly = FALSE)
      return(gds)
    } else {
      print("no gds file, converting vcf to gds if vcf is suplied")
      SNPRelate::snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only")
      gds <- SNPRelate::snpgdsOpen(gds_file,readonly = FALSE)
      gdsfmt::add.gdsn(gds,"snp.id", paste0(read.gdsn(index.gdsn(gds, "snp.chromosome")),
                                            ":",
                                            read.gdsn(index.gdsn(gds, "snp.position"))) ,
                       replace = TRUE)
      return(gds)
    }
  }
}

#####getting SNP information from gds file (SNPrelate)#####

#Obtain level of heterozygosity per SNP
snpgdsSNPHet <- function(gds       = NULL,
                         sample.id = NULL,
                         snp.id    = NULL,
						             margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                        sample.id = sample.id,
                                                                        snp.id = snp.id),
                                                                        margin,
                                                                        #verbose = FALSE,
                                                                        function(y){length(which(y==1))/length(which(!is.na(y)))})}

#Get reference allele count (AC) per SNP
snpgdsSNPACRef <- function(gds       = NULL,
                           sample.id = NULL,
                           snp.id    = NULL,
	   					             margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          margin,
                                                                          #verbose = FALSE,
                                                                          function(y){(length(which(y==2))*2)+length(which(y==1))})}

#Get frequency alternative allele
snpgdsSNPACAlt <- function(gds       = NULL,
                           sample.id = NULL,
                           snp.id    = NULL,
						               margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          margin,
                                                                          #verbose = FALSE,
                                                                          function(y){(length(which(y==0))*2)+length(which(y==1))})}

#Get alternative allele count (AC) per SNP
snpgdsSNPRef_freq <- function(gds    = NULL,
                           sample.id = NULL,
                           snp.id    = NULL,
						               margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          margin,
                                                                          #verbose = FALSE,
                                                                          function(y){((length(which(y==2))*2)+length(which(y==1)))/(length(which(!is.na(y)))*2)})}

#Get frequency reference allele
snpgdsSNPAlt_freq <- function(gds    = NULL,
                           sample.id = NULL,
                           snp.id    = NULL,
						               margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          margin,
                                                                          #verbose = FALSE,
                                                                          function(y){((length(which(y==0))*2)+length(which(y==1)))/(length(which(!is.na(y)))*2)})}

#Get frequency missing data
snpgdsSNPMissing <- function(gds       = NULL,
                             sample.id = NULL,
                             snp.id    = NULL,
						                 margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                            sample.id = sample.id,
                                                                            snp.id = snp.id),
                                                                            margin,
                                                                            #verbose = FALSE,
                                                                            function(y){length(which(is.na(y)))/length(y)})}




#SNP summary table from gds SNPRelate object
#function returns data frame with following statistics:
# Chromosome                             (CHR)
# Linkage group                          (LG) i.e. numerical variant to CHR
# snp.id                                 (LOC)
# snp.position                           (POS)
# Reference allele                       (REF)
# Alternative allele                     (ALT)
# fraction of missing data               (missing)
# minimum allele frequency               (MAF)
# heterozygosity                         (HET)
# Allele count reference allele          (AC_ref)
# Allele count alternative allele        (AC_alt)
# Frequency reference allele             (REF_freq)
# Frequency alternative allele           (ALT_freq)
# p-value for Hardy-Weinberg equilibrium (HWE_pval)
snpgdsSNPsum <- function(gds_snprelate = NULL,
                         sample.id     = NULL,
                         snp.id        = NULL,
                         extended      = FALSE){SNP_info <- data.frame(CHR = read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
                                                                       LG  = as.numeric( str_remove_all(read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")), "[:alpha:]|[:punct:]" )),
                                                                       LOC = paste0(read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
                                                                                    ":",
                                                                                    read.gdsn(index.gdsn(gds_snprelate, "snp.position"))),
                                                                       POS = read.gdsn(index.gdsn(gds_snprelate, "snp.position")),
                                                                       REF = stringr::str_extract(read.gdsn(index.gdsn(gds_snprelate, "snp.allele")),"^\\w"),
                                                                       ALT = stringr::str_extract(read.gdsn(index.gdsn(gds_snprelate, "snp.allele")),"\\w$"))

                              if(!is.null(snp.id)){SNP_info <- SNP_info[which(SNP_info$LOC %in% snp.id),]}
                              
                              if(extended){
                              SNP_info$MISSING  = SNPRelate::snpgdsSNPRateFreq(gds_snprelate, sample.id = sample.id, snp.id = snp.id)$MissingRate
                              SNP_info$MAF      = SNPRelate::snpgdsSNPRateFreq(gds_snprelate, sample.id = sample.id, snp.id = snp.id)$MinorFreq
                              SNP_info$HET      = snpgdsSNPHet(gds = gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$AC_ref   = snpgdsSNPACRef(gds = gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$AC_alt   = snpgdsSNPACAlt(gds = gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$REF_freq = snpgdsSNPRef_freq(gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$ALT_freq = snpgdsSNPAlt_freq(gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$MISS     = snpgdsSNPMissing(gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$HWE_pval = SNPRelate::snpgdsHWE(gds_snprelate, sample.id = sample.id, snp.id = snp.id)
                              }
                              return(SNP_info)
}

#same as function above but it requires a vector indicating populations
#this function will provide information about all SNPs per population
snpgdsSNPsum_byPOP <- function(gds       = NULL,
                               sample.id = NULL,
                               snp.id    = NULL,
                               pop.id    = NULL){ if(is.null(sample.id)){sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))}
  
                                                  tmp_df1 <- data.frame(matrix(ncol = 16, nrow = 0))
                                                  colnames(tmp_df1) <- c("CHR","LG","LOC","POS","REF","ALT","MISSING","MAF","HET","AC_ref","AC_alt","REF_freq","ALT_freq","MISS","HWE_pval","POP")
                                                  
                                                  pops <- unique(pop.id)
                                                  populations <- split(sample.id, pop.id)
                                                  
                                                  for (i in pops) {
                                                    tmp_df2 <- snpgdsSNPsum(gds, sample.id = populations[[i]], snp.id = snp.id, extended = TRUE)
                                                    tmp_df2$POP <- i
                                                    tmp_df1 <- rbind(tmp_df1,tmp_df2)}
                                                  return(tmp_df1)
}

#get either maf or mac per populations, the last column in the returned df indicates the number of populatios that have the alternative allele
#this analyses assumes that the alternative allele is the minar allele, i.e. always uses the ALT column for the vcf/gds file
#this this way maf/mac is always giving for the same allele per POP.
#function needs:
# gds    = gds file
# pop.id = population vector
# mac    = logical operator the return allele counts when true

snpgdsSNPmaf_byPOP <- function(gds       = NULL,
                               mac       = FALSE,
                               pop.id    = NULL){sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))
                               populations <- split(sample.id, pop.id)
                               
                               N <- lengths(populations)
                               P <- names(populations)
                               Npop <- length(N)
                               
                               tmp_df <- data.frame(LOC = read.gdsn(index.gdsn(gds, "snp.id")))
                               #get frequencies of the reference alleles
                               freq_REF <- snpgdsSNPList(gds)
                               #determine whether the REF or ALT allele has the highest freq
                               # 0 = REF allele has highest freq
                               # 1 = ALT allele has highest freq
                               # this sets the same reference allele for all populations regardless of the frequency
                               major_allele <- rep(0,nrow(freq_REF))
                               major_allele[which(freq_REF$afreq < 0.5)] <- 1
                               
                               major_freq <- abs( major_allele - freq_REF$afreq )
                               
                               for (i in P) {
                                 #if mac is TRUE alternative allele counts are given instead of frequency   
                                 if(!mac){
                                   X1 <- snpgdsSNPList(gds, sample.id = populations[[i]])
                                   X2 <- abs( major_allele - X1$afreq )
                                   X3 <- 1-X2
                                 } else {
                                   X1 <- length(populations[[i]])*2*major_allele
                                   X2 <- snpgdsSNPACAlt(gds = gds, sample.id = populations[[i]])
                                   X3 <- abs(X2-X1)
                                 }
                                 tmp_df <- cbind(tmp_df,X3)
                               }
                               
                               X <- apply(tmp_df[,-1],1,function(z){length( which(z > 0))})
                               tmp_df <- cbind(tmp_df,X)
                               colnames(tmp_df) <- c("LOC",paste0(P,"(",N,")"),"N")
                               
                               return(tmp_df)
}


#####getting sample information from gds file (SNPrelate)#####
#Obtain level of heterozygosity per sample
snpgdsINDHet <- function(gds       = NULL,
                         sample.id = NULL,
                         snp.id    = NULL){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          1,
                                                                          #verbose = FALSE,
                                                                          function(y){length(which(y==1))/length(which(!is.na(y)))})}



snpgdsINDMiss <- function(gds       = NULL,
                          sample.id = NULL,
                          snp.id    = NULL){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                           sample.id = sample.id,
                                                                           snp.id = snp.id),
                                                                           1,
                                                                           #verbose = FALSE,
                                                                           function(y){length(which(is.na(y)))/length(y)})}






#####estimating pairwise FST between sample locations#####

#estimate pwFSTs using SNPrelate
# 2 parametes required: gds & pop.id
# requires gds object generated via SNPrelate (not SeqArray!)
# vector for populations should reflect the order of samples in the gds object!!!
snpgdsPWFST <- function(gds = NULL, pop.id = NULL, snp.id = NULL, sample.id = NULL, out = NULL, reps = 1000, signf = 0.05,
                        exclude.monomorphic = TRUE, autosome.only = FALSE, FST.method = "W&C84", maf = 0.05, missing.rate = 0.95){
  #set vector of samples included FST estimates
  if(is.null(sample.id)){sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))}
  #extract population information 
  pops <- as.character(unique(pop.id))
  N_pops <- length(pops)
  populations <- split(sample.id,  pop.id)
  pop_combs <- as.data.frame(t(utils::combn(pops,2)))
  pop_combs_df <- as.data.frame(t(utils::combn(pops,2)))
  pop_combs_df <- pop_combs %>% mutate(merge = paste(V1,V2,sep = "."))
  N_comb <- nrow(pop_combs)
  pop_df <- data.frame(IND = sample.id,
                       POP = pop.id)
  #print population information
  print("starting pairwise FST estimation", quote = F)
  print("analyses includes:", quote = F)
  print(paste(length(sample.id),"individuals"), quote = F)
  print(paste(N_pops, "populations"), quote = F)
  print(paste(N_comb,"population combinations"), quote = F)
  #create empty data frame to print output
  pwFST_df <- data.frame(matrix(ncol = 16,nrow = 0))
  #loop over population combinations
  for(i in 1:N_comb){
    #get populations
    pop1       <- pop_combs_df$V1[i] 
    pop2       <- pop_combs_df$V2[i]
    indo_2pops <- pop_df %>% dplyr::filter(POP %in% c(pop1, pop2))
    #print information
    print("#######################################################", quote = F)
    print(paste("estimating",i,"out of",N_comb))
    print(paste("population1:",pop1))
    print(paste("population2:",pop2))
    #estimate PWFst 
    fst_2pop <- snpgdsFst(gdsobj         = gds,
                          sample.id      = indo_2pops$IND,
                          population     = factor(indo_2pops$POP,levels = c(pop1,pop2)),
                          method         = FST.method,
                          autosome.only  = autosome.only,
                          remove.monosnp = exclude.monomorphic,
                          maf            = maf,
                          missing.rate   = missing.rate)
    N_SNPs <- length(fst_2pop$FstSNP)
    # bootstrapping for CI over non-weighted FST estimates form all SNPs
    confid <- 1-signf
    # internal function to obtain the mean
    Bmean <- function(data,indices){
      d <- data[indices] # allows boot to select sample 
      return(mean(d))}
    results <- boot(data= fst_2pop$FstSNP, statistic=Bmean, R=reps)
    #create QQplots (only when out is specified)
    if(!is.null(out)){
      dir.create(paste0(dirname(out),"/QQplots"), recursive = TRUE, showWarnings = FALSE)
      png(paste0(dirname(out),"/QQplots/",basename(out),"_qqplot_",pop1,"_",pop2,"_",reps,"bootstrap.png"), width = 14, height = 7, units = "in", res = 300)
      plot(results)   
      dev.off()}
    
    #obtain pval
    Ttest <- t.test(results$t)
    pval  <- Ttest$p.value
    
    #obtain confidence intervals
    CIs <- boot::boot.ci(results, type=c("norm", "basic", "perc"), conf = confid)
    CI_vec <- c(CIs$normal[2],CIs$normal[3],CIs$basic[4],CIs$basic[5],CIs$percent[4],CIs$percent[5],pval)
    
    data <- as.character(c(pop1,pop2,fst_2pop$Fst,summary(fst_2pop$FstSNP),CI_vec,N_SNPs))
    pwFST_df <- pwFST_df %>% mutate_all(as.character)
    pwFST_df <- rbind(pwFST_df,data)
  }
  #assign column names
  colnames(pwFST_df) <- c("pop1", "pop2","FST_weighted", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max","CI95%_min_norm","CI95%_max_norm","CI95%_min_basic","CI95%_max_basic","CI95%_min_perc","CI95%_max_perc", "pval","N_snps")
  pwFST_df[,3:15] <- sapply(pwFST_df[,3:15],as.numeric)
  pwFST_df[,3:15] <- pwFST_df[,3:15] %>% rowwise() %>% round(digits = 7)
  pwFST_df <- pwFST_df %>% mutate(CI_range = paste(`CI95%_min_norm`,`CI95%_max_norm`, sep = " - "))
  #print pwFST table
  
  if(!is.null(out)){
    write_tsv(pwFST_df,paste0(out,"_pwFST.tsv"))
  }
  return(pwFST_df)
}

#convert output to matrix for plotting with pheatmap
snpgdsFSTdf2mat <- function(pwFST_df=NULL, order=NULL){
                            ### convert pwFST data frame to square matrix via distance matrix
                            ## set up storage matrix
                            # get names for row and columns
                            nameVals <- sort(unique(unlist(pwFST_df[1:2])))
                            # set factor for pops
                            if(!is.null(order)){
                              print("setting matrix to custom order")
                              nameVals <- nameVals[order]
                            }
                            # construct 0 matrix of correct dimensions with row and column names
                            myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
                            # fill in the matrix with matrix indexing on row and column names
                            myMat[as.matrix(pwFST_df[c("pop1", "pop2")])] <- pwFST_df[["FST_weighted"]]
                            myMat <- myMat + t(myMat)
                            return(myMat)
}


#####converting gds to bed file#####
snpgdsSNPbim <- function(gds_snprelate = NULL,
                         sample.id     = NULL,
                         snp.id        = NULL){SNP_info <- data.frame(chr = read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
                                                                      id  = paste0(read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
                                                                                   ":",
                                                                                   read.gdsn(index.gdsn(gds_snprelate, "snp.position"))),
                                                                      posg  = 0,
                                                                      pos = read.gdsn(index.gdsn(gds_snprelate, "snp.position")),
                                                                      ref = stringr::str_extract(read.gdsn(index.gdsn(gds_snprelate, "snp.allele")),"^\\w"),
                                                                      alt = stringr::str_extract(read.gdsn(index.gdsn(gds_snprelate, "snp.allele")),"\\w$"))
                         return(SNP_info)
}


snpgdsgds2bed <- function(gds = NULL, pop_df = NULL, out = NULL, exclude = NULL ){
  library(genio)
  
  ### get bim file ###
  bim <- snpgdsSNPbim(gds = gds)
  
  #check if CHR is numeric - change if not
  if(is.numeric(class(bim$chr))){
    print("CHR is in numeric format, no modifications made to CHR info field")
  } else {
    print("CHR is non-numeric, will replace CHR names with numeric sequence, original CHR name is still indicated in the variant identifier")
    chr_old <- unique(bim$chr)
    chr_new <- seq(length(chr_old))
    bim$chr <- chr_new[match(bim$chr, chr_old)]
  }
  
  ### get fam file ###
  if(!is.null(pop_df)){
    print("pop_file provided, generating fam file")
    print("looking for IND and POP fields")
    
    #vector if samples in the order they occr in the vcf/gds
    INDS_ordered <- read.gdsn(index.gdsn(gds, "sample.id"))
    
    if(!is.null(exclude)){
      print("sample to exclude were provided, removing now...")
      INDS_ordered <- setdiff(INDS_ordered,exclude)
      length(INDS_ordered)
    }
    
    #read pop file and filter the file for names present in vcf/gds
    pop_df <- pop_df %>% filter(IND %in% INDS_ordered) %>% dplyr::arrange(IND)
    #order samples in the order they occur in the vcf/gds
    pop_df[match(INDS_ordered, pop_df$IND),]
    #create df
    fam <- data.frame(fam   =  pop_df$POP,
                      id    =  pop_df$IND,
                      pat   =  0,
                      mat   =  0,
                      sex   =  0,
                      pheno = -9)
    
    #write_tsv(fam,file = glue("{out_ext}.fam"), col_names = F)  
  } else {
    INDS_ordered <- read.gdsn(index.gdsn(gds, "sample.id"))
    fam <- data.frame(fam   =  INDS_ordered,
                      id    =  INDS_ordered,
                      pat   =  0,
                      mat   =  0,
                      sex   =  0,
                      pheno = -9)
  }
  
  ### get bed file ###
  GT <- snpgdsGetGeno(gds)
  bed  <- t(GT)
  GT_cols <- which(read.gdsn(index.gdsn(gds, "sample.id")) %in% INDS_ordered)
  bed  <- bed[,GT_cols]
  
  ### write files ###
  #write new plink
  genio::write_plink(file = out, X = bed, fam = fam, bim = bim)
  
  
}

#collect chromosome information
#date frame requires following columns:
# Chromosome   as "CHR"
# SNP position as "POS"
chr_info <- function(SNP_df = NULL){chr_info <- SNP_df %>% group_by(CHR) %>% summarise(Length = max(POS)) %>% 
                                                arrange(factor(CHR, levels = unique(SNP_df$CHR)))         %>% 
                                                mutate(tot = cumsum(Length)-Length, 
                                                       LG  = c(1:length(unique(SNP_df$CHR))))
                                                return(chr_info)
}

#create data frame that contains information for where chromosome start and stop with in plot
#requires: 
# CHR
# LG
# BPcum
axisdf  <- function(SNP_df){ axisdf <- SNP_df %>% group_by(CHR) %>% summarize(LG     =  first(LG),
                                                                              center = (max(BPcum) + min(BPcum) ) / 2, 
                                                                              start  =  min(BPcum),
                                                                              stop   =  max(BPcum),
                                                                              ymin   =  -Inf,
                                                                              ymax   =  Inf) %>% arrange(LG)
                                                                              return(axisdf)
}
#create x-scale
x_scale <- function(axisdf){ x_scale <- scale_x_continuous(label=axisdf$CHR,breaks = axisdf$center,expand=c(0,0))
                            return(x_scale)
}


#create rectangles for background in ggplot
bg_rect <- function(axisdf){ bg_rect <- annotate("rect", xmin = axisdf$start,xmax = axisdf$stop, 
                                                         ymin = axisdf$ymin, ymax = axisdf$ymax, 
                                                         alpha = .5, 
                                                         fill = rep(c("darkgrey","white"),(nrow(axisdf)+1)/2)[1:nrow(axisdf)])
                                                         return(bg_rect)
}

#make manhattan plot
manhattan_outliers <- function(snp_df=NULL,outlier_df=NULL,out="manhattan_plot.png"){
  #internal functions:
  chr_info <- function(SNP_df = NULL){chr_info <- SNP_df %>% group_by(CHR) %>% summarise(Length = max(POS)) %>% 
    arrange(factor(CHR, levels = unique(SNP_df$CHR)))         %>% 
    mutate(tot = cumsum(Length)-Length, 
           LG  = c(1:length(unique(SNP_df$CHR))))
  return(chr_info)}
  axisdf  <- function(SNP_df){ axisdf <- SNP_df %>% group_by(CHR) %>% summarize(LG     =  first(LG),
                                                                                center = (max(BPcum) + min(BPcum) ) / 2, 
                                                                                start  =  min(BPcum),
                                                                                stop   =  max(BPcum),
                                                                                ymin   =  -Inf,
                                                                                ymax   =  Inf) %>% arrange(LG)
  return(axisdf)}
  x_scale <- function(axisdf){ x_scale <- scale_x_continuous(label=axisdf$CHR,breaks = axisdf$center,expand=c(0,0))
  return(x_scale)}
  bg_rect <- function(axisdf){ bg_rect <- annotate("rect", xmin = axisdf$start,xmax = axisdf$stop, 
                                                   ymin = axisdf$ymin, ymax = axisdf$ymax, 
                                                   alpha = .5, 
                                                   fill = rep(c("darkgrey","white"),(nrow(axisdf)+1)/2)[1:nrow(axisdf)])
  return(bg_rect)}
  #get chromosome information 
  chr_info <- chr_info(snp_df)    
  #add chr info to SNPinfo
  snp_df <- left_join(snp_df,chr_info[,c("CHR","LG","tot")]) %>% 
    arrange(LG,POS)                               %>% 
    mutate(BPcum = tot + POS)
  #ggplot objects
  axisdf <- axisdf(snp_df)
  bg_rect <- bg_rect(axisdf)
  x_scale <- x_scale(axisdf)
  
  #merge dataframes
  outlier_df$snp <- "nonsignf"
  outlier_df$snp[which(outlier_df$OutlierFlag)] <- "signf"
  outlier_df$snp[which(outlier_df$selected)] <- "select"
  outlier_df$snp <- factor(outlier_df$snp, levels = c("nonsignf","signf","select"))
  snp_df      <- left_join(snp_df,outlier_df) %>% filter(!is.na(FST))
  N_nignf     <- length(which(snp_df$OutlierFlag))
  N_selected  <- length(which(snp_df$selected))
  #make plot
  plot <- ggplot(snp_df,aes(x=BPcum,y=FST,color=snp)) + 
    bg_rect +
    geom_point() + 
    labs(title = "Manhattan plot",
         subtitle = paste(N_nignf,"significant SNPs,",N_selected,"selected"),
         x = "Chromosome",
         y = expression(italic(F)[ST])) +
    xlab("Chromosome")+
    coord_cartesian(ylim = c(min(snp_df$FST,na.rm =T)
                             ,max(snp_df$FST,na.rm =T)))+
    scale_color_manual(values = c("darkgrey","yellow","red")) +
    theme_bw() + 
    x_scale + 
    theme(panel.grid  = element_blank(),
          plot.margin = unit(c(0,0,0.0,0), units = "in"),
          axis.text.x = element_text(angle = 0, size = 12))
  ggsave(plot = plot, filename = out, width = 16, height = 4, dpi = 300, units = "in")
  #return(plot)
}

### distance
#function that estimates the ditance betrween locations taking into account boundaries between terrestrial and marine habitat
#required are:
# a shapefile of the terrestrial habitat of the area of interest
# extent of the area of interest in a vector as: c(min_long, max_long, min_lath, max_lath)
# GPS coordinates of samples or sample locations
# if interested in marine habitats set marine = TRUE! this flips the matrix
# if plots = TRUE, plot showing extent of area is shown
distance <- function(name       = NULL,
                     long       = NULL,
                     lath       = NULL,
                     shapefile  = NULL,
                     extent     = NULL,
                     gridsize   = 0.1,
                     directions = 8,
                     marine     = FALSE,
                     plot       = FALSE){#collect GPS information samples
                                         GPS_info <- data.frame(NAME = name,
                                                                LONG = long,
                                                                LATH = lath)  
                                          
                                         GPS_info <- GPS_info %>% group_by(NAME) %>% summarise(LONG = mean(LONG), LATH = mean(LATH))
                                          
                                         #create raster
                                         r <- raster(extent(extent), res=gridsize)
                                         rr <- rasterize(shapefile, r)
                                         
                                         #if marine is TRUE flip raster  
                                         if(marine){
                                            roce <- rr
                                            roce[is.na(roce[])] <- -99
                                            roce[roce[]>=1] <- NA
                                            roce[roce[]==-99] <- 1
                                          } else { roce <- rr }
                                         
                                         if(plot){
                                            plot(roce,ext=c(extent))
                                         }
                                         
                                         #create transition object
                                         troce <- transition(roce, mean, directions = directions) #less than 1 min
                                         #corrext for local distances
                                         troceC <- geoCorrection(troce, "c", scl = F)
                                         #put GPS coordinates in matrix
                                         pC_region <- as.matrix(GPS_info[c("LONG","LATH")])
                                         
                                         pC_region <- as.matrix(GPS_info[,-1])
                                         #calculate least-cost distance
                                         cosDist_reg <- costDistance(troceC, pC_region)
                                         #convert to Km and round to the nearest Km
                                         cosDist_reg <- cosDist_reg/1000
                                         cosDist_reg <- round(cosDist_reg,0)
                                         #put distances in square-matrix
                                         cosDist_reg_mat <- as.data.frame(as.matrix(cosDist_reg))
                                         colnames(cosDist_reg_mat) <- GPS_info$NAME
                                         rownames(cosDist_reg_mat) <- GPS_info$NAME
                                         cosDist_reg_mat <- as.matrix(cosDist_reg_mat)
                                         
                                         N_inf <-  length(which(is.infinite(cosDist_reg_mat[1,])))
                                         if(N_inf != 0){
                                           print(paste("warning,",N_inf,"individual(s) have GPS coordinates located on excluded grids"))
                                           print(paste("please look at individual(s):", row.names(cosDist_reg_mat)[  which(is.infinite(cosDist_reg_mat[1,]))]))
                                           print("this can be caused by wrong GPS coordinates, or coordinates too close to borders of of the grid edges")
                                           print("for now, slightly move GPS coorddinates so individuals are analysed correctly")
                                         }
                                         
                                         
                                         cosDist_reg_df <- reshape2::melt(as.matrix(cosDist_reg_mat),varnames = c("pop1","pop2"))
                                         colnames(cosDist_reg_df) <- c("pop1","pop2","distance")
                                         
                                         
                                         
                                         distance_list <- list()
                                         distance_list$matrix <- cosDist_reg_mat
                                         distance_list$dataframe <- cosDist_reg_df
                                         distance_list$raster <- roce
                                         return(distance_list)
}

### convert data frame to square matrix via distance matrix
df2squarematrix <- function(pop1 = NULL, pop2 = NULL, value = NULL){
  df <- data.frame(pop1 = pop1, pop2 = pop2, value = value)
  ## set up storage matrix
  # get names for row and columns
  nameVals <- sort(unique(unlist(df[1:2])))
  # construct 0 matrix of correct dimensions with row and column names
  myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
  # fill in the matrix with matrix indexing on row and column names
  myMat[as.matrix(df[c("pop1", "pop2")])] <- df[["value"]]
  #myMat <- myMat + t(myMat) #doubled all my values should remove after double checking
  return(myMat)
}


