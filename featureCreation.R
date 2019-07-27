
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(caret)
library(DT)
library(MLmetrics)
library(ggplot2)
library(pROC)

#############################################################################
# initial data formatting
# for one variant the VAF was missing.  We assigned VAF to be exp(-200)

# the code refers to the following variables:
#    patientID :  patient ID
#    variantID :  ID for the variants
#    VAF :        variant allele frequency (mutation frequency)
#    impact_severity:   values are LOW, MEDIUM, HIGH
#    polyphen_score
#    sift_score
#    gerp_bp_score
#    cadd_scaled
#    clinPred
#    cosmic_ids
#    cancerCategory
#    type:        values are "indel" "SNP"
#    end :        end position
#    start :      start position.  Together with "end", calculate length
#    alt :        alternative allele
#    ref :        reference allele

# the code below refers to several datasets:
#    ind.data.subset : a dataset with all patients analyzesd, and multiple rows per patient 
		       if they had more than one variant (one row per variant).
			Contains patientID and variantID.
#    variant.data : annotation information for each variant.
			Contains variantID 
#    one.data: Just the rows for one patient, extracted from ind.data.subset

#########################################
# Values of -Inf were converted to zeros.
# Individuals with no variants at all were assigned 0 scores for most features
#############################################################################


############################################################################
# this is a generic function to summarize a feature across all variants
#identified in the same individual.  A function "functioncall" is provided
#specifying the desired transformation.  "vlabel" provides a name for the
#result.
############################################################################
perpersonsummary <- function(fn.ind.data, functioncall, vlabel) {
  # var.data global parameter
  patientids <- table(fn.ind.data$patientID)
  for (ii in (1:length(patientids))) {
    one.data <- fn.ind.data[fn.ind.data$patientID==names(patientids)[ii],]
    variants <- match(one.data$variantID, var.data$variant_id)
    sone.data <- data.frame(patientID = names(patientids)[ii],
                            v1 = functioncall(one.data,
                                              var.data[variants,]))
    if (ii==1) {summar.data <- sone.data}
    if (ii>1) {summar.data <- rbind(summar.data, sone.data)}
  }
  colnames(summar.data)[2] <- vlabel
  return(summar.data)
}
###########################################################################


##########################################################################
# Below are many specific feature creation functions.  They are
# each passed into perpersonsummary, and merged into an dataset to be
# used for analysis.
#######################################################################


fn.entropy <- function(one.data, variant.data) {
  return(sum(log(one.data$VAF)))  }
fn.entropy.high <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='HIGH']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(log(tmp)))}
}
fn.entropy.med <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='MED']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(log(tmp)))}
}
fn.entropy.low <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='LOW']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(log(tmp)))}
}
fn.1mentropy <- function(one.data, variant.data) {
  return(sum(log(1.00001-one.data$VAF)))  }
fn.mnvaf <- function(one.data, variant.data) {
  return(mean(one.data$VAF)) }
fn.maxvaf <- function(one.data, variant.data) {return(max(one.data$VAF))}
fn.max.polyphen <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$polyphen_score), na.rm=TRUE)
  if (is.na(tmp))  tmp <- 0
  return(tmp)}
fn.max.sift <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$sift_score), na.rm=TRUE)
  if (is.na(tmp)) tmp<-0
  return(tmp)}
fn.max.gerp_bp <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$gerp_bp_score), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}
fn.max.cadd <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$cadd_scaled), na.rm=TRUE)
  if (is.na(tmp)) tmp<-0
  return(tmp)}
fn.max.clinPred <- function(one.data, variant.data)  {
  tmp <- max(as.numeric(variant.data$clinPred), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0
  return(tmp) }
fn.mn.polyphen <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$polyphen_score)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }
fn.1me.polyphen <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$polyphen_score) *
               log(1.00001-one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.e.polyphen <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$polyphen_score) *
               log(one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.maxvaf.polyphen <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$polyphen_score) * one.data$VAF
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)}
fn.mn.sift <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$sift_score)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }
fn.1me.sift <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$sift_score) *
               log(1.00001-one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.e.sift <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$sift_score) *
               log(one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.maxvaf.sift <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$sift_score) * one.data$VAF, na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp)}
fn.mn.gerp_bp <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$gerp_bp_score)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }
fn.1me.gerp_bp <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$gerp_bp_score) *
               log(1.00001-one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.e.gerp_bp <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$gerp_bp_score) *
               log(one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.maxvaf.gerp_bp <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$gerp_bp_score) * one.data$VAF,na.rm=TRUE)
  return(tmp)}
fn.mn.cadd <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$cadd_scaled)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }
fn.1me.cadd <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$cadd_scaled) *
               log(1.00001-one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.e.cadd <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$cadd_scaled) *
               log(one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.maxvaf.cadd <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$cadd_scaled)*one.data$VAF, na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp)}
fn.mn.clinPred <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$clinPred)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }
fn.1me.clinPred <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$clinPred) *
               log(1.00001-one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.e.clinPred <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$clinPred) *
               log(one.data$VAF), na.rm=TRUE)
  if (is.na(tmp)) tmp <- 0; return(tmp)}
fn.maxvaf.clinPred <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$clinPred) * one.data$VAF, na.rm=TRUE)
  if (is.na(tmp)) tmp <-0
  return(tmp)}

fn.number <- function(one.data, variant.data)  {return(nrow(one.data))}
fn.number.high <- function(one.data, variant.data) {
  return(nrow(one.data[variant.data$impact_severity=='HIGH',])) }
fn.number.med <- function(one.data, variant.data) {
  return(nrow(one.data[variant.data$impact_severity=='MED',])) }
fn.number.low <- function(one.data, variant.data) {
  return(nrow(one.data[variant.data$impact_severity=='LOW',])) }
fn.diag <- function(one.data, variant.data) {return(one.data$tumourType[1])}
fn.tumorType <- function(one.data, variant.data) {return(one.data$tumourSubtype[1])}
fn.miss <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$type=='indel']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(length(tmp))}
}

# now summarize for each person using functions below
t1 <- perpersonsummary(ind.data.subset, fn.entropy, "entropy")
t.all <- t1
t2 <- perpersonsummary(ind.data.subset, fn.entropy.high, "entropyHIGH")
t.all <- merge(t1, t2, by.x=1, by.y=1)
t3 <- perpersonsummary(ind.data.subset, fn.entropy.med, "entropyMED")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
t4 <- perpersonsummary(ind.data.subset, fn.entropy.low, "entropyLOW")
t.all <- merge(t.all, t4, by.x=1, by.y=1)
t5 <- perpersonsummary(ind.data.subset, fn.1mentropy, "OME")
t.all <- merge(t.all, t5, by.x=1, by.y=1)
t6 <- perpersonsummary(ind.data.subset, fn.mnvaf, "MnVAF")
t.all <- merge(t.all, t6, by.x=1, by.y=1)
rm(t1); rm(t2); rm(t3); rm(t4); rm(t5); rm(t6);

t5 <- perpersonsummary(ind.data.subset, fn.number, "NmVnts")
t.all <- merge(t.all, t5, by.x=1, by.y=1)
t6 <- perpersonsummary(ind.data.subset, fn.number.high, "NmVntsHIGH")
t.all <- merge(t.all, t6, by.x=1, by.y=1)
t7 <- perpersonsummary(ind.data.subset, fn.number.med, "NmVntsMED")
t.all <- merge(t.all, t7, by.x=1, by.y=1)
t8 <- perpersonsummary(ind.data.subset, fn.number.low, "NmVntsLOW")
t.all <- merge(t.all, t8, by.x=1, by.y=1)
t9 <- perpersonsummary(ind.data.subset, fn.diag,"Diag")
t.all <- merge(t.all, t9, by.x=1, by.y=1)
#t10 <-perpersonsummary(ind.data.subset, fn.type,"Type")
#t.all <- merge(t.all, t10, by.x=1, by.y=1)
rm(t5); rm(t6); rm(t7); rm(t8);  rm(t9); rm(t10)

t1 <- perpersonsummary(ind.data.subset, fn.maxvaf, "maxvaf")
t.all <- merge(t.all, t1, by.x=1, by.y=1)
t2 <- perpersonsummary(ind.data.subset, fn.max.polyphen, "max.polyphen")
t.all <- merge(t.all, t2, by.x=1, by.y=1)
t3 <- perpersonsummary(ind.data.subset, fn.max.sift, "max.sift")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
t4 <- perpersonsummary(ind.data.subset, fn.max.gerp_bp, "max.gerp_bp")
t.all <- merge(t.all, t4, by.x=1, by.y=1)
t5 <- perpersonsummary(ind.data.subset, fn.max.cadd, "max.cadd")
t.all <- merge(t.all, t5, by.x=1, by.y=1)
#t6 <- perpersonsummary(ind.data.subset, fn.max.gerp_pval, "max.gerp_pval")
#t.all <- merge(t.all, t6, by.x=1, by.y=1)
t7 <- perpersonsummary(ind.data.subset, fn.max.clinPred, "max.clinPred")
t.all <- merge(t.all, t7, by.x=1, by.y=1)
rm(t1);rm(t2); rm(t3); rm(t4); rm(t5); rm(t7); #rm(t6); 

t1 <- perpersonsummary(ind.data.subset, fn.maxvaf.polyphen, "maxvaf.polyphen")
t.all <- merge(t.all, t1, by.x=1, by.y=1)
t2 <- perpersonsummary(ind.data.subset, fn.1me.polyphen,"OME.polyphen")
t.all <- merge(t.all, t2, by.x=1, by.y=1)
t3 <- perpersonsummary(ind.data.subset, fn.mn.polyphen,"mnvaf.polyphen")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
t4 <- perpersonsummary(ind.data.subset, fn.e.polyphen,"E.polyphen")
t.all <- merge(t.all, t4, by.x=1, by.y=1)
rm(t1); rm(t2); rm(t3); rm(t4)

t1 <- perpersonsummary(ind.data.subset, fn.maxvaf.sift, "maxvaf.sift")
t.all <- merge(t.all, t1, by.x=1, by.y=1)
t2 <- perpersonsummary(ind.data.subset, fn.1me.sift,"OME.sift")
t.all <- merge(t.all, t2, by.x=1, by.y=1)
t3 <- perpersonsummary(ind.data.subset, fn.mn.sift,"mnvaf.sift")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
t4 <- perpersonsummary(ind.data.subset, fn.e.sift,"E.sift")
t.all <- merge(t.all, t4, by.x=1, by.y=1)
rm(t1); rm(t2); rm(t3); rm(t4)

t1 <- perpersonsummary(ind.data.subset, fn.maxvaf.gerp_bp, "maxvaf.gerp_bp")
t.all <- merge(t.all, t1, by.x=1, by.y=1)
t2 <- perpersonsummary(ind.data.subset, fn.1me.gerp_bp,"OME.gerp_bp")
t.all <- merge(t.all, t2, by.x=1, by.y=1)
t3 <- perpersonsummary(ind.data.subset, fn.mn.gerp_bp,"mnvaf.gerp_bp")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
t4 <- perpersonsummary(ind.data.subset, fn.e.gerp_bp,"E.gerp_bp")
t.all <- merge(t.all, t4, by.x=1, by.y=1)
rm(t1); rm(t2); rm(t3); rm(t4)

t1 <- perpersonsummary(ind.data.subset, fn.maxvaf.cadd, "maxvaf.cadd")
t.all <- merge(t.all, t1, by.x=1, by.y=1)
t2 <- perpersonsummary(ind.data.subset, fn.1me.cadd,"OME.cadd")
t.all <- merge(t.all, t2, by.x=1, by.y=1)
t3 <- perpersonsummary(ind.data.subset, fn.mn.cadd,"mnvaf.cadd")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
t4 <- perpersonsummary(ind.data.subset, fn.e.cadd,"E.cadd")
t.all <- merge(t.all, t4, by.x=1, by.y=1)
rm(t1); rm(t2); rm(t3); rm(t4)

t1 <- perpersonsummary(ind.data.subset, fn.maxvaf.clinPred, "maxvaf.clinPred")
t.all <- merge(t.all, t1, by.x=1, by.y=1)
t2 <- perpersonsummary(ind.data.subset, fn.1me.clinPred,"OME.clinPred")
t.all <- merge(t.all, t2, by.x=1, by.y=1)
t3 <- perpersonsummary(ind.data.subset, fn.mn.clinPred,"mnvaf.clinPred")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
t4 <- perpersonsummary(ind.data.subset, fn.e.clinPred,"E.clinPred")
t.all <- merge(t.all, t4, by.x=1, by.y=1)
t5 <- perpersonsummary(ind.data.subset, fn.tumorType,"TumorSubtype")
t.all <- merge(t.all, t5, by.x=1, by.y=1)

rm(t1); rm(t2); rm(t3); rm(t4); rm(t5)


fn.cosm.id <- function(one.data, variant.data) {
  return( as.numeric( nrow(one.data) - nrow(one.data[variant.data$cosmic_ids=='None',])>0) ) }

t7 <- perpersonsummary(ind.data.subset, fn.cosm.id,"cosm.id")
t.all <- merge(t.all, t7, by.x=1, by.y=1)


####################################
#A count #indels
count.indel <- function(one.data, variant.data) {
  f<-sum(variant.data$type=="indel") 
  if(is.na(f))f<-0
  return(f)
}

t1 <- perpersonsummary(ind.data.subset, count.indel, "count.indel")

t.all <- merge(t.all, t1, by.x=1, by.y=1)
print(dim(t.all))
####################################
#B count #SNPs
count.snp <- function(one.data, variant.data) {
  f<-sum(variant.data$type=="snp") 
  if(is.na(f))f<-0
  return(f)
}

t2 <- perpersonsummary(ind.data.subset, count.snp, "count.snp")
t.all <- merge(t.all, t2, by.x=1, by.y=1)
print(dim(t.all))

####################################
#C count total #alterrations
count.all <- function(one.data, variant.data) {
  f<-sum(variant.data$type %in% c("snp" ,"indel" ) ) 
  if(is.na(f))f<-0
  return(f)
}

t3 <- perpersonsummary(ind.data.subset, count.all, "count.all")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
print(dim(t.all))

####################################
#VAF weighted A
VAF.count.indel <- function(one.data, variant.data) {
  f<-sum( (variant.data$type=="indel") * one.data$VAF) 
  if(is.na(f))f<-0
  return(f)
}

t1 <- perpersonsummary(ind.data.subset, VAF.count.indel, "VAF.count.indel")
t.all <- merge(t.all, t1, by.x=1, by.y=1)
print(dim(t.all))

####################################
#VAF weighted B
VAF.count.snp <- function(one.data, variant.data) {
  f<-sum( (variant.data$type=="snp")* one.data$VAF )
  if(is.na(f))f<-0
  return((f))
}

t1 <- perpersonsummary(ind.data.subset, VAF.count.snp, "VAF.count.snp")
t.all <- merge(t.all, t1, by.x=1, by.y=1)
print(dim(t.all))

####################################
#VAF weighted c
VAF.count.all <- function(one.data, variant.data) {
  f<-sum( (variant.data$type %in% c("snp" ,"indel" ) ) * one.data$VAF )
  if(is.na(f))f<-0
  return(f)
}

t3 <- perpersonsummary(ind.data.subset, VAF.count.all, "VAF.count.all")
t.all <- merge(t.all, t3, by.x=1, by.y=1)
print(dim(t.all))

####################################
#maximum (gerp * vaf *A)

max.vaf.gerp.indel <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type=="indel") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, max.vaf.gerp.indel, "max.vaf.gerp.indel")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)
print(dim(t.all))
####################################
#maximum (gerp * vaf *c)

max.vaf.gerp.all <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type%in%c("indel","snp") ) 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, max.vaf.gerp.all, "max.vaf.gerp.all")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)
print(dim(t.all))

####################################
#mean (gerp * vaf *A)

mn.vaf.gerp.indel <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type=="indel") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, mn.vaf.gerp.indel, "mn.vaf.gerp.indel")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)
print(dim(t.all))
####################################
#mean (gerp * vaf *c)

mn.vaf.gerp.all <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type%in%c("indel","snp") ) 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, mn.vaf.gerp.all, "mn.vaf.gerp.all")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)
print(dim(t.all))
####################################

#maximum (polyphen score * vaf * B)
max.vaf.poly.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$polyphen_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, max.vaf.poly.snp, "max.vaf.poly.snp")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)
print(dim(t.all))

####################################
max.vaf.sift.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$sift_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, max.vaf.sift.snp, "max.vaf.sift.snp")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)

print(dim(t.all))

#I have minus infinity here.

####################################
max.vaf.gerp.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, max.vaf.gerp.snp, "max.vaf.gerp.snp")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)

####################################
max.vaf.cadd.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$cadd_scaled) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, max.vaf.cadd.snp, "max.vaf.cadd.snp")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)


####################################
max.vaf.clin.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$clinPred) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, max.vaf.clin.snp, "max.vaf.clin.snp")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)


####################################

mean.vaf.poly.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$polyphen_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, mean.vaf.poly.snp, "mean.vaf.poly.snp")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)


#I have minus infinity here.
####################################
mean.vaf.sift.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$sift_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, mean.vaf.sift.snp, "mean.vaf.sift.snp")
#I have minus infinity here.
t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)

####################################
mean.vaf.gerp.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, mean.vaf.gerp.snp, "mean.vaf.gerp.snp")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)


####################################
mean.vaf.cadd.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$cadd_scaled) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, mean.vaf.cadd.snp, "mean.vaf.cadd.snp")

t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)


####################################
mean.vaf.clin.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$clinPred) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t3 <- perpersonsummary(ind.data.subset, mean.vaf.clin.snp, "mean.vaf.clin.snp")
t3[which(t3[,2]=="-Inf"),2] <- 0
t.all <- merge(t.all, t3, by.x=1, by.y=1)


####################################

no.variant <- function(one.data, variant.data) {
  tmp <- as.numeric(sum(!is.na(one.data$variantID)) )
  return(tmp)  
}

#t3 <- perpersonsummary(ind.data, no.variant, "no.variant")

t3 <- perpersonsummary(ind.data.subset, no.variant, "no.variant")
t.all <- merge(t.all, t3, by.x=1, by.y=1)

##
#now add the features of entropy

cor.entropy <-  function(one.data, variant.data) {
  return(sum(one.data$VAF*log(one.data$VAF)))  }
t4 <- perpersonsummary(ind.data.subset, cor.entropy, "cor.entropy")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


cor.entropy.high <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='HIGH']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(tmp*log(tmp)))}
  if(is.na(tmp)==1){return(0)}
}
t4 <- perpersonsummary(ind.data.subset, cor.entropy.high, "cor.entropy.high")
t.all <- merge(t.all, t4, by.x=1, by.y=1)

cor.entropy.med <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='MED']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(tmp*log(tmp)))}
}
t4 <- perpersonsummary(ind.data.subset, cor.entropy.med, "cor.entropy.med")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


cor.entropy.low <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='LOW']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(tmp*log(tmp)))}
}
t4 <- perpersonsummary(ind.data.subset, cor.entropy.low, "cor.entropy.low")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


cor.1mentropy <- function(one.data, variant.data) {
  return(sum((1.00001-one.data$VAF)*log(1.00001-one.data$VAF)))  }
t4 <- perpersonsummary(ind.data.subset, cor.1mentropy, "cor.1mentropy")
t.all <- merge(t.all, t4, by.x=1, by.y=1)

####################################################################################
#the characteristics of size of base

max.size.base <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$end - variant.data$start), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}

t4 <- perpersonsummary(ind.data.subset, max.size.base, "max.size.base")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


sum.size.base <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$end - variant.data$start), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}

t4 <- perpersonsummary(ind.data.subset, sum.size.base, "sum.size.base")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


max.size.bp.change <- function(one.data, variant.data) {
  tmp <- max(-( nchar(variant.data$alt) - nchar(variant.data$ref)), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}

t4 <- perpersonsummary(ind.data.subset, max.size.bp.change, "max.size.bp.change")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


sum.size.bp.change <- function(one.data, variant.data) {
  tmp <- sum(-( nchar(variant.data$alt) - nchar(variant.data$ref)), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}

t4 <- perpersonsummary(ind.data.subset, sum.size.bp.change, "sum.size.bp.change")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


####################################################################################
#VAF weighted size of change
mean.VAF.size.change  <- function(one.data, variant.data) {
  tmp <- as.numeric( ( nchar(variant.data$alt) - nchar(variant.data$ref)) ) * 
    one.data$VAF 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

t4 <- perpersonsummary(ind.data.subset, mean.VAF.size.change, "mean.VAF.size.change")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


mean.VAF.size.change.snp  <- function(one.data, variant.data) {
  tmp <- as.numeric( ( nchar(variant.data$alt) - nchar(variant.data$ref)) ) * 
    (one.data$VAF) * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

fn.category <- function(one.data, variant.data) {return(one.data$CancerCategory[1])}

fn.sum.log.vaf <- function(one.data, variant.data) {
  return(sum( log(one.data$VAF)) ) }

fn.sum.log.om.vaf <- function(one.data, variant.data) {
  return(sum( log(1-one.data$VAF)) ) }

t4 <- perpersonsummary(ind.data.subset, fn.sum.log.vaf, "sum.log.vaf")
t.all <- merge(t.all, t4, by.x=1, by.y=1)


t4 <- perpersonsummary(ind.data.subset, fn.sum.log.om.vaf, "sum.log.om.vaf")
t.all <- merge(t.all, t4, by.x=1, by.y=1)




