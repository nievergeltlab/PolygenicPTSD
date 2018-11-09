
<!-- README.md is generated from README.Rmd. Please edit that file -->
polygenicPTSD
=============

The goal of polygenicPTSD is to make an easy-to-use package for multi-center sites in the PGC-PTSD Physical Health Workgroup to conduct analysis on the interaction between PTSD and BP Polygenic Risk Score (PRS) in predicting various BP outcomes. It includes functions to ensure data has appropriate variables required for analysis, creates new variables which categorizes continuous BP into categorical hypertension outcomes based on old and new AHA criteria, stratifies data by ancestry, runs linear regressions on continuous BP outcomes and generalized linear models on categorical hypertension outcomes, and more.

Installation
------------

You can install polygenicPTSD from github with:

``` r
# install.packages("devtools")
devtools::install_github("nievergeltlab/PolygenicPTSD")
```

Example
-------

This is a basic example which shows you how to run the analysis pipeline on simulated data:

``` r
##To user: Replace STUDYPI below with last name of Principle Investigator (For example "Sumner")
study <- "StudyPI"


#To user: if you need to install R packages, this is how:
#install.packages("VGAM", "dplyr") 

#Load polygenicPTSD package and packages it depends on 
 library(polygenicPTSD)
 library(dplyr)
 library(VGAM)

#Load regression model parameters
 data("model_types") 

## STEP 1a Load and merge
## Load and merge data

# merged_example is an example dataset. 
# Your data structure should resemble this:
# data("merged_example")
# Variables should be named as they are in this dataset!

#To user: Load your data as a .csv file
#To user: Write in data file name where it says study_data.csv

#Note: Do not put PRS scores in your file manually! 
#They will be merged into the data in a subsequent step

 d1 <-  read.csv('mrs_bp_sumner.csv',header=T,stringsAsFactors=F,na.strings=c(NA,"#N/A","-9"))

#Merge polygenic risk scores with data itself
 sbp_prs <- read.table('output/mrsc_sbpr.all.score_prs',header=T,stringsAsFactors=F)
 dbp_prs <- read.table('output/mrsc_dbpr.all.score_prs',header=T,stringsAsFactors=F)
 htn_prs <- read.table('output/mrsc_htnr.all.score_prs',header=T,stringsAsFactors=F)

 d1a <- merge(d1,sbp_prs,by=c("FID","IID"),all=TRUE)
 d1b <- merge(d1a,dbp_prs,by=c("FID","IID"),all=TRUE)
 dat <- merge(d1b,htn_prs,by=c("FID","IID"),all=TRUE)
 

## STEP 1b Categorization
## Create dataframe with 7 new categorical variables for hypertension 
## classification, continuous BP adjustment for antihypertension use, 
## and ancestry (necessary for analysis)

classed <- htn_outcome_classify(dat)
table(classed$htn_aha_new)

## Create ancestry stratification indices

index_eur = grepl(0, classed$ancestry)
index_afr = grepl(1, classed$ancestry)


## STEP 2 Summary statistics

#Save summary statistics, stratified by ancestry
variable.summary.stratified(classed)

## Fisher's test on old vs new AHA guideline classifications stratified by ancestry
x = paste(study, "ancestry_strat_test", sep = "_") #creates study-specific name using input from earlier in script 
assign(x, ancestry.strat.test(classed), envir = .GlobalEnv) #uses custom function to make contingency tables of AHA old vs new, stratified by ancestry, and tests those tables with Fisher's exact
save(list = x, file=paste(x, ".RData", sep="_")) #Saves test outputs and contingency tables as an R list

## Fisher's test on old vs new AHA guideline classifications 
## stratified by ancestry

## Below functions all save output in working directory.
## We recommend you create a folder "results" and change your working directory to that location.

setwd("results/")


######################################################################
#Step 3: Test outcomes using functions from polygenicPTSD
######################################################################
#Each of the below uses functions from the polygenicPTSD package to run the required models for the continous and categorical outcomes.

#For now, PRS are based only on European ancestry data
#So we'll subset the data to just that. Updates will come later for AAMs!

classed_eur <- subset(classed,bestpop == "eur")
#Regression analysis 

#User:Please specify a comma delimited list of additional covariates, with quotations around the names. 
#By default example, this will be education and lifetime trauma count
#If you dont have these, just write: covars=NA
#Be sure you spell the names correctly! Otherwise all covariates will be dropped from the model

covars=c("educ","trauma_count_lt")




for (i in c(1:8))
{
 polygenicPTSD(model_types[i,],dframe=classed_eur,pop="eur",covlist=covars)
}


## For any questions on documentation refer to full pipeline document 
## and/or type: ?function to see help documentation

?merged_example
?polygenicPTSD_logit
```
