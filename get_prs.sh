#Download and install PRSice

https://choishingwan.github.io/PRSice/

#Get the variant rs-id names and stuff from UKBB
 wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

#Split the variant file into just the rs-ids, position, chr, ref and coded variants
  zcat variants.tsv.bgz | awk '{print $1,$6,$2,$3,$4,$5}' > variants_short.tsv

  #Get SBP
  wget https://www.dropbox.com/s/cs6sm1a8wnle966/4080_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4080_raw.gwas.imputed_v3.both_sexes.tsv.bgz
  #Get DBP
  wget https://www.dropbox.com/s/xq7i6wdvw4ov6c8/4079_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4079_raw.gwas.imputed_v3.both_sexes.tsv.bgz
  #Get HT
  wget https://www.dropbox.com/s/m8qlfp0cjnn4ka7/20002_1065.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20002_1065.gwas.imputed_v3.both_sexes.tsv.bgz
#Only take bi-allelic non-indel variants with MAF > 5%
#Note: Alt allele is the A1 allele. In other words, caution: the minor allele column in the file is not always the coded allele!!
join variants_short.tsv <(zcat 4080_raw.gwas.imputed_v3.both_sexes.tsv.bgz)  | grep -v NaN    | awk '{ if (NR == 1 || (length($5) == 1 && length($6) == 1)) print $2,$3,$4,$6,$5,$8,$10,$13,$14,$16}' | sed 's/^X/23/g' | awk '{if (NR == 1 || ($6 > 0.05 && $6 <0.95)) print}' > SBP.tsv
join variants_short.tsv <(zcat 4079_raw.gwas.imputed_v3.both_sexes.tsv.bgz)  | grep -v NaN    | awk '{ if (NR == 1 || (length($5) == 1 && length($6) == 1)) print $2,$3,$4,$6,$5,$8,$10,$13,$14,$16}' | sed 's/^X/23/g' | awk '{if (NR == 1 || ($6 > 0.05 && $6 <0.95)) print}' > DBP.tsv
#Note, since this is a case/control trait, shift some columns over by 1
join variants_short.tsv <(zcat 20002_1065.gwas.imputed_v3.both_sexes.tsv.bgz)  | grep -v NaN  | awk '{ if (NR == 1 || (length($5) == 1 && length($6) == 1)) print $2,$3,$4,$6,$5,$8,$11,$14,$15,$17}' | sed 's/^X/23/g' | awk '{if (NR == 1 || ($6 > 0.05 && $6 <0.95)) print}' > HT.tsv


#fix kripke code to note this..


#This is going to run PRSIce. Mainly, it's to get the scores themselves. But also, if you supply valid phenotyping, it will give a basic association analysis. 
#This is useful for making sure that the scores are predictive in the correct direction in the data.

study=mrsc
phenofile=mrs1_bp.pheno.txt #write name of phenotype file. Should be formatted like a PLINK phenotype file with FID IID then phenotype columns
phenoname=SBP #Write column name of phenotype to be analyzed

mkdir output # I put the outputs in a directory called output. 

Rscript PRSice.R \
--prsice PRSice_linux \
--base SBP.tsv \
--target genotypes/"$study"_bgn_eur_chr"#" \
--thread 16 \
--full \
--nonfounders \
--snp rsid --stat beta --se se  --A1 alt --A2 ref  --pvalue pval \
--pheno-file $phenofile --pheno-col $phenoname \
--binary-target F \
--upper 1 \
--interval 0.01 \
--all-score \
--out output/"$study"_sbpr

#There is an "$study"_sbpr.all.score file that has each subject's PRS at each risk score threshold.
#We'll change the column names from just p-value thresholds to also include the type of PRS.
#This will let the data read into the R script more easily..

head -n1  output/"$study"_sbpr.all.score | sed 's/\<[0-9]\>/sbp_prs_&/g' > sbp_header.txt

cat sbp_header.txt  <(tail -n+2 output/"$study"_sbpr.all.score) > output/"$study"_sbpr.all.score_prs


###Now just do the same thing for the other two blood pressure traits...

#DBP
phenoname=DBP 

Rscript PRSice.R \
--prsice PRSice_linux \
--base DBP.tsv \
--target genotypes/"$study"_bgn_eur_chr"#" \
--thread 16 \
--full \
--nonfounders \
--snp rsid --stat beta --se se  --A1 alt --A2 ref  --pvalue pval \
--pheno-file $phenofile --pheno-col $phenoname \
--binary-target F \
--upper 1 \
--interval 0.01 \
--all-score \
--out output/"$study"_dbpr

head -n1  output/"$study"_dbpr.all.score | sed 's/\<[0-9]\>/dbp_prs_&/g' > dbp_header.txt
cat dbp_header.txt  <(tail -n+2 output/"$study"_dbpr.all.score) > output/"$study"_dbpr.all.score_prs


#Hypertension
phenoname=HT

Rscript PRSice.R \
--prsice PRSice_linux \
--base HT.tsv \
--target genotypes/"$study"_bgn_eur_chr"#" \
--thread 16 \
--full \
--nonfounders \
--snp rsid --stat beta --se se  --A1 alt --A2 ref  --pvalue pval \
--pheno-file $phenofile --pheno-col $phenoname \
--binary-target T \
--upper 1 \
--interval 0.01 \
--all-score \
--out output/"$study"_htnr

head -n1  output/"$study"_htnr.all.score | sed 's/\<[0-9]\>/htn_prs_&/g' > htn_header.txt
cat htn_header.txt  <(tail -n+2 output/"$study"_htnr.all.score) > output/"$study"_htnr.all.score_prs

##That's it. You now have PRS that you can use in the pipeline

