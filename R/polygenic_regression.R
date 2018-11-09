polygenicPTSD <- function(model,dframe,pop="eur"){

  bp_outcome <- model[1]
  prs_varname <- model[2]
  #antihtn <- model[3]
  modeltype <- model[4]
  
  #Get the columns related to PRS. They may be incorrectly ordered depending on R
  prs_vars <- grep(prs_varname, colnames(dframe),value=TRUE)


  #First create list of PRS variables 
  ptsd_vars <- colnames(dframe[grepl("ptsd", colnames(dframe))]) #selecting the four ptsd variable names
  
  
  #I will stratify the data by ancestry prior to this step so we can just rerun AAM later on.

#Loop runs through PRS variables within the correct subset specified above based on outcome
   for (ptsd in c(ptsd_vars)){ #Loop runs through PTSD variables
    for (test_type in c("*","+")){ #Loop runs through main effect/interaction
     if (test_type == "+"){ #Define model name for file title (Main Effects model)
       effect = "ME"
     }
     if (test_type == "*"){ #Define model name for file title (Interaction model)
       effect = "Int"
     }
     for (age_choice in c("1", "age")){ #Loop runs through age and age*age (still figuring out how to make good file name for this, currently it is age_1 and age_age)
      for (gender in c("all", "male", "female")){ #Loop subsets based on gender
       if (gender =="all"){
         dat_gen = dframe
         #But only add a sex covariate if there are really both sexes represented!
         if(dim(subset(dframe,sex==1))[1] > 0 & dim(subset(dframe,sex==0))[1] > 0)
         {
          sex= "+ sex"
         } 
          else { sex=" "}
       }
       if (gender =="male"){
         dat_gen = subset(dframe,sex==1)
         sex=" "
       }
       if (gender =="female" ){
         dat_gen = subset(dframe,sex==0)
         sex=" "
       } 
       
       #For the length of PRS, generate a list of object to store model results
       for (covar in c("base", "adjust")){ #Loop does original base with sex, pcs, age, then adds extra covars
        if (covar == "base"){
          covars = ""
        }
        #check if these covariates exist before perform this adjustment?
        if (covar == "adjust"){
          covars = "+ educ + trauma_count_lt"
        }
       
        #Make an age squared variable if the data choice was age squared...
        if(age_choice == "age")
        {
         dat_gen$agesq <- dat_gen$age^2
         agevar="+ agesq"
        } else {agevar = "+1"}
        
        if(dim(dat_gen)[1] > 0) #Only do a sex stratified analysis if the data frame has subjects to analyze 
        {
         modelname = paste(study, pop, gender, bp_outcome, prs_varname, ptsd, effect,  "age", age_choice,  "covar",covar, sep = ".") 
         mouts <- vector("list", length(prs_vars)+1)
         mouts_keller <- vector("list", length(prs_vars)+1)
         i=1
         for (prs in c(prs_vars)){ 
         
          modelformula = paste(bp_outcome, "~  P1 + P2 + P3 ", sex, "+ age",agevar, "+" ,covars, "+", ptsd, test_type, prs, sep = "") #'antihtn' is taken out for now, its not used in any model directly

          modeltype <- as.character(modeltype)
          print(modelformula)
          print(modelname)
          mouts[[i]] <- tryCatch(
                         summary(vglm(as.formula(modelformula), family = modeltype, data=dat_gen))@coef3 # , envir = .GlobalEnv
                        ,error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
          if (test_type == "*"){ #If this is an interaction test, also do the keller model that has cov x e (PTSD) interactions as well as cov x G (PRS)
          
            modelformula2 = paste(bp_outcome, "~  (P1 + P2 + P3 ", sex, "+ age",agevar, covars, ")*", ptsd, "+ (P1 + P2 + P3 ", sex, "+ age",agevar, covars, "+", ptsd,")*", prs, sep = "")
                      mouts_keller[[i]] <- tryCatch(
                         summary(vglm(as.formula(modelformula2), family = modeltype, data=dat_gen))@coef3 # , envir = .GlobalEnv
                        ,error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                        print(mouts_keller[[i]])
          }               
                        
          
          i = i+1
         } 
          
        }
        #Save model outputs as an R object
        save(mouts,file=paste(modelname,'r',sep='.'))
        save(mouts_keller,file=paste(modelname,"keller",'r',sep='.'))
        }
       }
      }
     }
    }
}
