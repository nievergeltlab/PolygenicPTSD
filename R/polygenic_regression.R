polygenicPTSD <- function(model,dframe,pop="eur"){

  bp_outcome <- model[1]
  prs_varname <- model[2]
  antihtn <- model[3]
  modeltype <- model[4]
  
  #Get the columns related to PRS. They may be incorrectly ordered depending on R
  prs_vars <- grep(prs_varname, colnames(dframe),value=TRUE)


  #First create list of PRS variables 
  ptsd_vars <- colnames(dframe[grepl("ptsd", colnames(dframe))]) #selecting the four ptsd variable names
  
  
  #I will stratify the data by ancestry prior to this step so we can just rerun AAM later on.

#Loop runs through PRS variables within the correct subset specified above based on outcome
   for (ptsd in c(ptsd_vars)){ #Loop runs through PTSD variables
    for (test_type in c("+", "*")){ #Loop runs through main effect/interaction
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
       if (gender =="female"){
         dat_gen = subset(dframe,sex==0)
         sex=" "
       } #Will fail if there are no individuals of a given sex unless using trycatch. Maybe use the subset command when runnign the lm, and use that under try??
       
       #For the length of PRS, generate a list of object to store model results
       for (covar in c("base", "adjust")){ #Loop does original base with sex, pcs, age, then adds extra covars
        if (covar == "base"){
          covars = ""
        }
        #check if these covariates exist before perform this adjustment?
        if (covar == "adjust"){
          covars = "+ educ + trauma_count_lt"
        }

        modelname = paste(study, pop, gender, bp_outcome, prs_varname, ptsd, effect,  "age", age_choice,  "covar",covar, sep = ".") 
        mouts <- vector("list", length(prs_vars))
        i=1
        for (prs in c(prs_vars)){ 
         modelformula = paste(bp_outcome, "~  P1 + P2 + P3 ", sex, "+ age*",age_choice, "+","+" ,covars, "+", ptsd, test_type, prs, sep = "") #'antihtn' is taken out for now, its not used in any model directly
         print(modelformula)
         modeltype <- as.character(modeltype)
         print(modelname)
         mouts[[i]] <- tryCatch(
                        vglm(as.formula(modelformula), family = modeltype, data=dat_gen) # , envir = .GlobalEnv
                       ,error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
         
         i = i+1
          
        }
         #print(modelname)
        #Save model outputs as an R object
        save(mouts,file=paste(modelname,'r',sep='.'))
        }
       }
      }
     }
    }
}
