library(metafor)

#Loop over all analysis to get lowest p-value association for a given ananylsis
#Load a model output'

 for (pop in c("eur","aam")) {
  for (gender in c("all","male","female")) {
   for (bp_outcome in c("htn_aha_new","htn_aha_old","htn_aha_new_bi","htn_aha_old_bi","DBP_meas_adj_5","DBP_meas_adj_10","SBP_meas_adj_10","SBP_meas_adj_05")) {
   
    if(bp_outcome %in% c("htn_aha_new","htn_aha_old","htn_aha_new_bi","htn_aha_old_bi")){
     prs_varname = "htn_prs"
    } else if( bp_outcome %in% c("DBP_meas_adj_5","DBP_meas_adj_10")) {
      prs_varname = "dbp_prs"
     } else if( bp_outcome %in% c("SBP_meas_adj_10","SBP_meas_adj_05")) {
       prs_varname = "sbp_prs"
     } else { prs_varname= "ERROR"; print ("no handling for bp_outcome - prs varname!!!")}
      for (ptsd in c("ptsd_diag_cur","ptsd_diag_lt","ptsd_sx_cur","ptsd_sx_lt")){ #Loop runs through PTSD variables. May have to harmonize names, i.e. just copy the files with new names 
       for (effect in c("ME","Int")) {
        for (age in c("1","age")) {
         for (covar in c("base","adjust")) {
         
          pop="eur"
          gender="all"
          bp_outcome="SBP_meas_adj_10"
          prs_varname="sbp_prs"
          ptsd="ptsd_sx_cur"
          effect="Int"
          age="age"
          covar="base"

          #If keller model, add a .keller in name.. or just do keller only.
          if(effect == "ME") {
           modelname = paste(pop, gender, bp_outcome, prs_varname, ptsd, effect,  "age", age,  "covar",covar, "r",sep = ".") 
          } else if(effect == "Int") {
           modelname = paste(pop, gender, bp_outcome, prs_varname, ptsd, effect,  "age", age,  "covar",covar,"keller", "r", sep = ".") 
          }
          
          #Have to handle that some analysis dont exist... read from a file parameter list maybe??

          if(effect == "ME"){
           load(paste('results/gali/gali.',modelname,sep=""))
            mouts_gali <- mouts
            rm(mouts)
           load(paste('results/mrsc/mrsc.',modelname,sep=""))
            mouts_mrsc <- mouts
            rm(mouts)
           load(paste('results/chpt/chpt.',modelname,sep=""))
            mouts_chpt <- mouts 
            rm(mouts)
           load(paste('results/dcsr/dcsr.',modelname,sep=""))
            mouts_dcsr <- mouts 
            rm(mouts) 
           load(paste('results/logue/results/Miller.',modelname,sep=""))
            mouts_miller <- mouts 
            rm(mouts)
          }
                       
          if(effect == "Int") {
           load(paste('results/gali/gali.',modelname,sep=""))
            mouts_gali <- mouts_keller
            rm(mouts_keller)  
           load(paste('results/mrsc/mrsc.',modelname,sep=""))
            mouts_mrsc <- mouts_keller
            rm(mouts_keller)
           load(paste('results/chpt/chpt.',modelname,sep=""))
            mouts_chpt <- mouts_keller 
            rm(mouts_keller)
           load(paste('results/dcsr/dcsr.',modelname,sep=""))
            mouts_dcsr <- mouts_keller 
            rm(mouts_keller)             
           load(paste('results/logue/results/Miller.',modelname,sep=""))
            mouts_miller <- mouts_keller 
            rm(mouts_keller)
           }
             
             


#Find the most significant term

#Regular models are called mouts
#Interaction keller models are called mouts_keller

#interaction should be the last (max) row 

 
           #Get the best Z score overall
           zscoremat <- rep(NA,101)
           
           for (i in c(1:101)){
            #The last term should be the interaction
            #PTSD will be second ot last in non interaction models.
            galidim <- dim(mouts_gali[[101]])[1]
            mrscdim <- dim(mouts_mrsc[[101]])[1]
            dcsrdim <- dim(mouts_dcsr[[101]])[1]
            chptdim <- dim(mouts_chpt[[101]])[1]
            millerdim <- dim(mouts_miller[[101]])[1]  

            galiweight <-sqrt (300)
            mrscweight <- sqrt(2333)
            millerweight <- sqrt(387)
            chptweight <-sqrt(96) # consider redoing, different studies??
            dcsrweight <- sqrt(49)
            
            weightsz <- c(galiweight,mrscweight,millerweight,dcsrweight,chptweight)
            models <- rbind(mouts_gali[[i]][galidim,] ,  mouts_mrsc[[i]][mrscdim,],  mouts_miller[[i]][millerdim,],mouts_dcsr[[i]][dcsrdim,],mouts_chpt[[i]][chptdim,]) #   

            weightdenom <- sqrt(sum(weightsz^2))
            weightnumer <- sum(models[,3]*weightsz)
            zscoremat[i] <- weightnumer/weightdenom
           }
           
           write.csv(zscoremat,file=paste('results_meta/meta.',modelname,'.csv',sep=''))
       }
      }
     } 
    }
   }
  }
 }
 
 
 #rma(yi=models[,1],sei=models[,2],method="FE")
 
ptsdrow <- max(
 grep("ptsd",row.names(mouts_keller[[1]])) 
  )

statmatrix <- as.data.frame(matrix(nrow=101,ncol=3))
names(statmatrix) <- c("Beta","SE","z")

#add in the row number because it corresopnds directly to prs 
for (i in 1:101)
{
 statmatrix[i,] <- mouts_keller[[i]][ptsdrow,]
 ptsdrowname <- row.names(mouts_keller[[i]])[ptsdrow]
  
 row.names(statmatrix)[i] <- ptsdrowname
}

statmatrix$P <- pchisq(statmatrix$z^2,1,lower.tail=F)

head(statmatrix[order(statmatrix$P),] )

#Extract the PTSD row from each meta analysis, make a matrix of beta se and p??


#The next step would be to meta-analyze these and pick the best out of the meta-analysis. 


 library(metafor)
 
 rma(yi=c(326,238),sei=c(183,238),method="FE")
 
 #Do I need ot select same pT for both?? seems like it for the thing to be remotely interpretable...@
 #for non dx variabes, i also need measure scale the same
 