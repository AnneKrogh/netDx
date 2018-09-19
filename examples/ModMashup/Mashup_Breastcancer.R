#' @author Guodong Xu, Shraddha Pai
rm(list=ls())

# Change this to a local directory where you have write permission
outDir <- sprintf("%s/TCGA_BRCA_mashup",getwd())
cat(sprintf("All intermediate files are stored in:\n%s\n",outDir))

numCores 	<- 8L  	# num cores available for parallel processing
GMmemory 	<- 4L  	# java memory in Gb
cutoff		<- 9L  	# score cutoff for feature-selected networks
TRAIN_PROP <- 0.67 	# fraction of samples to use for training

if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir, recursive=TRUE)

logFile <- sprintf("%s/log.txt",outDir)
sink(logFile,split=TRUE)
cat(sprintf("Logs will be stored in:\n%s\n",logFile))

# import the required packages
require("netDxmashup")
require(netDx.examples)
data(TCGA_BRCA)

#Split the train and test
subtypes<- c("LumA")
pheno$STATUS[which(!pheno$STATUS %in% subtypes)] <- "other"
subtypes <- c(subtypes,"other") # add residual

pheno$TT_STATUS <- splitTestTrain(pheno,
                                  pctT = TRAIN_PROP,setSeed = 42)

# Create similairty network
pheno_FULL	<- pheno
xpr_FULL 	<- xpr
pheno		<- subset(pheno,TT_STATUS %in% "TRAIN")
xpr			<- xpr[,which(colnames(xpr)%in% pheno$ID)]


## Pathway
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
                    path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
head(pathwayList)

## Gene data networks
profDir <- sprintf("%s/profiles",outDir)
netDir <- sprintf("%s/networks",outDir)


netList <- makePSN_NamedMatrix(xpr, rownames(xpr), 
                               pathwayList, profDir, verbose=FALSE, 
                               numCores=numCores, writeProfiles=FALSE, 
                               sparsify = TRUE)

netList <- unlist(netList)
head(netList)

# Feature selection for each class
top_net_file <- list()
mashup_tally <- list()
for (g in subtypes) {
  pDir <- sprintf("%s/%s",outDir,g)
  if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
  
  dir.create(pDir)
  
  cat(sprintf("\n******\nSubtype %s\n",g))
  pheno_subtype <- pheno
  
  ## label patients not in the current class as a residual
  pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)] <- "nonpred"
  ## sanity check
  print(table(pheno_subtype$STATUS,useNA="always"))
  
  Mashup_resDir <- sprintf("%s/Mashup_results",pDir)
  ## query for feature selection comprises of training 
  ## samples from the class of interest
  trainPred <- pheno$ID[which(pheno$STATUS %in% g)]
  
  # Cross validation for mashup
  # remember to set keyword write_query = FALSE if you want to uncomment 
  # GeneMANIA algorithm,
  # which indicates mashup will use query file from genemania instead of 
  # generating query files by itself, so the query files are shared 
  # between genemania and 
  # mashup for further comparation.
  # mashup_res <- mashup_runCV_featureSet(profDir, Mashup_resDir, pheno_subtype, trainID_pred = trainPred,
  mashup_res <- mashup_runCV_featureSet(profDir, GM_resDir, 
                                        pheno_subtype, trainID_pred = trainPred,
                                        write_query = TRUE, smooth = TRUE, verbose=T, 
                                        numCores = numCores, cut_off = cutoff)
  
  # List of selected top networks name
  mashup_tally[[g]] <- mashup_res$tally
  # Selected top networks txt file name
  top_net_file[[g]] <- mashup_res$top_net
  cat(sprintf("Mashup-%s: %i networks\n",g,length(mashup_tally[[g]])))
}

# Rank test patients using trained model
pheno <- pheno_FULL
predRes_mashup <- list()
for (g in subtypes) {
  pDir <- sprintf("%s/%s",outDir,g)
  # get feature selected net names
  profDir_mashup <- sprintf("%s/profiles_mashup",pDir)
  
  # prepare nets for net mashup db
  tmp <- makePSN_NamedMatrix(xpr_FULL, rownames(xpr_FULL), 
                             pathwayList[which(names(pathwayList)%in% mashup_tally[[g]])],
                             profDir_mashup,verbose=FALSE,
                             numCores=numCores,writeProfiles=FALSE, sparsify = TRUE)
  
  # Delete existent result file in case of conflicts.
  redundant_result_file <- list.files(path = sprintf("%s", pDir),
                                      pattern = "query")
  unlink(paste0(pDir, "/",redundant_result_file))
  
  # query of all training samples for this class
  qSamps <- pheno$ID[which(pheno$STATUS %in% g & pheno$TT_STATUS%in%"TRAIN")]
  qFile <- sprintf("%s/%s_query",pDir,g)
  GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)
  
  ### TODO  - mention that you're activating label prop mode for runMashup()
  ### by setting ranking=TRUE. THis is important.
  
  # Running patient ranking for mashup
  mashup_resFile <- runMashup(profDir_mashup, qFile, pheno, 
                              top_net = top_net_file[[g]], ranking = TRUE, 
                              smooth = TRUE)
  # Save the reresult.
  predRes_mashup[[g]] <- GM_getQueryROC(mashup_resFile, pheno, g, 
                                        plotIt=TRUE)
}

# Stats for Mashup patients prediction result.
predClass_mashup <- GM_OneVAll_getClass(predRes_mashup)
cat("Start Print result of mashup..")
both <- merge(x=pheno,y=predClass_mashup,by="ID")
print(table(both[,c("STATUS","PRED_CLASS")]))
pos <- (both$STATUS %in% "LumA")
tp <- sum(both$PRED_CLASS[pos]=="LumA")
fp <- sum(both$PRED_CLASS[!pos]=="LumA")
tn <- sum(both$PRED_CLASS[!pos]=="other")
fn <- sum(both$PRED_CLASS[pos]=="other")
cat(sprintf("Accuracy = %i of %i (%i %%)\n",tp+tn,nrow(both),
            round(((tp+tn)/nrow(both))*100)))
cat(sprintf("PPV = %i %%\n", round((tp/(tp+fp))*100)))
cat(sprintf("Recall = %i %%\n", round((tp/(tp+fn))*100)))  
