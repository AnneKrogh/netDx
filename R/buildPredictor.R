#' Run nested cross-validation on data 
#' @return symmetric matrix of size ncol(dat) (number of patients) containing
#' pairwise patient similarities
#' 
#' @details wrapper function to run netDx with nested cross-validation, 
#' with an inner loop of X-fold cross-validation and an outer loop of different
#' random splits of data into train and blind test. The user needs to supply
#' a custom function to create PSN, see createPSN_MultiData(). This wrapper
#' provides flexibility for designs with one or several heterogeneous data
#' types, and one or more ways of defining patient similarity. 
#' For example, designs it handles includes
#' 1) Single datatype, single similarity metric: Expression data -> pathways
#' 2) Single datatype, multiple metrics: Expression data -> pathways
#'	(Pearson corr) and single gene networks (normalized difference)
#' 3) Multiple datatypes, multiple metrics: Expression -> Pathways; 
#'	Clinical -> single or grouped nets
#' @param pheno (data.frame) sample metadata, must have ID and STATUS columns
#' @param dataList (list) keys are datatypes; values contain patient data
#' for the corresponding datatype. e.g. dataList[["rna"]] contains expression
#' matrix. Rows are units (e.g. genes, individual clinical variables) and 
#' columns are patients
#' @param groupList (list) keys are datatypes and values are lists indicating
#' how units for those datatypes are to be grouped. Keys must match those 
#' in dataList. Each entry of groupList[[k]] will generate a new PSN.
#'  e.g. groupList[["rna"]] could be a list of pathway definitions. 
#' So keys(groupList[["rna"]]) would have pathway names, generating one PSN
#' per pathways, and values(groupList[["rna"]]) would be genes that would be
#' grouped for the corresponding pathwayList.
#' @param makeNetFunc (function) user-defined function for creating the set
#' of input PSN provided to netDx. See createPSN_MultiData()::customFunc.
#' @param outDir (char) directory where results will be stored. If this 
#' directory exists, its contents will be overwritten
#' @param trainProp (numeric 0 to 1) Percent samples to use for training
#' @param useMonteCarlo (logical) if TRUE use Monte Carlo for resampling
#'      and if FALSE use cross-validation in inner loop
#' @param innerSplitsMC (integer) Number of splits in Monte Carlo resampling in inner loop
#' @param featScoreMax (integer) determine size of training ((featScoreMax-1)/featScoreMax)
#' and test (1/featScoreMax) splits in inner loop. If CV is used for resampling featureScoreMax
#' is the number of folds in inner loop
#' @param numSplits (integer) number of train/blind test splits 
#' (i.e. iterations of outer loop)
#' @param numCores (integer) number of CPU cores for parallel processing
#' @param JavaMemory (integer) memory in (Gb) used for each fold of CV
#' @param featSelCutoff (integer) cutoff for inner-fold CV to call 
#' feature-selected in a given split
#' @param keepAllData (logical) if TRUE keeps all intermediate files, even
#' those not needed for assessing the predictor. Use very cautiously as for
#' some designs, each split can result in using 1Gb of data.
#' @param startAt (integer) which of the splits to start at (e.g. if the
#' job aborted part-way through)
#' @param preFilter (logical) if TRUE uses lasso to prefilter dataList within 
#' cross-validation loop. Only variables that pass lasso get included. The
#' current option is not recommended for pathway-level features as most genes
#' will be eliminated by lasso. Future variations may allow other prefiltering
#' options that are more lenient.
#' @param preFilterGroups (char) vector with subset of names(dataList)
#' to which prefiltering needs to be limited. Allows users to indicate
#' which data layers should be prefiltered using regression and which
#' are to be omitted from this process. Prefiltering uses regression, which
#' omits records with missing values. Structured missingness can result in
#' empty dataframes if missing values are removed from these, which in turn
#' can crash the predictor. To impute missing data, see the 'impute' and 
#' 'imputeGroups' parameters. 
#' @param preFilterSeed (integer) seed for lasso based pre-selection
#' @param impute (logical) if TRUE applies imputation by median within CV
#' @param imputeGroups (char) If impute set to TRUE, indicate which groups you 
#' want imputed. 
#' @param logging (char) level of detail with which messages are printed. 
#' Options are: 1) none: turn off all messages; 2) all: greatest level of 
#' detail (recommended for advanced users, or for debugging); 3) default: 
#' print key details (useful setting for most users)
#' @import glmnet
#' @importFrom stats median na.omit coef
#' @importFrom utils read.delim write.table
#' @return No value. Side effect of writing log messages to <outDir>/log.txt
#' and generating predictor-related data in <outDir>.
#' @export
#' @examples
#' load(sprintf("%s/extdata/buildPred_input.rda",
#'              path.package("netDx.examples")))
#' 
#' # Custom function to tell netDx how to build features (PSN) from input data.
#' # Each datatype is in a different entry of the dataList object.
#' # The sets into which variables should be grouped (e.g. genes into pathways)i
#' # is indicated in groupList.
#' # names(dataList) should correspond to names(groupList)
#' KIRC_makeNets <- function(dataList, groupList, netDir,...) {
#'     netList <- c()
#'     # make RNA nets: group by pathway, use default similarity metric (pearson corr)
#'     if (!is.null(groupList[["rna"]])) {
#'     netList <- makePSN_NamedMatrix(dataList$rna,
#'                     rownames(dataList$rna),
#'                 groupList[["rna"]],netDir,verbose=FALSE,
#'                 writeProfiles=TRUE,...)  
#'     netList <- unlist(netList)
#'     }
#' 
#'     # make clinical nets, one net for **each variable**.
#'   # use custom similarity metric
#'     netList2 <- c()
#'     if (!is.null(groupList[["clinical"]])) {
#'     netList2 <- makePSN_NamedMatrix(dataList$clinical,
#'         rownames(dataList$clinical),
#'         groupList[["clinical"]],netDir,
#'         simMetric="custom",customFunc=netDx::normDiff, # custom function
#'         writeProfiles=FALSE,
#'         sparsify=TRUE,verbose=TRUE,append=TRUE,...)
#'     }
#'     netList2 <- unlist(netList2)
#' 
#'     netList <- c(netList,netList2)  # concatenate net names as output
#'     return(netList)
#' }
#' 
#' # build predictor
#' # uncomment to run
#' # buildPredictor(pheno,
#' #   dataList=dats,groupList=groupList,
#' #   makeNetFunc=KIRC_makeNets, ### custom function defined above
#' #   outDir=sprintf("%s/pred_output",getwd()), ## absolute path
#' #   numCores=1L,featScoreMax=2L, featSelCutoff=1L,numSplits=2L)
buildPredictor <- function(pheno,dataList,groupList,outDir,makeNetFunc,
	trainProp=0.8,numSplits=10L,
        useMonteCarlo=TRUE,innerSplitsMC=10L,featScoreMax=10L,
        numCores,JavaMemory=4L,featSelCutoff=9L,
	keepAllData=FALSE,startAt=1L, preFilter=FALSE,preFilterSeed=123,
	impute=FALSE,
	preFilterGroups=NULL, imputeGroups=NULL,
	logging="default") { 

verbose_default <- TRUE
verbose_runQuery <- FALSE	  # messages when running individual queries
verbose_compileNets <- FALSE  # message when compiling PSN into database
verbose_runFS <- TRUE		  # runFeatureSelection() 
verbose_predict <- TRUE
verbose_compileFS <- FALSE
verbose_makeFeatures <- FALSE

if (logging == "all") {
	verbose_runQuery <- TRUE
	verbose_compileNets <- TRUE 
	verbose_compileFS <- TRUE
	verbose_makeFeatures <- TRUE
} else if (logging=="none") {
	verbose_runFS<-FALSE
	verbose_default <- FALSE
	verbose_predict <- FALSE
}

# Check input
if (missing(dataList)) stop("dataList must be supplied.\n")
if (missing(groupList)) stop("groupList must be supplied.\n")
tmp <- unlist(lapply(groupList,class))
not_list <- sum(tmp == "list")<length(tmp)
names_nomatch <- !identical(sort(names(groupList)),sort(names(dataList)))
if (!is(groupList,"list") || not_list || names_nomatch ) 
	stop("groupList must be a list of lists. Names must match those in dataList, and each entry should be a list of networks for this group.")
if (!is(dataList,"list")) 
	stop("dataList must be a list of data.frames or matrices, with names corresponding to groupList")
patnames <- unlist(lapply(dataList,colnames))
if (any(!(patnames %in% pheno$ID))) 
	stop("one or more patient IDs in dataList do not match those listed in pheno$ID. Check.")
if (trainProp <= 0 | trainProp >= 1) 
		stop("trainProp must be greater than 0 and less than 1")
if (startAt > numSplits) stop("startAt should be between 1 and numSplits")

megaDir <- outDir
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

# set aside for testing within each split
pheno_all <- pheno; 

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
cat("Predictor started at:\n")
print(Sys.time())
tryCatch({

# run featsel once per subtype
subtypes <- unique(pheno$STATUS)

if (verbose_default){
	cat(sprintf("-------------------------------\n"))
	cat(sprintf("# patients = %i\n", nrow(pheno)))
	cat(sprintf("# classes = %i { %s }\n", length(subtypes),
		paste(subtypes,collapse=",")))
	cat("Sample breakdown by class\n")
	print(table(pheno$STATUS))
	cat(sprintf("%i train/test splits",numSplits))
	cat(sprintf("Feature selection cutoff = %i of %i\n",
		featSelCutoff,featScoreMax))
	cat(sprintf("Datapoints:\n"))
	for (nm in names(dataList)) {
		cat(sprintf("\t%s: %i units\n", nm, nrow(dataList[[nm]])))
	}
}

# create master list of possible networks
netFile <- sprintf("%s/inputNets.txt", megaDir)

cat("NetType\tNetName\n",file=netFile)
for (nm in names(groupList)) {
	curNames <- names(groupList[[nm]])
	for (nm2 in curNames) {
		cat(sprintf("%s\t%s\n",nm,nm2),file=netFile,append=TRUE)
	}
}


if (verbose_default) {
	cat("\n\nCustom function to generate input nets:\n")
	print(makeNetFunc)
	cat(sprintf("-------------------------------\n\n"))
}

for (rngNum in startAt:numSplits) {
	if (verbose_default) {
	cat(sprintf("-------------------------------\n"))
	cat(sprintf("RNG seed = %i\n", rngNum))
	cat(sprintf("-------------------------------\n"))
	}
	outDir <- sprintf("%s/rng%i",megaDir,rngNum)
	dir.create(outDir)

	pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
								setSeed=rngNum*5,verbose=verbose_default)
	pheno <- subset(pheno_all, TT_STATUS %in% "TRAIN")
	dats_train <- lapply(dataList,function(x) { 
						 x[,which(colnames(x) %in% pheno$ID)]})

	if (impute) {
	if (verbose_default) cat("**** IMPUTING ****\n")
	if (is.null(imputeGroups)) imputeGroups <- names(dats_train)
	if (!any(imputeGroups %in% names(dats_train))) stop("imputeGroups must match names in dataList")

	
	dats_train <- sapply(names(dats_train), function(nm) {
		x <- dats_train[[nm]]
		if (nm %in% imputeGroups) {
			missidx <- which(rowSums(is.na(x))>0) 
			for (i in missidx) {
				na_idx <- which(is.na(x[i,]))
				x[i,na_idx] <- median(x[i,],na.rm=TRUE) 
			}
		} 
		x
	})
	}

	# prefilter with lasso
	if (preFilter) {
	if (is.null(preFilterGroups)) preFilter <- names(dats_train)
	if (!any(preFilterGroups %in% names(dats_train))) stop("preFilterGroups must match names in dataList")
	
	set.seed(preFilterSeed)
	cat("Prefiltering enabled\n")
	for (nm in preFilterGroups) {
		cat(sprintf("%s: %i variables\n",nm,nrow(dats_train[[nm]])))
		if (nrow(dats_train[[nm]])<2)  # only has one var, take it.
			vars <- rownames(dats_train[[nm]])
		else { 
			newx <- na.omit(dats_train[[nm]])
			tmp <- pheno[which(pheno$ID %in% colnames(newx)),]
			tryCatch( {
			fit <- cv.glmnet(x=t(newx),
					y=factor(tmp$STATUS), family="binomial", 
					alpha=1,maxit=1000000) # lasso
			}, error=function(ex) {
				print(ex)
				cat("*** You may need to set impute=TRUE for prefiltering ***\n")
			},finally={
			})
			wt <- abs(coef(fit,s="lambda.min")[,1])
			vars <- setdiff(names(wt)[which(wt>.Machine$double.eps)],
				"(Intercept)")
			}
		if (length(vars)>0) {
			tmp <- dats_train[[nm]]
			tmp <- tmp[which(rownames(tmp) %in% vars),,drop=FALSE]
			dats_train[[nm]] <- tmp
		} else {
			# leave dats_train as is, make a single net
		} 
		cat(sprintf("rngNum %i: %s: %s pruned\n",rngNum,nm,length(vars)))
		}
	}
	
	if (verbose_default) {
	cat("# values per feature (training)\n")
	for (nm in names(dats_train)) {
		cat(sprintf("\tGroup %s: %i values\n", 
			nm,nrow(dats_train[[nm]])))
	}
	}

	netDir <- sprintf("%s/networks",outDir)
	if (verbose_default) cat("** Creating features\n")
	createPSN_MultiData(dataList=dats_train,groupList=groupList,
			netDir=netDir,customFunc=makeNetFunc,numCores=numCores,
			verbose=verbose_makeFeatures)
	if (verbose_default) cat("** Compiling features\n")
	dbDir	<- compileFeatures(netDir, pheno$ID, outDir, numCores=numCores, 
				verbose=verbose_compileNets)
	

	if (verbose_default) cat("\n** Running feature selection\n")
	for (g in subtypes) {
	    pDir <- sprintf("%s/%s",outDir,g)
	    if (file.exists(pDir)) unlink(pDir,recursive=TRUE);
			dir.create(pDir)
			if (verbose_default) cat(sprintf("\tClass: %s",g))
			pheno_subtype <- pheno
			pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)] <- "nonpred"
			trainPred <- pheno_subtype$ID[which(pheno_subtype$STATUS %in% g)]
			if (verbose_default) 
				print(table(pheno_subtype$STATUS,useNA="always"))
		
			# Cross validation
			resDir <- sprintf("%s/GM_results",pDir)
			cat(sprintf("\tScoring features\n"))
			runFeatureSelection(trainPred, 
				outDir=resDir, dbPath=dbDir$dbDir, 
				nrow(pheno_subtype),verbose=verbose_runFS,
				useMonteCarlo=useMonteCarlo,nrOfSplits=innerSplitsMC,
				numCores=numCores, verbose_runQuery=TRUE, # verbose_runQuery,
				featScoreMax=featScoreMax,JavaMemory=JavaMemory)
	
	  	# Compute network score
			nrank <- dir(path=resDir,pattern="NRANK$")
			if (verbose_default) cat("\tCompiling feature scores\n")
			pTally		<- compileFeatureScores(paste(resDir,nrank,sep="/"),
				verbose=verbose_compileFS)
			tallyFile	<- sprintf("%s/%s_pathway_CV_score.txt",resDir,g)
			write.table(pTally,file=tallyFile,sep="\t",col=TRUE,row=FALSE,
				quote=FALSE)
		if (verbose_default) cat("\n")
		if (verbose_default) cat("\n")
	}
	
	## Class prediction for this split
	if (verbose_default) cat("\n** Predicting labels for test\n")
	pheno <- pheno_all
	predRes <- list()
	for (g in subtypes) {
		if (verbose_default) cat(sprintf("%s\n",g))
		pDir <- sprintf("%s/%s",outDir,g)
		pTally <- read.delim(
			sprintf("%s/GM_results/%s_pathway_CV_score.txt",pDir,g),
			sep="\t",header=TRUE,as.is=TRUE)
		idx <- which(pTally[,2]>=featSelCutoff)

		pTally <- pTally[idx,1]
		pTally <- sub(".profile","",pTally)
		pTally <- sub("_cont","",pTally)
		if (verbose_default)
			cat(sprintf("\t%i feature(s) selected\n",length(pTally)))
		netDir <- sprintf("%s/networks",pDir)

		dats_tmp <- list()
		for (nm in names(dataList)) {
			passed <- rownames(dats_train[[nm]])
			tmp <- dataList[[nm]]
			# only variables passing prefiltering should be used to make PSN
			dats_tmp[[nm]] <- tmp[which(rownames(tmp) %in% passed),] 
		}		

		# ------
		# Impute test samples if flag set
		# impute
		if (impute) {
		train_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% "TRAIN")]
		test_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% "TEST")]
		dats_tmp <- sapply(names(dats_tmp), function(nm) {
			x <- dats_tmp[[nm]]
			if (nm %in% imputeGroups) {
				missidx <- which(rowSums(is.na(x))>0) 
				train_idx <- which(colnames(x) %in% train_samp)
				test_idx <- which(colnames(x) %in% test_samp)
				for (i in missidx) {
					# impute train and test separately
					na_idx <- intersect(which(is.na(x[i,])),train_idx)
					na_idx1 <- na_idx
					x[i,na_idx] <- median(x[i,train_idx],na.rm=TRUE) 
		
					na_idx <- intersect(which(is.na(x[i,])),test_idx)
					na_idx2 <- na_idx
					x[i,na_idx] <- median(x[i,test_idx],na.rm=TRUE) 
				}
			}
			x
		})
		#alldat_tmp <- do.call("rbind",dats_tmp)
		}

		if (verbose_default) cat(sprintf("\tCreate & compile features\n",g))
		if (length(pTally)>=1) {
		createPSN_MultiData(dataList=dats_tmp,groupList=groupList,
			netDir=sprintf("%s/networks",pDir),
			customFunc=makeNetFunc,numCores=numCores,
			filterSet=pTally,verbose=verbose_default)
		dbDir <- compileFeatures(netDir,pheno$ID,pDir,numCores=numCores,
			verbose=verbose_compileNets)

		# run query for this class
		qSamps <- pheno$ID[which(pheno$STATUS %in% g & pheno$TT_STATUS%in%"TRAIN")]
		qFile <- sprintf("%s/%s_query",pDir,g)
		writeQueryFile(qSamps,"all",nrow(pheno),qFile)
		if (verbose_default) cat(sprintf("\t** %s: Compute similarity\n",g))
		resFile <- runQuery(dbDir$dbDir,qFile,resDir=pDir,
			JavaMemory=JavaMemory, numCores=numCores,
			verbose=verbose_runQuery)
		predRes[[g]] <- getPatientRankings(sprintf("%s.PRANK",resFile),pheno,g)
		} else {
			predRes[[g]] <- NA
		}
	}
	if (verbose_default) cat("\n")
	
	if (sum(is.na(predRes))>0 & verbose_default) {
		cat(sprintf("RNG %i : One or more classes have no selected features. Not classifying\n", rngNum))
	} else {
		if (verbose_default) cat("** Predict labels\n")
		predClass <- predictPatientLabels(predRes,
			verbose=verbose_predict)
		out <- merge(x=pheno_all,y=predClass,by="ID")
		outFile <- sprintf("%s/predictionResults.txt",outDir)
		write.table(out,file=outFile,sep="\t",col=TRUE,row=FALSE,quote=FALSE)
		
		acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
		if (verbose_default)
			cat(sprintf("Split %i: ACCURACY (N=%i test) = %2.1f%%\n",
			rngNum, nrow(out), acc*100))
	}
        
	if (!keepAllData) {
    system(sprintf("rm -r %s/dataset %s/tmp %s/networks",                       
        outDir,outDir,outDir))                                                  
	for (g in subtypes) {
    system(sprintf("rm -r %s/%s/dataset %s/%s/networks",
        outDir,g,outDir,g))
	}
	}# endif !keepAllData
	if (verbose_default) {
		cat("\n----------------------------------------\n")
	}
	}
}, error=function(ex){
	print(ex)
}, finally={
	cat("Predictor completed at:\n")
	print(Sys.time())
	sink(NULL)
})

}

