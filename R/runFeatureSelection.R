#' Run GeneMANIA cross-validation with a provided subset of networks
#'
#' @details Creates query files, runs GM for 10-fold cross validation.
#' @param trainID_pred (char) vector with universe of predictor class
#' patients (ie all that can possibly be included in the query file
#' @param outDir (char) directory to store query file and GM results
#' @param dbPath (char) path to GeneMANIA generic database with
#'	training population
#' @param numTrainSamps (integer) number of training samples in total
#' leave blank to use 5 training samples in order to save memory
#' @param incNets (char) vector of networks to include in this analysis
#' (features/pathway names). Useful for subset-based feature selection
#' @param orgName (char) organism name for GeneMANIA generic database.
#' The default value will likely never need to be changed.
#' @param fileSfx (char) file suffix
#' @param verbose (logical) print messages
#' @param useMonteCarlo (logical) if TRUE use Monte Carlo for resampling
#'      and if FALSE use cross-validation
#' @param nrOfSplits (integer) Number of times to run query in Monte Carlo resampling
#' @param featScoreMax (integer) determine size of training ((featScoreMax-1)/featScoreMax)
#' and test (1/featScoreMax) splits. In cross-validation featScoreMax is also the number of
#' Number of times to run query, usually equal to the max score for features in
#' the design (e.g. if featScoreMax=10, then this value is 10).
#' @param numCores (logical) num parallel threads for cross-validation
#' @param JavaMemory (integer) memory for GeneMANIA run, in Gb.
#' @param seed_queryResample (integer) RNG seed for inner cross validation loop.
#' Makes deterministic samples held-out for each GeneMANIA query (see
#' makeQueries())
#' @param verbose_runQuery (logical) print messages for runQuery()
#' @param ... args for \code{makeQueries()}
#' @return No value. Side effect of generating feature scores.
#' @examples
#' data(MB.pheno)
#' dbPath <- sprintf("%s/extdata/dbPath",path.package("netDx"))
#' runFeatureSelection(MB.pheno$ID[which(MB.pheno$STATUS%in% "WNT")],
#'	"~/tmp",dbPath,103L)
#' @export
runFeatureSelection <- function(trainID_pred,outDir,dbPath,numTrainSamps = NULL,
	incNets="all",orgName="predictor",fileSfx="CV",verbose=FALSE,
	useMonteCarlo=TRUE,featScoreMax=10L,nrOfSplits=10L,
	numCores=2L,JavaMemory=6L,seed_queryResample=42L,
	verbose_runQuery=FALSE,...) {

	#TODO if results already exist, what do we do? Delete with a warning?
	if (!file.exists(outDir)) dir.create(outDir)

	# get query names
	if (verbose) cat("\tWriting queries:\n")
	qSamps <- makeQueries(trainID_pred,verbose=verbose,
		useMonteCarlo=useMonteCarlo,nrOfSplits=nrOfSplits,featScoreMax=featScoreMax,
		setSeed=seed_queryResample,...)

	# write query files
	for (m in 1:length(qSamps)) {
		qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx,m)

		if(is.null(numTrainSamps)){
			numTrainSamps = 5
			cat("Memory saver option: using 5 training samples for CV")
		}

		writeQueryFile(qSamps[[m]], incNets, numTrainSamps,
						  qFile,orgName)
	}
	qFiles <- list()
	for (m in 1:length(qSamps)) {
		qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx, m)
		qFiles <- append(qFiles, qFile)
	}

	runQuery(dbPath, qFiles, outDir, JavaMemory=JavaMemory, 
			verbose=verbose_runQuery,
	       numCores=numCores)

}
