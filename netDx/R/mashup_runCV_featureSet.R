#' Run Mashup cross-validation with a provided networks.
#'
#' @details Creates query files if choosen, runs Mashup for 10-fold cross
#' validation. 
#' @param Mashup_db (char) path to directory with network .profile files
#' @param queries_dir (char) directory where a list of query file stored 
#' @param trainID_pred (list) vector of training sample ID,matching pheno$ID.
#' STATUS columns
#' @param true_pheno (data.frame) phenotype table with ID and STATUS columns.
#' @param incNets (char) vector of networks to include in this analysis.
#' (features/pathway names). Useful for subset-based feature selection.
#' @param smooth (logical) perform smooth in the network or not.
#' value?
#' @param cut_off (integer) cutoff to eliminate redundant network through, 
#' network tally, default and recommendation is 9/10 * num_cv that is 9.
#' @param orgName (char) organism name for Mashup generic database.
#' The default value will likely never need to be changed.
#' @param write_query (logical) set to TRUE to have mashup generate its own 
#' query files. If FALSE looks for .query files in queries_dir
#' @param fileSfx (char) file suffix.
#' @param verbose (logical) print messages.
#' @param numCores (logical) num parallel threads for cross-validation.
#' @param seed_CVqueries (integer) RNG seed for inner cross validation loop.
#' Makes deterministic samples held-out for each mashup query (see
#' makeCVqueries())
#' @param ... args for \code{makeCVqueries()}
#' @examples
#' mashup_runCV_featureSet(Mashup_db, queries_dir, trainID_pred, true_pheno, 
#' 	  incNets = "all", smooth = TRUE, cut_off = 9, orgName = "predictor", 
#'		write_query = TRUE, fileSfx = "CV", verbose = FALSE, numCores = 2L, 
#			seed_CVqueries = 42L, ...) 
#' @export
mashup_runCV_featureSet <- function (Mashup_db, queries_dir, trainID_pred, 
		true_pheno, incNets = "all", smooth = TRUE, cut_off = 9, 
		orgName = "predictor", write_query = TRUE, fileSfx = "CV", 
		verbose = FALSE, numCores = 2L, seed_CVqueries = 42L, ...) 
{
  num_train_samps <- length(true_pheno)
  if (!file.exists(queries_dir)) 
    dir.create(queries_dir)
  if (write_query){
    if (verbose) cat("\tWriting GM queries: ")

    qSamps <- makeCVqueries(trainID_pred, verbose = verbose, 
                            setSeed = seed_CVqueries, ...)
    for (m in 1:length(qSamps)) {
      if (verbose) cat(sprintf("%i ", m))
      qFile <- sprintf("%s/%s_%i.query", queries_dir, fileSfx, m)
      GM_writeQueryFile(qSamps[[m]], incNets, num_train_samps, 
      	qFile, orgName)
    }
  }

### when ranking=FALSE switches to the 'feature selection' mode.
  runMashup(Mashup_db, queries_dir, true_pheno, trainID_pred = trainID_pred,
		smooth = smooth, ranking = FALSE, cut_off = cut_off, verbose = verbose)
  top_net = sprintf("top_networks")
  if(smooth)
    top_net <- sprintf("%s/smooth_result/top_networks.txt", queries_dir)
  if(!smooth)
    top_net <- sprintf("%s/no_smooth_result/top_networks.txt", queries_dir)
  mashupTally <- read.delim(top_net, header = FALSE)
  mashupTally <- sub("_cont", "", mashupTally[[1]])
  return(list(top_net = top_net, tally = mashupTally))
}

