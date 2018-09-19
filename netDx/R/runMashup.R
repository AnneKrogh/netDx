#' Run a Mashup feature selection or patients ranking query.
#'
#' @author Guodong Xu, Shraddha Pai
#' @param Mashup_db (char) path to directory with Mashup generic database
#' @param queries (char) for ranking, path to query file. for selection, path to query directory.
#' @param true_pheno (data.frame) sample information with ID and STATUS columns.
#' @param trainID_pred (char) vector of training sample ID,matching pheno$ID.
#' @param smooth (logical) perform smooth in the network or not.
#' @param verbose (logical) print messages
#' @param ranking (logical) rank patients or run feature selection.  *** 
#' changes the mode of this function. If ranking=TRUE, "queries" must be a 
#' single query file. If FALSE, "queries" must be a directory with queries 
#' for each cross-validation.
#' @param top_net (char) a file stores selected top networks for the 
#' interested type.
#' @param cut_off (integer) cutoff to eliminate redundant network through 
#' network tally, range from 0 to numCV.
#' @return if used for patients ranking, return path to Mashup PRANK file.
#' @examples
#' runMashup(Mashup_db, queries, true_pheno, trainID_pred = NULL,smooth = TRUE,
#' 	verbose = TRUE, ranking = TRUE, top_net = NULL, cut_off = 9) 
#' @export
runMashup <- function (Mashup_db, queries, true_pheno, trainID_pred = NULL, 
		smooth = TRUE, verbose = TRUE, 
    ranking = TRUE, top_net = NULL, cut_off = 9) 
{
  # write id and labels file.
  if (!is.null(trainID_pred) ){
    true_pheno$STATUS[which(true_pheno$ID %in% trainID_pred)] <- 1
    true_pheno$STATUS[which(!true_pheno$ID %in% trainID_pred)] <- -1
    labels_file <- sprintf("%s/labels.txt", queries)
    if (verbose)
      cat(labels_file)
    write.table(true_pheno[c("ID", "STATUS")], 
             file = labels_file, col.names = FALSE, row.names = FALSE, 
							quote = FALSE)
  }
  id <- sprintf("%s/ids.txt", dirname(queries))
  if (verbose)
    cat(id)
  write.table(true_pheno["ID"], 
              file = id, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Check if want to smooth the similarity network.
  smooth_str <- ifelse(smooth, "true", "false")
  
  # Default value for cmd
  cmd <- sprintf("julia")
  mashup_julia <- sprintf("%s/julia/mashup.jl", path.package("netDxmashup"))
  
# ---------------------------------------------
# Mode 1: Label propagation for patient ranking
# ---------------------------------------------
  # If runnning for pantients ranking.
  if (ranking){
    # If top_net is null, then we do patient ranking with all networks without 
		# selection.
    if (is.null(top_net)){
      stopifnot(!dir.exists(queries))
      cmd <- sprintf("julia %s ranking --net %s --id %s --CV_query %s --smooth %s --res_dir %s",
           mashup_julia, Mashup_db, id, queries, smooth_str, 
					 dirname(queries))
    }
    else{
      # If top_net is not null, then we do patient ranking with selected networks.
      stopifnot(!dir.exists(queries))
    cmd <- sprintf("julia %s ranking --top_net %s --net %s --id %s --CV_query %s --smooth %s --res_dir %s",
                   mashup_julia, top_net, Mashup_db, id, queries, smooth_str, 
									dirname(queries))
    }
  }
# ---------------------------------------------
# Mode 2: Perform cross-validation and feature selection
# ---------------------------------------------
  # If running for network selection.
  else{
    stopifnot(dir.exists(queries))
    cmd <- sprintf("julia %s selection --net %s --id %s --labels %s --CV_query %s --smooth %s --cut_off %d --res_dir %s",
                 mashup_julia, Mashup_db, id, labels_file, queries, smooth_str,
								 cut_off, queries)
  }
  
  print(cmd)
  cat("\nRunning command:\n")
  system.time(system(cmd, wait = TRUE, ignore.stdout = !verbose, 
			ignore.stderr = !verbose))
  #cat(sprintf("QueryRunner time taken: %1.1f s\n", Sys.time() - t0))
  if (!dir.exists(queries)){
    if (smooth){
      res_file <- sprintf("%s_smooth_mashup_PRANK.txt", queries)
      res_file
    }
    else{
      res_file <- sprintf("%s_no_smooth_mashup_PRANK.txt", queries)
      res_file
    }
  }
}
