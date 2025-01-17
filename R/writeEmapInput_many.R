#' Wrapper to generate multiple EnrichmentMaps (perhaps one per class)
#'
#' @param featScores (list) keys are classes, and values are data.frames of
#' network scores across cross-validation (output of getFeatScores()).
#' @param namedSets_valid (list) Grouped unit variables limited to the
#' units contained in the dataset. e.g. keys are pathways and values are
#' the genes measured in this dataset.
#' e.g.:
#' $`MISSPLICED_GSK3BETA_MUTANTS_STABILIZE_BETA-CATENIN`
#' [1] "PPP2R5E" "PPP2CB"  "APC"     "AXIN1"   "PPP2R1B" "PPP2R1A" "CSNK1A1"
#' [8] "PPP2R5D" "PPP2R5C" "PPP2R5B" "PPP2R5A" "PPP2CA"  "GSK3B"
#' @param netTypes (data.frame) "inputNets.txt" file
#' generated by NetDx. Dataframe has two columns, network type and
#' network  name. I.E:
#'  clinical                                          clinical
#'       rna GUANOSINE_NUCLEOTIDES__I_DE_NOVO__I__BIOSYNTHESIS
#'       rna                              RETINOL_BIOSYNTHESIS
#' @param outDir (char) path to output directory
#' @param ... parameters for writeEMap()
#' @return (list)
#' 1) <outPfx>.gmt file - for enrichment map
#' 2) <outPfx>_nodeAttr.txt (file) table with node properties net name, max
#' score and net type
#' @examples
#' data(featScores)
#' pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
#'           path.package("netDx.examples"))
#' pathwayList <- readPathways(pathFile)
#' pathwayList <- pathwayList[c(1:5)]
#' netInfoFile <- sprintf("%s/extdata/KIRC_output/inputNets.txt",
#'      path.package("netDx.examples"))
#' netTypes <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
#' outDir <- paste(tempdir(),"plots",sep="/")
#' if (!file.exists(outDir)) dir.create(outDir)
#' EMap_input <- writeEMapInput_many(featScores,pathwayList,
#'      netTypes,outDir=outDir)
#' @export
writeEMapInput_many <- function(featScores, namedSets_valid, netTypes,
	outDir,...){

gmt_attr_files <- list()
for (gp in names(featScores)) {
	cur_out_files <- writeEMapInput(featScores[[gp]],namedSets_valid,netTypes,
			outPfx=sprintf("%s/%s",outDir,gp),...)

  gmt_attr_files[[gp]] <- cur_out_files

}
return(gmt_attr_files)
}
