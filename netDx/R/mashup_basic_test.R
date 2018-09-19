
rm(list=ls())


# Change this to a local directory where you have write permission
outDir <- sprintf("%s/TCGA_BRCA",getwd())
cat(sprintf("All intermediate files are stored in:\n%s\n",outDir))

numCores 	<- 2L  	# num cores available for parallel processing
GMmemory 	<- 4L  	# java memory in Gb
cutoff		<- 9L  	# score cutoff for feature-selected networks
TRAIN_PROP <- 0.67 	# fraction of samples to use for training

if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

# import the required packages
# require("netDxmashup")
require(netDx)
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
                               pathwayList,profDir,verbose=TRUE,
                               numCores=numCores,writeProfiles=TRUE) # use Genemania to generate interaction networks

netList <- unlist(netList)
head(netList)
