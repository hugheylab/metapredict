#' @import ggplot2
#' @import methods
#' @importFrom foreach foreach
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom magrittr "%>%"
#' @importMethodsFrom AnnotationDbi mappedkeys
#' @importMethodsFrom Biobase experimentData
#' @importMethodsFrom Biobase exprs
#' @importMethodsFrom Biobase exprs<-
#' @importMethodsFrom Biobase featureData
#' @importMethodsFrom Biobase pData
#' @importMethodsFrom Biobase phenoData
#' @importMethodsFrom BiocGenerics as.list
NULL

globalVariables(c('studyName', 'coefSparse', 'validationStudyName', 'classLevel'))


fixCustomCdfGeneIds = function(geneIds) {
	return(sub('_at', '', geneIds))}


fixGeoSampleNames = function(sampleNames) {
	sampleNames = paste0(toupper(sampleNames), '_')
	regexResult = regexpr('^GSM[0-9]+[^0-9]', sampleNames)
	sampleNamesNew = mapply(function(sampleName, matchLength) substr(sampleName, 1, matchLength-1), sampleNames,
									attr(regexResult, 'match.length'))
	return(unname(sampleNamesNew))}


fixCelSampleNames = function(sampleNames) {
	sampleNamesNew = gsub('\\.cel(\\.gz)*$', '', sampleNames, ignore.case=TRUE)
	return(sampleNamesNew)}


getGeneProbeMappingAffy = function(mappingFilePath) {
	mapping = utils::read.table(mappingFilePath, sep='\t', header=TRUE, stringsAsFactors=FALSE)
	mappingUnique = unique(mapping[,c('Probe.Set.Name', 'Affy.Probe.Set.Name')])
	mappingUnique = mappingUnique[apply(mappingUnique, MARGIN=1, function(r) !any(is.na(r))),]
	rownames(mappingUnique) = NULL
	colnames(mappingUnique) = c('geneId', 'probeSet')
	mappingUnique[,'probeSet'] = as.character(mappingUnique[,'probeSet'])
	return(mappingUnique)}


getGeneProbeMappingDirect = function(featureDf, geneColname, probeColname='ID') {
	mapping = featureDf[,c(probeColname, geneColname)]
	mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
	mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
	colnames(mapping) = c('probeSet', 'geneId')
	return(mapping)}


getGeneProbeMappingAnno = function(featureDf, dbName, interName) {
	mappingProbeIntermediate = featureDf[!is.na(featureDf[,interName]) & featureDf[,interName]!='', c('ID', interName)]
	colnames(mappingProbeIntermediate) = c('probeSet', 'geneInter')

	mapTmp1 = eval(parse(text=sprintf('%s.db::%s', substr(dbName, 1, 9), dbName)))
	mapTmp2 = mappedkeys(mapTmp1)
	mapTmp3 = as.list(mapTmp1[mapTmp2])
	geneId = do.call(c, mapTmp3)

	geneInter = do.call(c, mapply(function(inter, len) rep_len(inter, len), names(mapTmp3), sapply(mapTmp3, length),
											SIMPLIFY=FALSE))
	if (dbName=='org.Hs.egUNIGENE2EG') {
		geneInter = sub('Hs.', '', geneInter, fixed=TRUE)}
	mappingIdInter = data.frame(geneId, geneInter, stringsAsFactors=FALSE)
	return(merge(mappingIdInter, mappingProbeIntermediate, by='geneInter', sort=FALSE))}


calcExprsByGene = function(eset, mapping) {
	geneIds = unique(mapping[,'geneId'])
	exprsByGene = matrix(nrow=length(geneIds), ncol=ncol(eset), dimnames=list(geneIds, Biobase::sampleNames(eset)))
	for (geneId in geneIds) {
		exprsTmp = exprs(eset)[mapping[mapping[,'geneId']==geneId, 'probeSet'],, drop=FALSE]
		if (nrow(exprsTmp)==1) {
			exprsByGene[geneId,] = exprsTmp
		} else {
			exprsByGene[geneId,] = matrixStats::rowMedians(t(exprsTmp), na.rm=TRUE)}}
	return(exprsByGene)}


#' Get the GPLs for microarray platforms that are currently supported.
#'
#' @return A character array of currently supported GPLs.
#'
#' @export
getSupportedPlatforms = function() {
	return(c('GPL180', 'GPL341', 'GPL571', 'GPL885', 'GPL887', 'GPL890', 'GPL962', 'GPL1053', 'GPL1073',
				'GPL1291', 'GPL1293', 'GPL1390', 'GPL1708', 'GPL3921', 'GPL5645', 'GPL6254', 'GPL6333',
				'GPL6480', 'GPL6865', 'GPL6880', 'GPL6884', 'GPL6885', 'GPL6887', 'GPL6947', 'GPL7015',
				'GPL7202', 'GPL8177', 'GPL10332', 'GPL10379', 'GPL10687', 'GPL13607', 'GPL15331', 'GPL15450',
				'GPL18721', 'GPL20769'))}


#' Get the GPLs for unsupported microarray platforms.
#'
#'\code{getUnsupportedPlatforms} checks if any studies in studyMetadata
#' whose \code{studyDataType=='series_matrix'} have \code{platformInfo}
#' that is not supported.
#'
#' @param studyMetadata data.frame of study metadata, with columns for \code{studyDataType} and \code{platformInfo}
#'
#' @export
getUnsupportedPlatforms = function(studyMetadata) {
	platformInfos = studyMetadata[studyMetadata[,'studyDataType']=='series_matrix', 'platformInfo']
	unsupportedPlatforms = sort(unique(platformInfos[!(platformInfos %in% getSupportedPlatforms())]))
	if (length(unsupportedPlatforms)==0) {
		cat("Whew, all microarray platforms for studies whose studyDataType=='series_matrix' are supported.\n")
	} else {
		warningStr = "Uh oh, not all microarray platforms for studies whose studyDataType=='series_matrix' are supported.
						  You will need to either drop those studies or extend the metapredict package.
						  See the vignette for details."
		warning(gsub('\\s+', ' ', warningStr), call.=FALSE)}
	return(unsupportedPlatforms)}


#' Get the gene expression data for one study.
#'
#' \code{getStudyData} loads the gene expression data for one study,
#' 	maps probes to Entrez Gene IDs, and normalizes and transforms
#' 	the expression values.
#'
#' @param parentFolderPath Path to the folder that contains the data.
#' @param studyName Name of dataset to load.
#' @param studyDataType Type of data.
#' @param platformInfo Microarray platform.
#'
#' @return An ExpressionSet, unless the platform or data type is not supported, then NA.
#'
#' @export
getStudyData = function(parentFolderPath, studyName, studyDataType, platformInfo) {
	cat(sprintf('Loading study %s...\n', studyName))
	if (studyDataType %in% c('affy_geo', 'affy_custom')) {
		require(platformInfo, character.only=TRUE)

		folderPath = file.path(parentFolderPath, studyName)
		if (dir.exists(folderPath)) {
			cwd = setwd(folderPath)
			eset = affy::justRMA(cdfname=platformInfo)
			setwd(cwd)
		} else {
			stop(sprintf('Folder %s does not exist.', folderPath))}

		Biobase::featureNames(eset) = fixCustomCdfGeneIds(Biobase::featureNames(eset))
		if (studyDataType=='affy_geo') {
			Biobase::sampleNames(eset) = fixGeoSampleNames(Biobase::sampleNames(eset))
		} else {
			Biobase::sampleNames(eset) = fixCelSampleNames(Biobase::sampleNames(eset))}

	} else if (studyDataType=='affy_series_matrix') {
		mapping = getGeneProbeMappingAffy(file.path(parentFolderPath, paste0(platformInfo, '_mapping.txt')))
		esetOrig = GEOquery::getGEO(filename=file.path(parentFolderPath, paste0(studyName, '_series_matrix.txt')))
		exprs(esetOrig)[exprs(esetOrig)<=0] = min(exprs(esetOrig)[exprs(esetOrig)>0])
		exprs(esetOrig) = log2(exprs(esetOrig))
		exprsByGene = calcExprsByGene(esetOrig, mapping)
		rownames(exprsByGene) = fixCustomCdfGeneIds(rownames(exprsByGene))
		colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
		eset = Biobase::ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig),
												experimentData=experimentData(esetOrig))

	} else if (studyDataType=='series_matrix') {
		supportedPlatforms = getSupportedPlatforms()

		if (!(platformInfo %in% supportedPlatforms)) {
			warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
			return(NA)}

		esetOrig = GEOquery::getGEO(filename=file.path(parentFolderPath, paste0(studyName, '_series_matrix.txt')))
		if (is.list(esetOrig) && length(esetOrig)==1) {
			esetOrig = esetOrig[[1]]}

		featureDf = pData(featureData(esetOrig))
		idx = sapply(featureDf, is.factor)
		featureDf[idx] = lapply(featureDf[idx], as.character)
		if (platformInfo=='GPL180') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='GENE_SYM')
		} else if (platformInfo=='GPL341') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='ENTREZ_GENE_ID')
		} else if (platformInfo=='GPL571') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='ENTREZ_GENE_ID')
		} else if (platformInfo=='GPL885') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL887') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL890') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL962') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egUNIGENE2EG', interName='UNIGENE')
		} else if (platformInfo=='GPL1053') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='GENE')
		} else if (platformInfo=='GPL1073') {
			featureDf[,'GenBank'] = sapply(featureDf[,'GB_ACC'], function(x) strsplit(x, split='.', fixed=TRUE)[[1]][1])
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Mm.egACCNUM2EG', interName='GenBank')
		} else if (platformInfo=='GPL1291') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez Gene ID')
		} else if (platformInfo=='GPL1293') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez Gene ID')
		} else if (platformInfo=='GPL1390') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_ACC')
		} else if (platformInfo=='GPL1708') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL3921') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='ENTREZ_GENE_ID')
		} else if (platformInfo=='GPL5645') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='Gene Name')
		} else if (platformInfo=='GPL6254') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egENSEMBL2EG', interName='ENSEMBL_GENE_ID')
		} else if (platformInfo=='GPL6333') {
			featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'], function(x) strsplit(x, split='.', fixed=TRUE)[[1]][1])
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Mm.egREFSEQ2EG', interName='RefSeq')
		} else if (platformInfo=='GPL6480') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL6865') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='rep_name')
		} else if (platformInfo=='GPL6880') {
			featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'], function(x) strsplit(x, split='.', fixed=TRUE)[[1]][1])
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Mm.egREFSEQ2EG', interName='RefSeq')
		} else if (platformInfo=='GPL6884') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez_Gene_ID')
		} else if (platformInfo=='GPL6885') {
			featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'], function(x) strsplit(x, split='.', fixed=TRUE)[[1]][1])
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Mm.egREFSEQ2EG', interName='RefSeq')
		} else if (platformInfo=='GPL6887') {
			featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'], function(x) strsplit(x, split='.', fixed=TRUE)[[1]][1])
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Mm.egREFSEQ2EG', interName='RefSeq')
		} else if (platformInfo=='GPL6947') {
			featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'], function(x) strsplit(x, split='.', fixed=TRUE)[[1]][1])
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='RefSeq')
		} else if (platformInfo=='GPL7015') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_LIST')
		} else if (platformInfo=='GPL7202') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL8177') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_ACC')
		} else if (platformInfo=='GPL10332') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_ACC')
		} else if (platformInfo=='GPL10379') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='EntrezGeneID')
		} else if (platformInfo=='GPL10687') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='EntrezGeneID')
		} else if (platformInfo=='GPL13607') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_ACC')
		} else if (platformInfo=='GPL15331') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_ACC')
		} else if (platformInfo=='GPL15450') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Dr.egREFSEQ2EG', interName='GB_ACC')
		} else if (platformInfo=='GPL18721') {
			featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'], function(x) strsplit(x, split='.', fixed=TRUE)[[1]][1])
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='RefSeq')
		} else if (platformInfo=='GPL20769') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='ORF')
		} else {
			warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
			return(NA)}

		exprsByGene = calcExprsByGene(esetOrig, mapping)
		if (any(is.na(exprsByGene))) {
			warning(sprintf('Imputing missing expression values for study %s.', studyName))
			resultImputed = impute::impute.knn(exprsByGene)
			exprsByGene = resultImputed$data}
		colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
		eset = Biobase::ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig),
												experimentData=experimentData(esetOrig))

	} else if (studyDataType=='eset_rds') {
		esetOrig = readRDS(file.path(parentFolderPath, paste0(studyName, '.rds')))

		featureDf = pData(featureData(esetOrig))
		if (platformInfo=='ready') {
			return(esetOrig)
		} else if (platformInfo=='rosetta') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='EntrezGene.ID', probeColname='probe')
		} else {
			warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
			return(NA)}

		exprsByGene = calcExprsByGene(esetOrig, mapping)
		if (any(is.na(exprsByGene))) {
			warning(sprintf('Imputing missing expression values for study %s.', studyName))
			resultImputed = impute::impute.knn(exprsByGene)
			exprsByGene = resultImputed$data}
		eset = Biobase::ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig),
												experimentData=experimentData(esetOrig))

	} else {
		warning(sprintf('Study %s not loaded, because data type %s is not currently supported.', studyName, studyDataType))
		eset = NA}
	return(eset)}


#' Get the gene expression data for multiple studies.
#'
#' \code{getStudyDataList} runs \code{getStudyData} for multiple studies.
#'
#' @param parentFolderPath Path to the folder that contains the data.
#' @param studyMetadata A data.frame where each row corresponds to one study.
#' 	rownames should correspond to study names, with columns for
#' 	studyDataType and platformInfo.
#'
#' @return A named list of ExpressionSets.
#'
#' @export
getStudyDataList = function(parentFolderPath, studyMetadata) {
	esetList = foreach(studyName=rownames(studyMetadata)) %do% {
		if (any(is.na(studyMetadata[studyName,]))) {
			NA
		} else {
			getStudyData(parentFolderPath, studyName, studyMetadata[studyName, 'studyDataType'],
							 studyMetadata[studyName, 'platformInfo'])}}
	names(esetList) = rownames(studyMetadata)
	return(esetList[!is.na(esetList)])}


#' Extract the expression matrices containing the desired samples.
#'
#' @param esetList list of ExpressionSets.
#' @param sampleMetadata data.frame with columns for study and sample names.
#'
#' @return A named list of matrices.
#'
#' @export
extractExpressionData = function(esetList, sampleMetadata) {
	ematList = foreach(studyName=names(esetList)) %do% {
		keepIdx = colnames(esetList[[studyName]]) %in% sampleMetadata[sampleMetadata[,'study']==studyName, 'sample']
		exprs(esetList[[studyName]])[,keepIdx]}
	names(ematList) = names(esetList)
	return(ematList)}
