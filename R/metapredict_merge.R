makeMatchSampleMapping = function(metadata, subStudyNames, matchColname) {
	metadataNow = metadata[metadata[,'study'] %in% subStudyNames, c('study', 'sample', matchColname)]
	metadataNow = metadataNow[order(metadataNow[,'study'], decreasing=is.unsorted(subStudyNames)), c('sample', matchColname)]
	headFunc = function(x) x[1]
	mappingDf = metadataNow %>%
		dplyr::group_by_(.dots=list(matchColname)) %>%
		dplyr::summarise_each(dplyr::funs(headFunc)) %>%
		data.frame(check.names=FALSE)
	mapping = mappingDf[,'sample']
	names(mapping) = mappingDf[,matchColname]
	return(mapping)}


#' Merge gene expression from different platforms that was measured on the same biological samples.
#'
#' @param ematAtomicList list of expression matrices.
#' @param studyMetadataAtomic data.frame for study metadata,
#' 	with columns for study and sample names.
#' @param sampleMetadataAtomic data.frame for sample metadata,
#' 	with rownames corresponding to sample names.
#' @param matchColname column in sampleMetadata used to match samples.
#' @param mergeFunc function to summarize multiple gene expression values.
#'
#' @return A named list.
#' \item{ematList}{Named list of expression matrices.}
#' \item{studyMetadata}{data.frame of study metadata.}
#' \item{sampleMetadata}{data.frame of sample metadata.}
#'
#' @export
mergeMatchStudyData = function(ematAtomicList, studyMetadataAtomic, sampleMetadataAtomic, matchColname,
										 mergeFunc=function(x) mean(x, na.rm=TRUE)) {
	ematList = list()
	sampleMetadataList = list()

	for (matchStudyName in unique(studyMetadataAtomic[,'matchStudy'])) {
		if (sum(studyMetadataAtomic[,'matchStudy']==matchStudyName)==1) {
			ematList[[matchStudyName]] = ematAtomicList[[matchStudyName]]
			sampleMetadataList[[matchStudyName]] = sampleMetadataAtomic[sampleMetadataAtomic[,'study']==matchStudyName,]

		} else if (sum(studyMetadataAtomic[,'matchStudy']==matchStudyName)>1) {
			atomicStudyNames = studyMetadataAtomic[studyMetadataAtomic[,'matchStudy']==matchStudyName, 'study']
			edfListNow = list()
			for (atomicStudyName in atomicStudyNames) {
				edf = data.frame(rownames(ematAtomicList[[atomicStudyName]]), ematAtomicList[[atomicStudyName]])
				rownames(edf) = NULL
				colnames(edf) = c('geneId', sampleMetadataAtomic[colnames(edf)[2:ncol(edf)], matchColname])
				edfListNow[[atomicStudyName]] = edf}

			edfMerged = suppressWarnings(dplyr::bind_rows(edfListNow)) %>%
				dplyr::group_by_(.dots=list('geneId')) %>%
				dplyr::summarise_each(dplyr::funs(mergeFunc)) %>%
				data.frame(check.names=FALSE)
			rownames(edfMerged) = edfMerged[,'geneId']
			edfMerged = edfMerged[,-1]

			mapping = makeMatchSampleMapping(sampleMetadataAtomic, atomicStudyNames, matchColname)
			colnames(edfMerged) = mapping[colnames(edfMerged)]
			ematList[[matchStudyName]] = as.matrix(edfMerged)

			idx = (sampleMetadataAtomic[,'study'] %in% atomicStudyNames) &
				(sampleMetadataAtomic[,'sample'] %in% colnames(edfMerged))
			sampleMetadataList[[matchStudyName]] = sampleMetadataAtomic[idx,]
			sampleMetadataList[[matchStudyName]][,'study'] = matchStudyName}}

	headFunc = function(x) x[1]
	studyMetadata = studyMetadataAtomic %>%
		dplyr::group_by_(.dots=list('matchStudy')) %>%
		dplyr::summarise_each(dplyr::funs(headFunc)) %>%
		data.frame(check.names=FALSE)
	studyMetadata = studyMetadata[,colnames(studyMetadata)!='study']
	colnames(studyMetadata)[colnames(studyMetadata)=='matchStudy'] = 'study'
	rownames(studyMetadata) = studyMetadata[,'study']
	studyMetadata[,'matchStudy'] = studyMetadata[,'study']

	sampleMetadata = suppressWarnings(dplyr::bind_rows(sampleMetadataList)) %>% data.frame(check.names=FALSE)
	rownames(sampleMetadata) = sampleMetadata[,'sample']

	result = list(ematList, studyMetadata, sampleMetadata)
	names(result) = c('ematList', 'studyMetadata', 'sampleMetadata')
	return(result)}


#' Merge gene expression data from multiple studies.
#'
#' \code{mergeStudyData} merges gene expression data from
#' 	multiple studies, keeping only genes measured in each
#' 	dataset, then optionally performs cross-study normalization
#' 	using ComBat.
#'
#' @param ematList list of expression matrices.
#' @param sampleMetadata data.frame for sample metadata, with rownames
#' 	corresponding to sample names.
#' @param batchColname column in sampleMetadata containing
#' 	batch information for ComBat.
#' @param covariateName column in sampleMetadata containing
#' 	additional covariates for ComBat besides batch.
#' @param batchCorrection TRUE indicates cross-study normalization
#' 	will be performed, FALSE indicates datasets will be merged
#' 	without performing normalization.
#' @param parPrior passed to ComBat's \code{par.prior} argument.
#'
#' @return A matrix of expression for genes by samples.
#'
#' @export
mergeStudyData = function(ematList, sampleMetadata, batchColname='study', covariateName=NA, batchCorrection=TRUE,
								  parPrior=TRUE) {
	sampleNames = do.call(c, lapply(ematList, colnames))
	if (!all(sampleNames %in% rownames(sampleMetadata))) {
		stop('sampleMetadata must have rownames corresponding to the colnames of each matrix in ematList.', call.=FALSE)}

	geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
	ematList2 = foreach(studyName=names(ematList)) %do% {ematNow = ematList[[studyName]][geneIds,]}
	if (batchCorrection) {
		# if both one-color and two-color data is present, ComBat can fail catastrophically, if data is not scaled beforehand
		ematListScaled = lapply(ematList2, function(emat) (emat - mean(emat)) / stats::sd(emat))
		ematMerged = do.call(cbind, ematListScaled)
		if (is.na(covariateName) || length(unique(sampleMetadata[colnames(ematMerged), covariateName]))<2) {
			covariateInfo = stats::model.matrix(~rep_len(1, ncol(ematMerged)))
		} else {
			covariateInfo = stats::model.matrix(~sampleMetadata[colnames(ematMerged), covariateName])}

		if (length(unique(sampleMetadata[colnames(ematMerged), batchColname]))>1) {
			ematMergedNorm = sva::ComBat(ematMerged, batch=as.character(sampleMetadata[colnames(ematMerged), batchColname]),
												  mod=covariateInfo, par.prior=parPrior)
		} else {
			ematMergedNorm = ematMerged}

		return(ematMergedNorm)
	} else {
		return(do.call(cbind, ematList2))}}
