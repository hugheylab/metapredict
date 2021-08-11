makeMatchSampleMapping = function(metadata, subStudyNames, matchSampleColname) {

  study = ..cols = ..matchSampleColname = NULL
  cols = c('study', 'sample')
  if (!(matchSampleColname %in% cols)) cols = c(cols, matchSampleColname)
  metadataNow = data.table(metadata)[study %in% subStudyNames, ..cols]

  order = if (is.unsorted(subStudyNames)) -1L else 1L
  setorderv(metadataNow, 'study', order = order)

  mappingDf = metadataNow[, .SD[1L], by = ..matchSampleColname]
  mapping = mappingDf$sample
  names(mapping) = mappingDf[[matchSampleColname]]
  return(mapping)}


#' Merge gene expression from different platforms that was measured on the same biological samples.
#'
#' @param ematAtomicList list of expression matrices.
#' @param studyMetadataAtomic data frame for study metadata.
#' @param matchStudyColname column in studyMetadataAtomic used to match studies.
#' @param sampleMetadataAtomic data frame for sample metadata.
#' @param matchSampleColname column in sampleMetadataAtomic used to match samples.
#' @param mergeFunc function to summarize multiple gene expression values.
#'
#' @return A named list.
#' \item{ematList}{Named list of expression matrices.}
#' \item{studyMetadata}{data.frame of study metadata.}
#' \item{sampleMetadata}{data.frame of sample metadata.}
#'
#' @export
mergeMatchStudyData = function(
  ematAtomicList, studyMetadataAtomic, matchStudyColname, sampleMetadataAtomic,
  matchSampleColname, mergeFunc = function(x) mean(x, na.rm = TRUE)) {

  study = geneId = ..matchStudyColname = ..colnamesKeep = NULL
  ematList = list()
  sampleMetadataList = list()
  sampleMetadataAtomicDT = data.table(sampleMetadataAtomic)

  for (matchStudyName in unique(studyMetadataAtomic[[matchStudyColname]])) {
    if (sum(studyMetadataAtomic[[matchStudyColname]] == matchStudyName) == 1) {
      ematList[[matchStudyName]] = ematAtomicList[[matchStudyName]]
      sampleMetadataList[[matchStudyName]] = sampleMetadataAtomicDT[study == matchStudyName]

    } else if (sum(studyMetadataAtomic[[matchStudyColname]] == matchStudyName) > 1) {
      atomicStudyNames = studyMetadataAtomic$study[studyMetadataAtomic[[matchStudyColname]] == matchStudyName]
      edfListNow = list()
      for (atomicStudyName in atomicStudyNames) {
        edf = data.frame(rownames(ematAtomicList[[atomicStudyName]]), ematAtomicList[[atomicStudyName]])
        rownames(edf) = NULL

        matchNames = mergeDataTable(colnames(edf)[2:ncol(edf)], sampleMetadataAtomicDT)[[matchSampleColname]]
        colnames(edf) = c('geneId', matchNames)
        edfListNow[[atomicStudyName]] = edf}

      edfMerged = data.frame(data.table(
        rbind(edfListNow)[[1L]])[, lapply(.SD, mergeFunc), by = geneId],
        check.names = FALSE)

      rownames(edfMerged) = edfMerged$geneId
      edfMerged = edfMerged[, -1L]

      mapping = makeMatchSampleMapping(sampleMetadataAtomic, atomicStudyNames, matchSampleColname)
      colnames(edfMerged) = mapping[colnames(edfMerged)]
      ematList[[matchStudyName]] = as.matrix(edfMerged)

      idx = (sampleMetadataAtomic$study %in% atomicStudyNames) &
        (sampleMetadataAtomic$sample %in% colnames(edfMerged))
      sampleMetadataList[[matchStudyName]] = sampleMetadataAtomic[idx, ]
      sampleMetadataList[[matchStudyName]]$study = matchStudyName}}

  colnamesKeep = setdiff(colnames(studyMetadataAtomic), matchStudyColname)

  studyMetadata = studyMetadataAtomic[
    , .SD[1L], by = ..matchStudyColname][
      , study := matchStudyColname][, ..colnamesKeep]

  sampleMetadata = suppressWarnings(rbindlist(sampleMetadataList))

  result = list(ematList, studyMetadata, sampleMetadata)
  names(result) = c('ematList', 'studyMetadata', 'sampleMetadata')
  return(result)}


#' Merge gene expression data from multiple studies.
#'
#' Merge gene expression data from multiple studies, keeping only genes measured
#' in each dataset, then optionally perform cross-study normalization using
#' [sva::ComBat()].
#'
#' @param ematList list of expression matrices.
#' @param sampleMetadata data.frame for sample metadata.
#' @param batchColname column in sampleMetadata containing
#' 	 batch information for ComBat.
#' @param covariateName column in sampleMetadata containing
#' 	 additional covariates for ComBat besides batch.
#' @param batchCorrection `TRUE` indicates cross-study normalization
#' 	 will be performed, `FALSE` indicates datasets will be merged
#' 	 without performing normalization.
#' @param parPrior passed to ComBat's `par.prior` argument.
#'
#' @return A matrix of expression for genes by samples.
#'
#' @export
mergeStudyData = function(
  ematList, sampleMetadata, batchColname = 'study', covariateName = NA,
  batchCorrection = TRUE, parPrior = TRUE) {

  studyName = NULL
  sampleNames = do.call(c, lapply(ematList, colnames))
  if (!all(sampleNames %in% sampleMetadata$sample)) {
    stop('sampleMetadata must have samples corresponding to the colnames of each matrix in ematList.',
         call. = FALSE)}

  ematList = ematList[sapply(ematList, ncol) > 0]
  geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
  ematList2 = foreach(studyName = names(ematList)) %do% {
    ematNow = ematList[[studyName]][geneIds, ]}

  if (isTRUE(batchCorrection)) {
    # if both one-color and two-color data is present and data is not scaled beforehand,
    # ComBat can fail catastrophically
    ematListScaled = lapply(ematList2, function(emat) (emat - mean(emat)) / stats::sd(emat))
    ematMerged = do.call(cbind, ematListScaled)
    sm = mergeDataTable(colnames(ematMerged), sampleMetadata)

    if (is.na(covariateName) || length(unique(sm[[covariateName]])) < 2) {
      covariateInfo = stats::model.matrix(~ rep_len(1, ncol(ematMerged)))
    } else {
      covariateInfo = stats::model.matrix(~ sm[[covariateName]])}

    if (length(unique(sm[[batchColname]])) > 1) {
      ematMergedNorm = sva::ComBat(
        ematMerged, batch = as.character(sm[[batchColname]]),
        mod = covariateInfo, par.prior = parPrior)
    } else {
      ematMergedNorm = ematMerged}

    return(ematMergedNorm)
  } else {
    return(do.call(cbind, ematList2))}}
