#' @import data.table
#' @import ggplot2
#' @import methods
#' @importFrom foreach foreach %do% %dopar%
#' @importMethodsFrom AnnotationDbi mappedkeys
#' @importMethodsFrom Biobase experimentData
#' @importMethodsFrom Biobase exprs
#' @importMethodsFrom Biobase exprs<-
#' @importMethodsFrom Biobase featureData
#' @importMethodsFrom Biobase pData
#' @importMethodsFrom Biobase phenoData
#' @importMethodsFrom BiocGenerics as.list
NULL

globalVariables(c(
  'studyName', 'coefSparse', 'validationStudyName', 'classLevel', '.', 'study',
  'ii', 'studyDataType', 'platformInfo', 'geneId', 'geneIdTmp', 'coefficient',
  'mapping', '..className', '..colnamesKeep', '..cols', '..geneColname',
  '..old', '..probeColname', 'idx', 'prob', 'probClass', 'probeSet', 'studyRow',
  'trueClass', 'trueClassProb', '..interName', '..matchSampleColname',
  '..matchStudyColname'))


fixCustomCdfGeneIds = function(geneIds) {
  return(sub('_at', '', geneIds))}


fixGeoSampleNames = function(sampleNames) {
  sampleNames = paste0(toupper(sampleNames), '_')
  regexResult = regexpr('^GSM[0-9]+[^0-9]', sampleNames)
  sampleNamesNew = mapply(
    function(sampleName, matchLength) substr(sampleName, 1, matchLength - 1),
    sampleNames, attr(regexResult, 'match.length'))
  return(unname(sampleNamesNew))}


fixCelSampleNames = function(sampleNames) {
  sampleNamesNew = gsub('\\.cel(\\.gz)*$', '', sampleNames, ignore.case = TRUE)
  return(sampleNamesNew)}


getGeneProbeMappingAffy = function(mappingFilePath) {
  mapping = fread(mappingFilePath)
  old = c('Probe.Set.Name', 'Affy.Probe.Set.Name')
  mappingUnique = unique(mapping[, ..old])
  mappingUnique = mappingUnique[apply(mappingUnique, MARGIN = 1, function(r) !any(is.na(r))), ]
  setnames(mappingUnique, old, c('geneId', 'probeSet'))
  mappingUnique[, probeSet := as.character(probeSet)]
  return(mappingUnique)}


getGeneProbeMappingDirect = function(featureDt, geneColname, probeColname = 'ID') {
  mapping = featureDt[, c(..probeColname, ..geneColname)]
  mapping = mapping[apply(mapping, MARGIN = 1, function(x) all(!is.na(x) & x!='')), ]
  setnames(mapping, c(probeColname, geneColname), c('probeSet', 'geneId'))
  mapping[, probeSet := as.character(probeSet)]
  mapping[, geneId := as.character(geneId)]
  return(mapping)}


getGeneProbeMappingAnno = function(featureDt, dbName, interName) {
  mappingProbeIntermediate = featureDt[
    !is.na(featureDt[[interName]]) & featureDt[[interName]]!='',
    c('ID', ..interName)]
  setnames(mappingProbeIntermediate, c('ID', interName), c('probeSet', 'geneInter'))

  mapTmp1 = eval(parse(text = sprintf('%s.db::%s', substr(dbName, 1, 9), dbName)))
  mapTmp2 = mappedkeys(mapTmp1)
  mapTmp3 = as.list(mapTmp1[mapTmp2])
  geneId = do.call(c, mapTmp3)

  geneInter = do.call(c, mapply(function(inter, len) rep_len(inter, len), names(mapTmp3),
                                sapply(mapTmp3, length), SIMPLIFY = FALSE))
  if (dbName == 'org.Hs.egUNIGENE2EG') {
    geneInter = sub('Hs.', '', geneInter, fixed = TRUE)}
  mappingIdInter = data.table(geneId, geneInter)
  mapping = merge(mappingIdInter, mappingProbeIntermediate, by = 'geneInter', sort = FALSE)
  mapping[, probeSet := as.character(probeSet)]
  return(mapping)}


calcExprsByGene = function(eset, mapping) {
  geneIds = unique(mapping$geneId)
  exprsByGene = matrix(nrow = length(geneIds), ncol = ncol(eset),
                       dimnames = list(geneIds, Biobase::sampleNames(eset)))
  for (geneIdNow in geneIds) {
    exprsTmp = exprs(eset)[mapping[geneId == geneIdNow]$probeSet, , drop = FALSE]
    exprsByGene[geneIdNow, ] = matrixStats::colMedians(exprsTmp, na.rm = TRUE)}
  return(exprsByGene)}


#' Get the GPLs for microarray platforms that are currently supported.
#'
#' @return A data.table of currently supported GPLs and relevant information.
#'
#' @export
getSupportedPlatforms = function() {

  supportedPlatformsPath = system.file('extdata', 'supported_platforms.csv', package = 'metapredict')
  supportedPlatformsDt = fread(supportedPlatformsPath)

  return(supportedPlatformsDt)}


#' Get the GPLs for unsupported microarray platforms.
#'
#' @param studyMetadata `data.frame` of study metadata, with columns for
#'   `studyDataType` and `platformInfo`.
#'
#' @return A vector of platforms for which `studyDataType` is 'series_matrix'
#'   and `platformInfo` is not supported.
#'
#' @export
getUnsupportedPlatforms = function(studyMetadata) {
  # unsupportedPlatforms = studyMetadata %>%
  #   dplyr::filter(studyDataType == 'series_matrix',
  #                 !(platformInfo %in% getSupportedPlatforms())) %>%
  #   .$platformInfo

  unsupportedPlatforms = data.table(studyMetadata)[
    studyDataType == 'series_matrix' & !(platformInfo %in% getSupportedPlatforms()$platform)]$platformInfo

  if (length(unsupportedPlatforms) == 0) {
    cat("Whew, all microarray platforms for studies whose studyDataType == 'series_matrix' are supported.\n")
  } else {
    warningStr = "Uh oh, not all microarray platforms for studies whose studyDataType == 'series_matrix' are supported.
    You will need to either drop those studies or extend the metapredict package.
    See the vignette for details."
    warning(gsub('\\s+', ' ', warningStr), call. = FALSE)}
  return(unsupportedPlatforms)}


#' Get the gene expression data for one study.
#'
#' Load the gene expression data for one study, map probes to Entrez Gene IDs,
#' then normalize and transform the expression values.
#'
#' @param parentFolderPath Path to the folder that contains the data.
#' @param studyName Name of dataset to load.
#' @param studyDataType Type of data.
#' @param platformInfo Microarray platform.
#'
#' @return An `ExpressionSet`, unless the platform or data type is not
#'   supported, then NA.
#'
#' @export
getStudyData = function(parentFolderPath, studyName, studyDataType, platformInfo) {
  cat(sprintf('Loading study %s...\n', studyName))
  if (studyDataType %in% c('affy_geo', 'affy_custom')) {
    require(platformInfo, character.only = TRUE)

    folderPath = file.path(parentFolderPath, studyName)
    if (dir.exists(folderPath)) {
      cwd = setwd(folderPath)
      eset = affy::justRMA(cdfname = platformInfo)
      setwd(cwd)
    } else {
      stop(sprintf('Folder %s does not exist.', folderPath))}

    Biobase::featureNames(eset) = fixCustomCdfGeneIds(Biobase::featureNames(eset))
    if (studyDataType == 'affy_geo') {
      Biobase::sampleNames(eset) = fixGeoSampleNames(Biobase::sampleNames(eset))
    } else {
      Biobase::sampleNames(eset) = fixCelSampleNames(Biobase::sampleNames(eset))}

  } else if (studyDataType == 'affy_series_matrix') {
    mapping = getGeneProbeMappingAffy(file.path(
      parentFolderPath, paste0(platformInfo, '_mapping.txt')))
    esetOrig = GEOquery::getGEO(filename = file.path(
      parentFolderPath, paste0(studyName, '_series_matrix.txt.gz')))
    exprs(esetOrig)[exprs(esetOrig) <= 0] = min(exprs(esetOrig)[exprs(esetOrig) > 0])
    exprs(esetOrig) = log2(exprs(esetOrig))
    exprsByGene = calcExprsByGene(esetOrig, mapping)
    rownames(exprsByGene) = fixCustomCdfGeneIds(rownames(exprsByGene))
    colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
    eset = Biobase::ExpressionSet(
      assayData = exprsByGene, phenoData = phenoData(esetOrig),
      experimentData = experimentData(esetOrig))

  } else if (studyDataType == 'series_matrix') {
    supportedPlatformsDt = getSupportedPlatforms()

    if (!(platformInfo %in% supportedPlatformsDt$platform)) {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.',
                      studyName, platformInfo))
      return(NA)}

    esetOrig = GEOquery::getGEO(filename = file.path(
      parentFolderPath, paste0(studyName, '_series_matrix.txt.gz')))
    if (is.list(esetOrig) && length(esetOrig) == 1) {
      esetOrig = esetOrig[[1]]}

    featureDf = pData(featureData(esetOrig))
    idx = sapply(featureDf, is.factor)
    featureDf[idx] = lapply(featureDf[idx], as.character)
    featureDt = data.table(featureDf)

    platformDt = supportedPlatformsDt[platform == platformInfo,]

    if (!(is.na(platformDt$splitColumn)) && platformDt$splitColumn != '') {
      featureDt[[platformDt$interName]] = sapply(
        featureDt[[platformDt$splitColumn]], function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
    }
    if (platformDt$mappingFunction == 'Anno') {
      mapping = getGeneProbeMappingAnno(
        featureDt, dbName = platformDt$dbName, interName = platformDt$interName)
    } else if (platformDt$mappingFunction == 'Direct')  {
      mapping = getGeneProbeMappingDirect(featureDt, geneColname = platformDt$geneColname)
    } else {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.',
                      studyName, platformInfo))
      return(NA)
    }

    exprsByGene = calcExprsByGene(esetOrig, mapping)
    if (any(is.na(exprsByGene))) {
      warning(sprintf('Imputing missing expression values for study %s.', studyName))
      resultImputed = impute::impute.knn(exprsByGene)
      exprsByGene = resultImputed$data}
    colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
    eset = Biobase::ExpressionSet(
      assayData = exprsByGene, phenoData = phenoData(esetOrig),
      experimentData = experimentData(esetOrig))

  } else if (studyDataType == 'eset_rds') {
    esetOrig = readRDS(file.path(parentFolderPath, paste0(studyName, '.rds')))

    featureDf = pData(featureData(esetOrig))
    featureDt = data.table(featureDf)
    if (platformInfo == 'ready') {
      return(esetOrig)
    } else if (platformInfo == 'rosetta') {
      mapping = getGeneProbeMappingDirect(
        featureDt, geneColname = 'EntrezGene.ID', probeColname = 'probe')
    } else {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.',
                      studyName, platformInfo))
      return(NA)}

    exprsByGene = calcExprsByGene(esetOrig, mapping)
    if (any(is.na(exprsByGene))) {
      warning(sprintf('Imputing missing expression values for study %s.', studyName))
      resultImputed = impute::impute.knn(exprsByGene)
      exprsByGene = resultImputed$data}
    eset = Biobase::ExpressionSet(
      assayData = exprsByGene, phenoData = phenoData(esetOrig),
      experimentData = experimentData(esetOrig))

  } else {
    warning(sprintf('Study %s not loaded, because data type %s is not currently supported.',
                    studyName, studyDataType))
    eset = NA}
  return(eset)}


#' Get the gene expression data for multiple studies.
#'
#' Runs [getStudyData()] for multiple studies.
#'
#' @param parentFolderPath Path to the folder that contains the data.
#' @param studyMetadata A `data.frame` where each row corresponds to one study.
#' 	 Should have columns for `study`, `studyDataType`, and `platformInfo`.
#'
#' @return A named list of `ExpressionSet`s.
#'
#' @export
getStudyDataList = function(parentFolderPath, studyMetadata) {
  # esetList = foreach(ii = 1:nrow(studyMetadata)) %do% {
  #   if (any(is.na(studyMetadata[ii,]))) {
  #     NA
  #   } else {
  #     getStudyData(parentFolderPath, studyMetadata$study[ii],
  #                  studyMetadata$studyDataType[ii],
  #                  studyMetadata$platformInfo[ii])}}
  esetList = foreach(studyRow = iterators::iter(studyMetadata, by = 'row')) %do% {
    if (any(is.na(studyRow))) {
      NA
    } else {
      getStudyData(parentFolderPath, studyRow$study, studyRow$studyDataType,
                   studyRow$platformInfo)}}
  names(esetList) = studyMetadata$study
  return(esetList[!is.na(esetList)])}


#' Extract the expression matrices containing the desired samples.
#'
#' @param esetList list of `ExpressionSets`.
#' @param sampleMetadata `data.frame` with columns for study and sample names.
#'
#' @return A named list of matrices.
#'
#' @export
extractExpressionData = function(esetList, sampleMetadata) {
  ematList = foreach(studyName = names(esetList)) %do% {
    # sampleNamesNow = dplyr::filter(sampleMetadata, study == studyName)$sample
    sampleNamesNow = data.table(sampleMetadata)[study == studyName]$sample
    keepIdx = colnames(esetList[[studyName]]) %in% sampleNamesNow
    exprs(esetList[[studyName]])[, keepIdx]}
  names(ematList) = names(esetList)
  return(ematList)}
