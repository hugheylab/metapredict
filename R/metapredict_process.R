#' @import ggplot2
#' @import methods
#' @import data.table
#' @importFrom foreach foreach
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importMethodsFrom AnnotationDbi mappedkeys
#' @importMethodsFrom Biobase experimentData
#' @importMethodsFrom Biobase exprs
#' @importMethodsFrom Biobase exprs<-
#' @importMethodsFrom Biobase featureData
#' @importMethodsFrom Biobase pData
#' @importMethodsFrom Biobase phenoData
#' @importMethodsFrom BiocGenerics as.list
NULL

globalVariables(c('studyName', 'coefSparse', 'validationStudyName', 'classLevel',
                  '.', 'study', 'ii', 'studyDataType', 'platformInfo', 'geneId', 'geneIdTmp',
                  'coefficient'))


#' Install custom CDF packages from Brainarray.
#'
#' Install Brainarray custom CDFs for processing raw Affymetrix data. See
#' <http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp>.
#'
#' @param pkgs character vector of package names, e.g., 'hgu133ahsentrezgcdf'
#' @param ver integer version number (25 as of 5 Jan 2021)
#'
#' @export
installCustomCdfPackages = function(pkgs, ver = 25) {
  for (pkg in pkgs) {
    pkgUrl = sprintf('http://mbni.org/customcdf/%d.0.0/entrezg.download/%s_%d.0.0.tar.gz',
                     ver, pkg, ver)
    utils::install.packages(pkgUrl, repos = NULL)}}


#' Download custom CDF mapping files from Brainarray.
#'
#' Download Brainarray custom CDF mapping files, which are used for mapping probes to genes
#' in datasets whose `studyDataType` is 'affy_series_matrix'. See
#' \url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp}.
#'
#' @param cdf `data frame` with columns `download` (e.g., 'Mouse4302_Mm_ENTREZ')
#'   and `rename` (e.g., 'mouse4302mmentrezgcdf').
#' @param path directory into which to download the files.
#' @param ver integer version number (25 as of 5 Jan 2021).
#'
#' @export
downloadCustomCdfMappings = function(cdf, path = '.', ver = 25) {
  if (!dir.exists(path)) {
    dir.create(path)}
  for (ii in 1:nrow(cdf)) {
    temp = tempfile()
    utils::download.file(sprintf('http://mbni.org/customcdf/%d.0.0/entrezg.download/%s_%d.0.0.zip',
                                 ver, cdf$download[ii], ver), temp)
    utils::unzip(temp, files = paste0(cdf$download[ii], '_mapping.txt'), exdir = path)
    file.rename(file.path(path, paste0(cdf$download[ii], '_mapping.txt')),
                file.path(path, paste0(cdf$rename[ii], '_mapping.txt')))
    unlink(temp)}}


fixCustomCdfGeneIds = function(geneIds) {
  return(sub('_at', '', geneIds))}


fixGeoSampleNames = function(sampleNames) {
  sampleNames = paste0(toupper(sampleNames), '_')
  regexResult = regexpr('^GSM[0-9]+[^0-9]', sampleNames)
  sampleNamesNew = mapply(function(sampleName, matchLength) substr(sampleName, 1, matchLength-1),
                          sampleNames, attr(regexResult, 'match.length'))
  return(unname(sampleNamesNew))}


fixCelSampleNames = function(sampleNames) {
  sampleNamesNew = gsub('\\.cel(\\.gz)*$', '', sampleNames, ignore.case = TRUE)
  return(sampleNamesNew)}


getGeneProbeMappingAffy = function(mappingFilePath) {
  mapping = utils::read.table(mappingFilePath, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  mappingUnique = unique(mapping[,c('Probe.Set.Name', 'Affy.Probe.Set.Name')])
  mappingUnique = mappingUnique[apply(mappingUnique, MARGIN = 1, function(r) !any(is.na(r))),]
  rownames(mappingUnique) = NULL
  colnames(mappingUnique) = c('geneId', 'probeSet')
  mappingUnique[,'probeSet'] = as.character(mappingUnique[,'probeSet'])
  return(mappingUnique)}


getGeneProbeMappingDirect = function(featureDf, geneColname, probeColname = 'ID') {
  mapping = featureDf[,c(..probeColname, ..geneColname)]
  mapping = mapping[apply(mapping, MARGIN = 1, function(x) all(!is.na(x) & x!='')),]
  mapping = data.frame(lapply(mapping, as.character), stringsAsFactors = FALSE)
  colnames(mapping) = c('probeSet', 'geneId')
  return(mapping)}


getGeneProbeMappingAnno = function(featureDf, dbName, interName) {
  # mappingProbeIntermediate = featureDf[!is.na(featureDf[,interName]) & featureDf[,interName]!='',
  #                                      c('ID', interName)]
  mappingProbeIntermediate = featureDf[!is.na(featureDf[[interName]]) & featureDf[[interName]]!='',
                                       c('ID', ..interName)]
  colnames(mappingProbeIntermediate) = c('probeSet', 'geneInter')

  mapTmp1 = eval(parse(text = sprintf('%s.db::%s', substr(dbName, 1, 9), dbName)))
  mapTmp2 = mappedkeys(mapTmp1)
  mapTmp3 = as.list(mapTmp1[mapTmp2])
  geneId = do.call(c, mapTmp3)

  geneInter = do.call(c, mapply(function(inter, len) rep_len(inter, len), names(mapTmp3),
                                sapply(mapTmp3, length), SIMPLIFY = FALSE))
  if (dbName == 'org.Hs.egUNIGENE2EG') {
    geneInter = sub('Hs.', '', geneInter, fixed = TRUE)}
  mappingIdInter = data.table(geneId, geneInter, stringsAsFactors = FALSE)
  mapping = merge(mappingIdInter, mappingProbeIntermediate, by = 'geneInter', sort = FALSE)
  mapping$probeSet = as.character(mapping$probeSet)
  return(mapping)}


calcExprsByGene = function(eset, mapping) {
  geneIds = unique(mapping[['geneId']])
  exprsByGene = matrix(nrow = length(geneIds), ncol = ncol(eset),
                       dimnames = list(geneIds, Biobase::sampleNames(eset)))
  foreach(geneIdTmp = geneIds) %do% {
    # exprsTmp = exprs(eset)[mappingDf[mappingDf[,'geneId'] == geneId, 'probeSet'],, drop = FALSE]
    exprsTmp = exprs(eset)[mapping[geneId == geneIdTmp, probeSet],, drop = FALSE]
    if (nrow(exprsTmp) == 1) {
      exprsByGene[geneIdTmp,] = exprsTmp
    } else {
      exprsByGene[geneIdTmp,] = matrixStats::rowMedians(t(exprsTmp), na.rm = TRUE)}}
  return(exprsByGene)}


#' Get the GPLs for microarray platforms that are currently supported.
#'
#' @return A character array of currently supported GPLs.
#'
#' @export
getSupportedPlatforms = function() {
  return(c('GPL180', 'GPL341', 'GPL571', 'GPL885', 'GPL887', 'GPL890', 'GPL962',
           'GPL1053', 'GPL1073', 'GPL1261', 'GPL1291', 'GPL1293', 'GPL1390',
           'GPL1708', 'GPL3921', 'GPL4133', 'GPL4372', 'GPL5645', 'GPL6104',
           'GPL6254', 'GPL6333', 'GPL6480', 'GPL6865', 'GPL6880', 'GPL6884',
           'GPL6885', 'GPL6887', 'GPL6947', 'GPL7015', 'GPL7202', 'GPL8177',
           'GPL10332', 'GPL10379', 'GPL10558', 'GPL10687', 'GPL13607',
           'GPL13730', 'GPL15331', 'GPL15450', 'GPL18721', 'GPL20769'))}


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

  unsupportedPlatforms = data.table(studyMetadata)[which(studyDataType == 'series_matrix' & !(platformInfo %in% getSupportedPlatforms())),platformInfo]

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
    mapping = getGeneProbeMappingAffy(file.path(parentFolderPath,
                                                paste0(platformInfo, '_mapping.txt')))
    esetOrig = GEOquery::getGEO(filename = file.path(parentFolderPath,
                                                   paste0(studyName, '_series_matrix.txt.gz')))
    exprs(esetOrig)[exprs(esetOrig)<= 0] = min(exprs(esetOrig)[exprs(esetOrig)>0])
    exprs(esetOrig) = log2(exprs(esetOrig))
    exprsByGene = calcExprsByGene(esetOrig, mapping)
    rownames(exprsByGene) = fixCustomCdfGeneIds(rownames(exprsByGene))
    colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
    eset = Biobase::ExpressionSet(assayData = exprsByGene, phenoData = phenoData(esetOrig),
                                  experimentData = experimentData(esetOrig))

  } else if (studyDataType == 'series_matrix') {
    supportedPlatforms = getSupportedPlatforms()

    if (!(platformInfo %in% supportedPlatforms)) {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.',
                      studyName, platformInfo))
      return(NA)}

    esetOrig = GEOquery::getGEO(filename = file.path(parentFolderPath,
                                                   paste0(studyName, '_series_matrix.txt.gz')))
    if (is.list(esetOrig) && length(esetOrig) == 1) {
      esetOrig = esetOrig[[1]]}

    # featureDf = pData(featureData(esetOrig))
    # idx = sapply(featureDf, is.factor)
    # featureDf[idx] = lapply(featureDf[idx], as.character)
    featureDf = setDT(pData(featureData(esetOrig)))
    changeCols = colnames(featureDf)[which(as.vector(featureDf[,lapply(.SD, class)]) == "factor")]
    featureDf[,(changeCols):= lapply(.SD, as.character), .SDcols = changeCols]
    if (platformInfo == 'GPL180') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egSYMBOL2EG',
                                        interName = 'GENE_SYM')
    } else if (platformInfo == 'GPL341') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'ENTREZ_GENE_ID')
    } else if (platformInfo == 'GPL571') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'ENTREZ_GENE_ID')
    } else if (platformInfo == 'GPL885') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'GENE')
    } else if (platformInfo == 'GPL887') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'GENE')
    } else if (platformInfo == 'GPL890') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'GENE')
    } else if (platformInfo == 'GPL962') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egUNIGENE2EG',
                                        interName = 'UNIGENE')
    } else if (platformInfo == 'GPL1053') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egSYMBOL2EG',
                                        interName = 'GENE')
    } else if (platformInfo == 'GPL1073') {
      featureDf[,'GenBank'] = sapply(featureDf[,'GB_ACC'],
                                     function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Mm.egACCNUM2EG',
                                        interName = 'GenBank')
    } else if (platformInfo == 'GPL1261') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'ENTREZ_GENE_ID')
    } else if (platformInfo == 'GPL1291') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'Entrez Gene ID')
    } else if (platformInfo == 'GPL1293') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'Entrez Gene ID')
    } else if (platformInfo == 'GPL1390') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'GB_ACC')
    } else if (platformInfo == 'GPL1708') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'GENE')
    } else if (platformInfo == 'GPL3921') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'ENTREZ_GENE_ID')
    } else if (platformInfo == 'GPL4133') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'GENE')
    } else if (platformInfo == 'GPL4372') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'EntrezGeneID')
    } else if (platformInfo == 'GPL5645') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egSYMBOL2EG',
                                        interName = 'Gene Name')
    } else if (platformInfo == 'GPL6104') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'Entrez_Gene_ID')
    } else if (platformInfo == 'GPL6254') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egENSEMBL2EG',
                                        interName = 'ENSEMBL_GENE_ID')
    } else if (platformInfo == 'GPL6333') {
      featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'],
                                    function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Mm.egREFSEQ2EG',
                                        interName = 'RefSeq')
    } else if (platformInfo == 'GPL6480') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'GENE')
    } else if (platformInfo == 'GPL6865') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'rep_name')
    } else if (platformInfo == 'GPL6880') {
      featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'],
                                    function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Mm.egREFSEQ2EG',
                                        interName = 'RefSeq')
    } else if (platformInfo == 'GPL6884') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'Entrez_Gene_ID')
    } else if (platformInfo == 'GPL6885') {
      featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'],
                                    function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Mm.egREFSEQ2EG',
                                        interName = 'RefSeq')
    } else if (platformInfo == 'GPL6887') {
      featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'],
                                    function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Mm.egREFSEQ2EG',
                                        interName = 'RefSeq')
    } else if (platformInfo == 'GPL6947') {
      featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'],
                                    function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'RefSeq')
    } else if (platformInfo == 'GPL7015') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'GB_LIST')
    } else if (platformInfo == 'GPL7202') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'GENE')
    } else if (platformInfo == 'GPL8177') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'GB_ACC')
    } else if (platformInfo == 'GPL10332') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'GB_ACC')
    } else if (platformInfo == 'GPL10379') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'EntrezGeneID')
    } else if (platformInfo == 'GPL10558') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'Entrez_Gene_ID')
    } else if (platformInfo == 'GPL10687') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'EntrezGeneID')
    } else if (platformInfo == 'GPL13607') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'GB_ACC')
    } else if (platformInfo == 'GPL13730') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'ENTREZ_GENE_ID')
    } else if (platformInfo == 'GPL15331') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'GB_ACC')
    } else if (platformInfo == 'GPL15450') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Dr.egREFSEQ2EG',
                                        interName = 'GB_ACC')
    } else if (platformInfo == 'GPL18721') {
      featureDf[,'RefSeq'] = sapply(featureDf[,'GB_ACC'],
                                    function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egREFSEQ2EG',
                                        interName = 'RefSeq')
    } else if (platformInfo == 'GPL20769') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName = 'org.Hs.egSYMBOL2EG',
                                        interName = 'ORF')
    } else {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.',
                      studyName, platformInfo))
      return(NA)}

    exprsByGene = calcExprsByGene(esetOrig, mapping)
    if (any(is.na(exprsByGene))) {
      warning(sprintf('Imputing missing expression values for study %s.', studyName))
      resultImputed = impute::impute.knn(exprsByGene)
      exprsByGene = resultImputed$data}
    colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
    eset = Biobase::ExpressionSet(assayData = exprsByGene, phenoData = phenoData(esetOrig),
                                  experimentData = experimentData(esetOrig))

  } else if (studyDataType == 'eset_rds') {
    esetOrig = readRDS(file.path(parentFolderPath, paste0(studyName, '.rds')))

    # featureDf = pData(featureData(esetOrig))
    featureDf = setDT(pData(featureData(esetOrig)))
    if (platformInfo == 'ready') {
      return(esetOrig)
    } else if (platformInfo == 'rosetta') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname = 'EntrezGene.ID',
                                          probeColname = 'probe')
    } else {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.',
                      studyName, platformInfo))
      return(NA)}

    exprsByGene = calcExprsByGene(esetOrig, mapping)
    if (any(is.na(exprsByGene))) {
      warning(sprintf('Imputing missing expression values for study %s.', studyName))
      resultImputed = impute::impute.knn(exprsByGene)
      exprsByGene = resultImputed$data}
    eset = Biobase::ExpressionSet(assayData = exprsByGene, phenoData = phenoData(esetOrig),
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
  esetList = foreach(studyRow = iterators::iter(studyMetadata, by = 'row'), .combine = c) %do% {
    if (any(is.na(studyRow))) {
      NA
    } else {
      getStudyData(parentFolderPath, studyRow$study,
                   studyRow$studyDataType,
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
    sampleNamesNow = data.table(sampleMetadata)[study == studyName,]$sample
    keepIdx = colnames(esetList[[studyName]]) %in% sampleNamesNow
    exprs(esetList[[studyName]])[,keepIdx]}
  names(ematList) = names(esetList)
  return(ematList)}
