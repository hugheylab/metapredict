#' Plot coefficients for a gene expression meta-analysis.
#'
#' Make a ggplot of the coefficients from a `glmnet` object trained on a gene
#' expression meta-analysis (either logistic or multinomial classification).
#'
#' @param fitResult glmnet object.
#' @param lambda value of lambda to use in `fitResult`.
#' @param classLevels Order of classes in the confusion matrix. If `NA` (default),
#'   then the function uses the order in `fitResult`.
#' @param decreasing logical indicating whether to sort genes by decreasing coefficient.
#' @param geneIdOrder Optional character array of Entrez Gene IDs specifying the order
#'   of genes. If `NA` (default), the order from [makeCoefDf()] is used.
#' @param org Name of package for mapping Entrez Gene IDs to gene symbols,
#'   passed to `data` argument of [annotate::lookUp()].
#'
#' @return A `ggplot` object.
#'
#' @export
plotCoefficients = function(fitResult, lambda, classLevels=NA, decreasing=FALSE,
                            geneIdOrder=NA, org='org.Hs.eg') {
  coefDf = makeCoefDf(fitResult, lambda, decreasing=decreasing, classLevels=classLevels)
  coefDf = coefDf[coefDf$geneId!='(Intercept)',]

  if (!is.na(geneIdOrder[1])) {
    rownames(coefDf) = coefDf$geneId
    coefDf = coefDf[geneIdOrder,]
    rownames(coefDf) = NULL}

  if (ncol(coefDf)==2) {
    geneSymbols = do.call(c, annotate::lookUp(coefDf$geneId, org, 'SYMBOL', load=TRUE))
    coefDf$geneId = factor(coefDf$geneId, levels=rev(coefDf$geneId),
                           labels=sprintf('%s (%s)', rev(geneSymbols), rev(coefDf$geneId)))
    p = ggplot(coefDf) + geom_bar(aes_string(x='geneId', y='coefficient'), stat='identity')

  } else {
    if (is.na(classLevels[1])) {
      classLevels = colnames(coefDf)[2:ncol(coefDf)]}

    # coefDfMolten = tidyr::gather(coefDf, key=class, value=coefficient, -geneId)
    coefDfMolten = melt(as.data.table(coefDf), id.vars = 'geneId', variable.name = 'class', variable.factor = FALSE, value.name = 'coefficient')
    coefDfMolten$class = factor(coefDfMolten$class, levels=classLevels)

    geneIds = coefDf$geneId
    geneSymbols = do.call(c, annotate::lookUp(coefDf$geneId, org, 'SYMBOL', load=TRUE))
    coefDfMolten$geneId = factor(coefDfMolten$geneId, levels=rev(geneIds),
                                 labels=sprintf('%s (%s)', rev(geneSymbols), rev(geneIds)))
    p = ggplot(coefDfMolten) + facet_wrap(~ class, ncol=ncol(coefDf)-1) +
      geom_bar(aes_string(x='geneId', y='coefficient', fill='class'), stat='identity') +
      guides(fill=FALSE)}

  return(p + coord_flip() + labs(x='Gene', y='Coefficient') + theme_light())}


#' Plot heatmap of merged gene expression data.
#'
#' Make a heatmap of gene expression from multiple datasets.
#'
#' @param fitResult `glmnet` object.
#' @param lambda value of lambda to use in `fitResult`.
#' @param ematMerged matrix of gene expression for genes by samples.
#' @param sampleMetadata data.frame of sample metadata.
#' @param annoLevels Named list used to make the `annotation`
#'   argument passed to [pheatmap::pheatmap()]. Each name must correspond
#'   to a column name in `sampleMetadata`, and each item in the
#'   list must be a vector of values found in that particular column.
#' @param annoColors Passed to `annotation_colors` argument of
#'   [pheatmap::pheatmap()].
#' @param clusterTogether logical indicating whether to cluster the samples
#'   from each dataset together or separately.
#' @param geneIdOrder Optional character array of Entrez Gene IDs specifying the order
#'   of genes. If `NA` (default), the order from [makeCoefDf()] is used.
#' @param className column in sampleMetadata containing values of the response variable.
#' @param classLevels Order of classes for the column annotations.
#' @param org Name of package for mapping Entrez Gene IDs to gene symbols,
#'   passed to `data` argument of [annotate::lookUp()].
#' @param maxVal Maximum absolute value of scaled and centered gene expression, used
#'   to control dynamic range of color in the heatmap.
#' @param ... Additional arguments passed to [pheatmap::pheatmap()].
#'
#' @return A `pheatmap` object.
#'
#' @export
plotExpressionHeatmap = function(fitResult, lambda, ematMerged, sampleMetadata, annoLevels,
                                 annoColors=NA, clusterTogether=FALSE, geneIdOrder=NA,
                                 className='class', classLevels=NA, org='org.Hs.eg',
                                 maxVal=3, ...) {
  coefDf = makeCoefDf(fitResult, lambda)
  geneIds = coefDf$geneId[coefDf$geneId!='(Intercept)']
  geneSymbols = do.call(c, annotate::lookUp(geneIds, org, 'SYMBOL', load=TRUE))
  geneTexts = sprintf('%s (%s)', geneSymbols, geneIds)
  names(geneTexts) = geneIds
  emat = ematMerged[geneIds,]

  # order the samples
  if (clusterTogether) {
    d = stats::dist(t(emat))
    co = cba::order.optimal(d, stats::hclust(d)$merge)
    emat = emat[,co$order]
  } else {
    if (is.na(classLevels[1])) {
      # sm = dplyr::filter(sampleMetadata, sample %in% colnames(ematMerged))
      sm = data.table(sampleMetadata)[which(sample %in% colnames(ematMerged)),]
      classLevels = unique(sm[[className]])}

    ematSmallList = foreach(classLevel=classLevels) %do% {
      # sampleNames1 = sampleMetadata %>%
      #   dplyr::filter_(lazyeval::interp(~ a==classLevel, a=as.name(className))) %>%
      #   .$sample
      sampleNames1 = sampleMetadata$sample[sampleMetadata[[className]] == classLevel]
      x = emat[, colnames(emat) %in% sampleNames1]
      d = stats::dist(t(x))
      co = cba::order.optimal(d, stats::hclust(d)$merge)
      x = x[, co$order]}
    emat = do.call(cbind, ematSmallList)}

  # order the genes
  if (is.na(geneIdOrder[1])) {
    d = stats::dist(emat)
    co = cba::order.optimal(d, stats::hclust(d)$merge)
    emat = emat[co$order,]
    rownames(emat) = geneTexts[co$order]
  } else {
    emat = emat[geneIdOrder,]
    rownames(emat) = geneTexts[geneIdOrder]}

  # scale the matrix
  emat = t(scale(t(emat)))
  emat[emat > maxVal] = maxVal
  emat[emat < (-maxVal)] = -maxVal

  # annotation = tibble::tibble(sample = colnames(ematMerged)) %>%
  #   dplyr::inner_join(sampleMetadata, by='sample') %>%
  #   dplyr::select(!!names(annoLevels)) %>%
  #   as.data.frame()
  cols = names(annoLevels)
  annotation = as.data.frame(merge(data.table(sample = colnames(ematMerged)), sampleMetadata, by = 'sample', sort = FALSE)[,..cols])
  rownames(annotation) = colnames(ematMerged)

  for (annoName in names(annoLevels)) {
    if (!is.na(annoLevels[[annoName]][1])) {
      annotation[[annoName]] = factor(annotation[[annoName]], levels=annoLevels[[annoName]])}}

  p = pheatmap::pheatmap(emat, color=grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name='RdBu')))(100),
                         breaks=seq(-maxVal, maxVal, length.out=101), cluster_rows=FALSE, cluster_cols=FALSE,
                         treeheight_row=0, treeheight_col=0, show_colnames=FALSE, border_color=NA,
                         annotation_col=annotation, annotation_colors=annoColors, ...)
  return(p)}


#' Plot class probabilities for samples from cross-validation.
#'
#' Make a gtable consisting of multiple `ggplot`s, one for each dataset used in
#' cross-validation. Within each class, samples are sorted by the probability of
#' the true class.
#'
#' @param cvFit `cv.glmnet` object from [metapredictCv()].
#' @param lambda value of lambda to use in `cvFit`.
#' @param ematMerged matrix of gene expression for genes by samples.
#' @param sampleMetadata data.frame of sample metadata.
#' @param className column in sampleMetadata containing values of the response
#'   variable.
#' @param classLevels Preferred order of classes.
#' @param size Point size, passed to [ggplot2::geom_point()].
#' @param ggplotArgs List of additional arguments to add to each ggplot.
#'
#' @return A `ggplot` object constructed by [cowplot::plot_grid()].
#'
#' @export
plotClassProbsCv = function(cvFit, lambda, ematMerged, sampleMetadata, className='class', classLevels=NA,
                            size=1.5, ggplotArgs=NA) {

  sampleNames = colnames(ematMerged)
  # sm = tibble::tibble(sample = sampleNames) %>%
  #   dplyr::inner_join(sampleMetadata, by='sample')
  sm = merge(data.table(sample = sampleNames), sampleMetadata, by = 'sample', sort = FALSE)
  studyNames = sort(unique(sm$study))

  if (is.na(classLevels[1])) {
    classLevels = sort(unique(sm[[className]]))}

  cvProbs = cvFit$fit.preval[,,which.min(abs(cvFit$lambda - lambda))]
  pList = list()
  for (studyName in studyNames) {
    # sampleNamesNow = dplyr::filter(sm, study==studyName)$sample
    sampleNamesNow = sm[which(study == studyName), sample]
    # df = tibble::as_tibble(cvProbs[sm$study==studyName,])
    df = as.data.frame(cvProbs[sm$study==studyName,])
    colnames(df) = names(cvFit$glmnet.fit$beta)
    df$study = studyName
    df$sample = sampleNamesNow
    dfDT = df
    # df$trueClass = factor(dplyr::filter(sm, sample %in% sampleNamesNow)[[className]], levels=classLevels)
    df$trueClass = factor(data.table(sm)[which(sample %in% sampleNamesNow),][[className]], levels=classLevels)

    df$trueClassProb = apply(df, MARGIN=1, function(x) as.numeric(x[x['trueClass']]))

    df = df[order(df$trueClass, -df$trueClassProb),]
    df = do.call(rbind, lapply(classLevels, function(x) df[df$trueClass==x,]))

    idxTmp = c()
    for (classLevel in classLevels) {
      if (any(df$trueClass==classLevel)) {
        idxTmp = c(idxTmp, 1:(sum(df$trueClass==classLevel)))}}
    df$idx = idxTmp
    # dfMolten = tidyr::gather(df, key='probClass', value='prob', !!classLevels)
    dfMolten = melt(as.data.table(df), measure.vars = classLevels, variable.name = 'probClass', variable.factor = FALSE, value.name = 'prob')

    p = ggplot(dfMolten) +
      facet_grid(study ~ trueClass, scales='free_x', space='free_x') +
      geom_point(aes_string(x='idx', y='prob', color='probClass', shape='probClass'), size=size) +
      labs(x='Sample', y='Probability') + theme_light() + theme(legend.title=element_blank())

    if (!is.na(ggplotArgs[1])) {
      for (ggplotArg in ggplotArgs) {
        p = p + ggplotArg}}
    pList[[studyName]] = p}

  # g = do.call(gridExtra::arrangeGrob, c(pList, list(nrow=length(studyNames))))
  p = cowplot::plot_grid(plotlist = pList, ncol = 1, align = 'h')
  return(p)}


#' Plot class probabilities for samples from multiple validation datasets.
#'
#' Make a gtable consisting of multiple ggplots, one for each validation
#' dataset. Within each class, samples are sorted by the probability of the true
#' class.
#'
#' @param predsList Result from [metapredict()].
#' @param sampleMetadata data.frame of sample metadata.
#' @param className column in sampleMetadata containing values of the response
#'   variable.
#' @param classLevels Preferred order of classes.
#' @param size Point size, passed to [ggplot2::geom_point()].
#' @param ggplotArgs List of additional arguments to add to each `ggplot`.
#'
#' @return A ggplot object constructed by [cowplot::plot_grid()].
#'
#' @export
plotClassProbsValidation = function(predsList, sampleMetadata, className,
                                    classLevels, size=1.5, ggplotArgs=NA) {
  pList = list()
  for (validationStudyName in names(predsList)) {
    df = data.frame(predsList[[validationStudyName]][,,1])
    # sm = tibble::tibble(sample = rownames(df)) %>%
    #   dplyr::inner_join(sampleMetadata, by='sample')

    sm = merge(data.table(sample = rownames(df)), sampleMetadata, by = 'sample', sort = FALSE)

    df$study = sm[,study]
    df$sample = rownames(df)
    df$trueClass = factor(sm[[className]], levels=classLevels)
    df$trueClassProb = apply(df, MARGIN=1, function(x) as.numeric(x[x['trueClass']]))

    df = df[order(df$trueClass, -df$trueClassProb),]
    df = do.call(rbind, lapply(classLevels, function(x) df[df$trueClass==x,]))

    idxTmp = c()
    for (classLevel in classLevels) {
      if (any(df$trueClass==classLevel)) {
        idxTmp = c(idxTmp, 1:(sum(df$trueClass==classLevel)))}}
    df$idx = idxTmp
    rownames(df) = NULL

    # dfMolten = tidyr::gather(df, key='probClass', value='prob', !!classLevels)
    dfMolten = melt(as.data.table(df), measure.vars = classLevels, variable.name = 'probClass', variable.factor = FALSE, value.name = 'prob')

    p = ggplot(dfMolten) +
      facet_grid(study ~ trueClass, scales='free_x', space='free_x') +
      geom_point(aes_string(x='idx', y='prob', color='probClass', shape='probClass'), size=size) +
      labs(x='Sample', y='Probability') + theme_light() + theme(legend.title=element_blank())
    if (!is.na(ggplotArgs[1])) {
      for (ggplotArg in ggplotArgs) {
        p = p + ggplotArg}}
    pList[[validationStudyName]] = p}

  # g = do.call(gridExtra::arrangeGrob, c(pList, list(nrow=length(predsList))))
  p = cowplot::plot_grid(plotlist = pList, ncol = 1, align = 'h')
  return(p)}
