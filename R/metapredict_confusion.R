#' Calculate confusion matrix for cross-validation.
#'
#' Calculate a confusion matrix based on predictions from cross-validation.
#'
#' @param cvFit cv.glmnet object from [metapredictCv()].
#' @param lambda value of lambda at which to use predictions.
#' @param ematMerged matrix of gene expression for genes by samples.
#' @param sampleMetadata data.frame of sample metadata.
#' @param className name of column in `sampleMetadata` containing the true
#'   labels.
#' @param classLevels Order of classes in the confusion matrix. If `NA`
#'   (default), then the function uses the order in `cvFit`.
#'
#' @return An object of class `table`.
#'
#' @export
calcConfusionCv = function(cvFit, lambda, ematMerged, sampleMetadata,
                           className='class', classLevels=NA) {
  if (is.na(classLevels[1])) {
    classLevels = names(cvFit$glmnet.fit$beta)}

  cvProbs = cvFit$fit.preval[,,which.min(abs(cvFit$lambda - lambda))]
  rownames(cvProbs) = colnames(ematMerged)
  colnames(cvProbs) = names(cvFit$glmnet.fit$beta)
  preds = colnames(cvProbs)[apply(cvProbs, MARGIN=1, function(x) which.max(x))]
  predictedClass = factor(preds, levels=classLevels)

  # classValues = tibble::tibble(sample = colnames(ematMerged)) %>%
  #   dplyr::inner_join(sampleMetadata) %>%
  #   .[[className]]

  classValues = merge(data.table(sample = colnames(ematMerged)), sampleMetadata)[[className]]
  trueClass = factor(classValues, levels=classLevels)
  return(table(trueClass, predictedClass))}


#' Calculate confusion matrices (or matrix) for validation datasets.
#'
#' Calculate confusion matrices based on predictions for validation datasets.
#'
#' @param predsList list of predictions from [metapredict()].
#' @param lambda value of lambda at which to use predictions.
#' @param sampleMetadata data.frame of sample metadata.
#' @param className name of column in `sampleMetadata` containing the true
#'   labels.
#' @param classLevels Order of classes in the confusion matrix. If `NA`
#'   (default), then the function uses the order in `cvFit`.
#' @param each logical indicating whether to calculate a confusion matrix for
#'   each validation dataset (default) or one confusion matrix including all
#'   datasets.
#'
#' @return If `each==TRUE`, a list of objects of class `table`. Otherwise, an
#'   object of class `table`.
#'
#' @export
calcConfusionValidation = function(predsList, lambda, sampleMetadata,
                                   className='class', classLevels=NA,
                                   each=TRUE) {
  if (is.na(classLevels[1])) {
    classLevels = colnames(predsList[[1]])}

  if (each) {
    confusion = list()
    for (validationStudyName in names(predsList)) {
      predsProb = predsList[[validationStudyName]][,,1]
      predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1,
                                             function(x) which.max(x))]
      predictedClass = factor(predsClass, levels=classLevels)

      # sm = tibble::tibble(sample = rownames(predsProb)) %>%
      #   dplyr::inner_join(sampleMetadata, by='sample')

      sm = data.table(sample = rownames(predsProb))[sampleMetadata, on = 'sample', nomatch = 0]
      trueClass = factor(sm[[className]], levels=classLevels)
      confusion[[validationStudyName]] = table(trueClass, predictedClass)}

  } else {
    predsProb = do.call(rbind, lapply(predsList, function(x) x[,,1]))
    predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1,
                                           function(x) which.max(x))]
    predictedClass = factor(predsClass, levels=classLevels)

    # sm = tibble::tibble(sample = rownames(predsProb)) %>%
    #   dplyr::inner_join(sampleMetadata, by='sample')

    sm = data.table(sample = rownames(predsProb))[sampleMetadata, on = 'sample', nomatch = 0]
    trueClass = factor(sm[[className]], levels=classLevels)
    confusion = table(trueClass, predictedClass)}

  return(confusion)}
