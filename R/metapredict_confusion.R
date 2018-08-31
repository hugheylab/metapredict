#' Calculate confusion matrix for cross-validation.
#'
#'\code{calcConfusionCv} calculates a confusion matrix of cross-validation
#' results from \code{metapredictCv}.
#'
#' @param cvFit cv.glmnet object from \code{metapredictCv}.
#' @param lambda value of lambda at which to use predictions.
#' @param ematMerged matrix of gene expression for genes by samples.
#' @param sampleMetadata data.frame of sample metadata.
#' @param className name of column in \code{sampleMetadata} containing the true
#' labels.
#' @param classLevels Order of classes in the confusion matrix. If \code{NA}
#' (default), then the function uses the order in \code{cvFit}.
#'
#' @return An object of class \code{table}.
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
  classValues = tibble(sample = colnames(ematMerged)) %>%
    inner_join(sampleMetadata) %>%
    .[[className]]
  trueClass = factor(classValues, levels=classLevels)
  return(table(trueClass, predictedClass))}


#' Calculate confusion matrices (or matrix) for validation datasets.
#'
#' \code{calcConfusionCv} calculates confusion matrices based on predictions
#' for validation datasets from \code{metapredict}.
#'
#' @param predsList list of predictions from \code{metapredict}.
#' @param lambda value of lambda at which to use predictions.
#' @param sampleMetadata data.frame of sample metadata.
#' @param className name of column in \code{sampleMetadata} containing the true
#' labels.
#' @param classLevels Order of classes in the confusion matrix. If \code{NA}
#' (default), then the function uses the order in \code{cvFit}.
#' @param each logical indicating whether to calculate a confusion matrix for
#' each validation dataset (default) or one confusion matrix including all
#' datasets.
#'
#' @return If \code{each==TRUE}, a list of objects of class \code{table}.
#' Otherwise, an object of class \code{table}.
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
      sm = tibble::tibble(sample = rownames(predsProb)) %>%
        dplyr::inner_join(sampleMetadata, by='sample')
      trueClass = factor(sm[[className]], levels=classLevels)
      confusion[[validationStudyName]] = table(trueClass, predictedClass)}

  } else {
    predsProb = do.call(rbind, lapply(predsList, function(x) x[,,1]))
    predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1,
                                           function(x) which.max(x))]
    predictedClass = factor(predsClass, levels=classLevels)
    sm = tibble::tibble(sample = rownames(predsProb)) %>%
      dplyr::inner_join(sampleMetadata, by='sample')
    trueClass = factor(sm[[className]], levels=classLevels)
    confusion = table(trueClass, predictedClass)}

  return(confusion)}
