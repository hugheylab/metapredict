test_that('makeGlmnetArgs', {
  local_seed(15)
  glmnetArgsTest = makeGlmnetArgs(
    sampleMetadata[sampleMetadata$study %in% discoveryStudyNames, ],
    foldidColname = foldidColname)
  expect_equal(glmnetArgsTest, glmnetArgsControl, check.attributes = FALSE)
})

test_that('metapredictCv', {
  local_seed(31)
  cvFitListTest = metapredictCv(
    ematDiscoveryControl, sampleMetadata, yName = className,
    weights = glmnetArgsControl$weights, foldid = glmnetArgsControl$foldid,
    alpha = alpha, family = family, lambda.min.ratio = 0.001, keep = TRUE)
  expect_equal(cvFitListTest, cvFitListControl, check.attributes = FALSE)
})

test_that('metapredict', {
  local_seed(63)
  predsListTest = metapredict(
    ematListControl, studyMetadata, sampleMetadata, discoveryStudyNames,
    alpha = alpha, lambda = lambda, weights = glmnetArgsControl$weights,
    family = family)
  expect_equal(predsListTest, predsListControl, check.attributes = FALSE)
})
