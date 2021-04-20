

test_that('calcConfusionValidation each = TRUE', {
  calcConfusionValidationTest = calcConfusionValidation(predsListControl, lambda = lambda, sampleMetadata = sampleMetadata,
                                          classLevels = classLevels)
  expect_equal(calcConfusionValidationTest, calcConfusionValidationControl, check.attributes = FALSE)
})

test_that('calcConfusionValidation each = FALSE', {
  calcConfusionValidationEachFalseTest = calcConfusionValidation(predsListControl, lambda = lambda, sampleMetadata = sampleMetadata,
                                                        classLevels = classLevels, each = FALSE)
  expect_equal(calcConfusionValidationEachFalseTest, calcConfusionValidationEachFalseControl, check.attributes = FALSE)
})
