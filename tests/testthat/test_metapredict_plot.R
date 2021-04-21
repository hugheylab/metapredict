

test_that('plotCoefficients', {
  plotCoefficientsTest = plotCoefficients(fitResult, lambda, classLevels = classLevels) +
    scale_fill_brewer(type = 'qual', palette = 3) +
    theme(axis.text.y = element_text(size = 8))
  expect_equal(plotCoefficientsTest$data, plotCoefficientsControl$data, check.attributes = FALSE)
})

test_that('plotExpressionHeatmap', {
  plotExpressionHeatmapTest = plotExpressionHeatmap(fitResult, lambda, ematDiscoveryControl, sampleMetadata, annoLevels,
                                                annoColors, classLevels = classLevels, fontsize_row = 7)
  expect_equal(plotExpressionHeatmapTest$data, plotExpressionHeatmapControl$data, check.attributes = FALSE)
})

test_that('plotClassProbsValidation', {
  plotClassProbsValidationTest = plotClassProbsValidation(predsListControl, sampleMetadata, className, classLevels, size = 1.5,
                                                      ggplotArgs = list(scale_color_brewer(type = 'qual', palette = 3)))
  expect_equal(plotClassProbsValidationTest$data, plotClassProbsValidationControl$data, check.attributes = FALSE)
})
