library('data.table')
library('ggplot2')
library('RColorBrewer')
library('withr')

local_file('RPlots.pdf')

parentFolderPath = 'test_data'
foldidColname = 'study'
className = 'class'
family = 'multinomial'
alpha = 0.9
discoveryStudyNames = c('Bhattacharjee', 'GSE11969', 'GSE29016')
classLevels = c('AD', 'SQ', 'SCLC')

studyMetadataPath = file.path(parentFolderPath, 'study_metadata.csv')
studyMetadata = read.csv(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path(parentFolderPath, 'sample_metadata.csv')
sampleMetadata = read.csv(sampleMetadataPath, stringsAsFactors = FALSE)

ematDiscoveryControl = readRDS(file.path(parentFolderPath, 'ematDiscovery.rds'))

glmnetArgsControl = readRDS(file.path(parentFolderPath, 'glmnetArgs.rds'))

cvFitListControl = readRDS(file.path(parentFolderPath, 'cvFitList.rds'))
cvFit = cvFitListControl[[1]]
fitResult = cvFit$glmnet.fit
lambda = cvFit$lambda.min

predsListControl = readRDS(file.path(parentFolderPath, 'predsList.rds'))

plotCoefficientsControl = readRDS(file.path(parentFolderPath, 'plotCoefficients.rds'))

annoNames = c('study', 'class')
annoLevels = list(discoveryStudyNames, classLevels)
names(annoLevels) = annoNames
annoColors = list(brewer.pal(3, 'Dark2'), brewer.pal(3, 'Paired'))
names(annoColors) = annoNames
names(annoColors[[1]]) = annoLevels[[1]]
names(annoColors[[2]]) = annoLevels[[2]]

plotExpressionHeatmapControl = readRDS(file.path(parentFolderPath, 'plotExpressionHeatmap.rds'))

plotClassProbsValidationControl = readRDS(file.path(parentFolderPath, 'plotClassProbsValidation.rds'))

test_that('plotCoefficients', {
  plotCoefficientsTest = plotCoefficients(fitResult, lambda, classLevels = classLevels) +
    scale_fill_brewer(type = 'qual', palette = 3) +
    theme(axis.text.y = element_text(size = 8))
  expect_true(all.equal(plotCoefficientsTest$data, plotCoefficientsControl$data, tolerance = 0.000001, check.attributes = FALSE))
})

test_that('plotExpressionHeatmap', {
  plotExpressionHeatmapTest = plotExpressionHeatmap(fitResult, lambda, ematDiscoveryControl, sampleMetadata, annoLevels,
                                                annoColors, classLevels = classLevels, fontsize_row = 7)
  expect_true(all.equal(plotExpressionHeatmapTest$data, plotExpressionHeatmapControl$data, check.attributes = FALSE))
})

test_that('plotClassProbsValidation', {
  plotClassProbsValidationTest = plotClassProbsValidation(predsListControl, sampleMetadata, className, classLevels, size = 1.5,
                                                      ggplotArgs = list(scale_color_brewer(type = 'qual', palette = 3)))
  expect_true(all.equal(plotClassProbsValidationTest$data, plotClassProbsValidationControl$data, check.attributes = FALSE))
})
