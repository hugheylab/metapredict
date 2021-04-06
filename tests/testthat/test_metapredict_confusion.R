library('ggplot2')
library('RColorBrewer')
library('data.table')
library('glmnet')

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

ematListControl = readRDS(file.path(parentFolderPath, 'ematList.rds'))

ematDiscoveryControl = readRDS(file.path(parentFolderPath, 'ematDiscovery.rds'))

glmnetArgsControl = readRDS(file.path(parentFolderPath, 'glmnetArgs.rds'))

cvFitListControl = readRDS(file.path(parentFolderPath, 'cvFitList.rds'))
cvFit = cvFitListControl[[1]]
fitResult = cvFit$glmnet.fit
lambda = cvFit$lambda.min

predsListControl = readRDS(file.path(parentFolderPath, 'predsList.rds'))

calcConfusionValidationControl = readRDS(file.path(parentFolderPath, 'calcConfusionValidation.rds'))

calcConfusionValidationEachFalseControl = readRDS(file.path(parentFolderPath, 'calcConfusionValidationEachFalse.rds'))

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
