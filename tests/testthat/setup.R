library('ggplot2')
library('RColorBrewer')
library('data.table')
library('glmnet')
library('withr')

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

esetListControl = readRDS(file.path(parentFolderPath, 'esetList.rds'))

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
