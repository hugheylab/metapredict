library('metapredict')
library('ggplot2')
library('RColorBrewer')
# library('foreach')
# library('data.table')
# library('dplyr')

if (!(grepl('vignettes', getwd()))) {
  setwd('vignettes')
}

parentFolderPath = 'expression_data'
expressionDataFile = 'expression_data.rds'
denovo = TRUE

foldidColname = 'study'
className = 'class'
family = 'multinomial'
alpha = 0.9
discoveryStudyNames = c('Bhattacharjee', 'GSE11969', 'GSE29016')
classLevels = c('AD', 'SQ', 'SCLC')

studyMetadataPath = system.file('extdata', 'study_metadata.csv', package = 'metapredict')
studyMetadata = read.csv(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = system.file('extdata', 'sample_metadata.csv', package = 'metapredict')
sampleMetadata = read.csv(sampleMetadataPath, stringsAsFactors = FALSE)

if (denovo) {
  esetListOrig = getStudyDataList(parentFolderPath, studyMetadata)
  saveRDS(esetList, file = expressionDataFile)
} else {
  esetListOrig = readRDS(expressionDataFile)}

ematListOrig = extractExpressionData(esetListOrig, sampleMetadata)
ematDiscoveryOrig = mergeStudyData(ematListOrig[discoveryStudyNames], sampleMetadata)

glmnetArgs = makeGlmnetArgs(sampleMetadata[sampleMetadata$study %in% discoveryStudyNames,],
                            foldidColname = foldidColname)

cvFitList = metapredictCv(ematDiscovery, sampleMetadata, yName = className,
                          weights = glmnetArgs$weights, foldid = glmnetArgs$foldid,
                          alpha = alpha, family = family, lambda.min.ratio = 0.001,
                          keep = TRUE)
cvFit = cvFitList[[1]]

plot(cvFit)

fitResult = cvFit$glmnet.fit
lambda = cvFit$lambda.min

plotCoefficients(fitResult, lambda, classLevels = classLevels) +
  scale_fill_brewer(type = 'qual', palette = 3) +
  theme(axis.text.y = element_text(size = 8))

annoNames = c('study', 'class')
annoLevels = list(discoveryStudyNames, classLevels)
names(annoLevels) = annoNames
annoColors = list(brewer.pal(3, 'Dark2'), brewer.pal(3, 'Paired'))
names(annoColors) = annoNames
names(annoColors[[1]]) = annoLevels[[1]]
names(annoColors[[2]]) = annoLevels[[2]]

plotExpressionHeatmap(fitResult, lambda, ematDiscovery, sampleMetadata, annoLevels,
                      annoColors, classLevels = classLevels, fontsize_row = 7)

predsList = metapredict(ematList, studyMetadata, sampleMetadata, discoveryStudyNames,
                        alpha = alpha, lambda = lambda, weights = glmnetArgs$weights,
                        family = family)
#> Found4batches
#> Adjusting for0covariate(s) or covariate level(s)
#> Standardizing Data across genes
#> Fitting L/S model and finding priors
#> Finding parametric adjustments
#> Adjusting the Data

calcConfusionValidation(predsList, lambda = lambda, sampleMetadata = sampleMetadata,
                        classLevels = classLevels)
#> $GSE30219
#>          predictedClass
#> trueClass AD SQ SCLC
#>      AD   85  0    0
#>      SQ    9 52    0
#>      SCLC  4  0   17

plotClassProbsValidation(predsList, sampleMetadata, className, classLevels, size = 1.5,
                         ggplotArgs = list(scale_color_brewer(type = 'qual', palette = 3)))
