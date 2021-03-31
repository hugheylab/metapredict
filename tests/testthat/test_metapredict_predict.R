library(data.table)
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

test_that('makeGlmnetArgs', {
  glmnetArgsTest = makeGlmnetArgs(sampleMetadata[sampleMetadata$study %in% discoveryStudyNames,],
                                  foldidColname = foldidColname)
  expect_true(all.equal(glmnetArgsTest, glmnetArgsControl, check.attributes = FALSE))
})

test_that('metapredictCv', {
  cvFitListTest = metapredictCv(ematDiscoveryControl, sampleMetadata, yName = className,
                                weights = glmnetArgsControl$weights, foldid = glmnetArgsControl$foldid,
                                alpha = alpha, family = family, lambda.min.ratio = 0.001,
                                keep = TRUE)
  expect_true(all.equal(cvFitListTest, cvFitListControl, check.attributes = FALSE))
})
