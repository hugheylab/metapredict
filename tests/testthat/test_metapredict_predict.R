library(data.table)
parentFolderPath = 'test_data'
foldidColname = 'study'
discoveryStudyNames = c('Bhattacharjee', 'GSE11969', 'GSE29016')

studyMetadataPath = file.path(parentFolderPath, 'study_metadata.csv')
studyMetadata = read.csv(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path(parentFolderPath, 'sample_metadata.csv')
sampleMetadata = read.csv(sampleMetadataPath, stringsAsFactors = FALSE)

glmnetArgsControl = readRDS(file.path(parentFolderPath, 'glmnetArgs'))

test_that('makeGlmnetArgs', {
  glmnetArgsTest = makeGlmnetArgs(sampleMetadata[sampleMetadata$study %in% discoveryStudyNames,],
                                  foldidColname = foldidColname)
  expect_true(all.equal(glmnetArgsTest, glmnetArgsControl, check.attributes = FALSE))
})
