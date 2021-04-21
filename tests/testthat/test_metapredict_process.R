

test_that('getSupportedPlatforms', {
  supportedPlatformsControl = fread(file.path(parentFolderPath, 'supportedPlatforms.csv'))
  supportedPlatformsTest = getSupportedPlatforms()
  expect_equal(supportedPlatformsTest$platform, supportedPlatformsControl$platform, check.attributes = FALSE)
})

test_that('getUnsupportedPlatforms', {
  unsupPlatformsFunc = getUnsupportedPlatforms(studyMetadata)
  expect_true(length(unsupPlatformsFunc) == 0)
  warnFrame = data.frame(study = as.character(NA), studyDataType = as.character(NA), platformInfo = as.character(NA))
  warnFrame$study = 'GSE11969'
  warnFrame$studyDataType = 'series_matrix'
  warnFrame$platformInfo = 'abc123'
  studyMetadataWarn = rbind(studyMetadata, warnFrame)
  unsupPlatformsFuncWarn = getUnsupportedPlatforms(studyMetadataWarn)
  expect_equal(warnFrame$platformInfo, unsupPlatformsFuncWarn)
})
test_that('getStudyDataList', {
  local_seed(1)
  esetListTest = getStudyDataList(parentFolderPath, studyMetadata)
  expect_equal(esetListTest, esetListControl, check.attributes = FALSE)
})

test_that('extractExpressionData', {
  local_seed(3)
  ematListTest = extractExpressionData(esetListControl, sampleMetadata)
  expect_equal(ematListTest, ematListControl, check.attributes = FALSE)
})
