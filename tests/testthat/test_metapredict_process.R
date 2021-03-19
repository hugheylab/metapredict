test_that('getSupportedPlatforms', {
  platforms = c('GPL180', 'GPL341', 'GPL571', 'GPL885', 'GPL887', 'GPL890', 'GPL962',
    'GPL1053', 'GPL1073', 'GPL1261', 'GPL1291', 'GPL1293', 'GPL1390',
    'GPL1708', 'GPL3921', 'GPL4133', 'GPL4372', 'GPL5645', 'GPL6104',
    'GPL6254', 'GPL6333', 'GPL6480', 'GPL6865', 'GPL6880', 'GPL6884',
    'GPL6885', 'GPL6887', 'GPL6947', 'GPL7015', 'GPL7202', 'GPL8177',
    'GPL10332', 'GPL10379', 'GPL10558', 'GPL10687', 'GPL13607',
    'GPL13730', 'GPL15331', 'GPL15450', 'GPL18721', 'GPL20769')
  platformsFunc = getSupportedPlatforms()
  expect_equal(platforms, platformsFunc)
})
