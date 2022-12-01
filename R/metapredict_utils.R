#' Install custom CDF packages from Brainarray.
#'
#' Install Brainarray custom CDFs for processing raw Affymetrix data. See
#' <http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp>.
#'
#' @param pkgs character vector of package names, e.g., 'hgu133ahsentrezgcdf'
#' @param ver integer version number (25 as of 5 Jan 2021)
#'
#' @export
installCustomCdfPackages = function(pkgs, ver = 25) {
  for (pkg in pkgs) {
    pkgUrl = sprintf(
      'http://mbni.org/customcdf/%d.0.0/entrezg.download/%s_%d.0.0.tar.gz',
      ver, pkg, ver)
    utils::install.packages(pkgUrl, repos = NULL)}}


#' Download custom CDF mapping files from Brainarray.
#'
#' Download Brainarray custom CDF mapping files, which are used for mapping
#' probes to genes in datasets whose `studyDataType` is 'affy_series_matrix'.
#' See
#' \url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp}.
#'
#' @param cdf `data frame` with columns `download` (e.g., 'Mouse4302_Mm_ENTREZ')
#'   and `rename` (e.g., 'mouse4302mmentrezgcdf').
#' @param path directory into which to download the files.
#' @param ver integer version number (25 as of 5 Jan 2021).
#'
#' @export
downloadCustomCdfMappings = function(cdf, path = '.', ver = 25) {
  if (!dir.exists(path)) {
    dir.create(path)}
  for (ii in seq_len(nrow(cdf))) {
    temp = tempfile()
    utils::download.file(
      sprintf('http://mbni.org/customcdf/%d.0.0/entrezg.download/%s_%d.0.0.zip',
              ver, cdf$download[ii], ver), temp)
    utils::unzip(temp, files = paste0(cdf$download[ii], '_mapping.txt'),
                 exdir = path)
    file.rename(file.path(path, paste0(cdf$download[ii], '_mapping.txt')),
                file.path(path, paste0(cdf$rename[ii], '_mapping.txt')))
    unlink(temp)}}


mergeDataTable = function(sample, sampleMetadata) {
  merge(data.table(sample = sample), sampleMetadata, by = 'sample', sort = FALSE)}
