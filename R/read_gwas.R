#' read_gwas.R
#'
#' Read 1 GWAS Catalog summary-statistics file into the current environment.
#'
#' Strategy: download each .tsv.gz to a local cache directory, then parse with
#' data.table::fread (which transparently handles gzip via local path). This is
#' more robust than streaming over HTTPS for files of this size (~1 GB each):
#' downloads can be resumed/retried, and re-reads are cheap once cached.
#'
#'
#'
#' @param url        Full https URL to a .tsv.gz file on the GWAS Catalog FTP.
#' @param cache_dir  Directory to keep downloaded files (default: tempdir()).
#'                   Use a persistent path (e.g. "~/gwas_cache") to avoid
#'                   re-downloading across R sessions.
#' @param overwrite  If TRUE, redownload even if a cached file exists.
#' @param max_tries  Number of download attempts before giving up.
#' @param fread_args Named list of extra args forwarded to data.table::fread
#'                   (e.g. list(select = c("variant_id", "p_value"))).
#' @return A data.table with the downloaded data
#'
#' @importFrom data.table fread
#' @importFrom curl curl_download
#'
#' @export
#' @examples
#' \donttest{
#' # Ensure internet is available before running the example
#' if (curl::has_internet()) {
#'   try({
#'     dat <- read_gwas_one(paste0("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
#'     "GCST005001-GCST006000/GCST005527/tsoi_2012_23143594_pso_efo0000676_1_ichip.sumstats.tsv.gz"))
#'     head(dat)
#'   })
#' }
#' }
read_gwas_one <- function(url,
                          cache_dir  = tempdir(),
                          overwrite  = FALSE,
                          max_tries  = 5,
                          fread_args = list()) {
  # check just 1 url
  stopifnot(is.character(url), length(url) == 1L)

  # data.table::fread relies on R.utils to decompress .gz files; it is a soft
  # (Suggests) dependency, so check for it up front rather than failing midway.
  if (grepl("\\.gz$", url, ignore.case = TRUE) &&
      !requireNamespace("R.utils", quietly = TRUE)) {
    stop("Reading gzipped (.gz) files requires the 'R.utils' package; ",
         "please install it with install.packages('R.utils').")
  }

  # create directory for download
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(cache_dir, basename(url))

  # check if exists
  need_download <- overwrite || !file.exists(dest) || file.info(dest)$size == 0
  # need to download
  if (need_download) {
    message(sprintf("Downloading %s -> %s", url, dest))
    # curl::curl_download handles HTTPS, redirects, and large files efficiently;
    # wrap in a retry loop because the EBI FTP mirror occasionally drops conns.
    attempt <- 0L
    repeat {
      attempt <- attempt + 1L
      ok <- tryCatch({
        curl::curl_download(url, destfile = dest, mode = "wb", quiet = FALSE)
        TRUE
      }, error = function(e) {
        message(sprintf("  attempt %d failed: %s", attempt, conditionMessage(e)))
        FALSE
      })
      if (ok) break
      if (attempt >= max_tries)
        stop(sprintf("Failed to download %s after %d attempts.", url, max_tries))
      Sys.sleep(2 ^ attempt)  # exponential backoff
    }
  } else {
    message(sprintf("Using cached file: %s (%.1f MB)",
                    dest, file.info(dest)$size / 1e6))
  }

  message(sprintf("Reading %s with data.table::fread ...", basename(dest)))
  do.call(data.table::fread, c(list(input = dest), fread_args))
}

#' Read 2-5 GWAS Catalog summary-statistics files into the current environment.
#'
#' @param urls  Character vector of 2-5 https URLs to .tsv.gz files.
#'              If the vector has names, those names are used in the returned
#'              list; otherwise entries are named gwas1, gwas2, ...
#' @param ...   Forwarded to read_gwas_one (cache_dir, overwrite, max_tries,
#'              fread_args).
#' @return Named list of data.tables, one per input URL.
#'
#' @export
#' @examples
#' \donttest{
#' # Ensure internet is available before running the example
#' if (curl::has_internet()) {
#'   try({
#'     dat <- read_gwas_set(c(
#'       paste0("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
#'       "GCST005001-GCST006000/GCST005527/tsoi_2012_23143594_pso_efo0000676_1_ichip.sumstats.tsv.gz"),
#'       paste0("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
#'       "GCST003001-GCST004000/GCST003044/CD_trans_ethnic_association_summ_stats_b37.txt.gz")
#'     ))
#'     gwas1 <- dat$gwas1
#'     gwas2 <- dat$gwas2
#'     head(gwas1)
#'   })
#' }
#' }
read_gwas_set <- function(urls, ...) {
  if (!is.character(urls))
    stop("'urls' must be a character vector.")
  n <- length(urls)
  if (n < 2L || n > 5L)
    stop(sprintf("'urls' must contain between 2 and 5 entries (got %d).", n))
  if (any(!nzchar(urls)) || any(is.na(urls)))
    stop("'urls' must not contain NA or empty strings.")

  # dataset names
  out_names <- names(urls)
  if (is.null(out_names) || any(!nzchar(out_names)))
    out_names <- paste0("gwas", seq_len(n))

  # repeatedly call read_gwas_one()
  out <- vector("list", n)
  names(out) <- out_names
  for (i in seq_len(n)) {
    message(sprintf("[%d/%d] %s", i, n, urls[[i]]))
    out[[i]] <- read_gwas_one(urls[[i]], ...)
  }
  out
}
