#' @include main.R
NULL

getTIMSMetaDataFile <- function(path) file.path(path, "analysis.tdf")
openTIMSMetaDBScope <- withr::local_(function(x, f) DBI::dbConnect(RSQLite::SQLite(), getTIMSMetaDataFile(f)),
                                     function(x) DBI::dbDisconnect(x))

getTIMSMetaTable <- function(db, name, cols)
{
    as.data.table(DBI::dbGetQuery(db, sprintf("SELECT %s FROM %s", paste0(cols, collapse = ","), name)))
}
