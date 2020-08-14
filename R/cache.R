getCacheMode <- function()
{
    ret <- getOption("patRoon.cache.mode", default = "both")
    if (!ret %in% c("both", "save", "load", "none"))
        stop(sprintf("Invalid cache mode (patRoon.cache.mode = \"%s\"), should be \"both\", \"save\", \"load\" or \"none\".", ret))
    return(ret)
}

getCacheFile <- function() getOption("patRoon.cache.fileName")
getMaxCacheEntries <- function() 100000 # UNDONE

makeHash <- function(...)
{
    args <- list(...)

    # strip DT self refs as they sometimes mess up hashing
    args <- recursiveApplyDT(args, function(dt) prepareDTForComparison(copy(dt)), sapply, simplify = FALSE)

    return(digest::digest(args, algo = "xxhash64"))
}

makeFileHash <- function(...) digest::digest(sapply(list(...), digest::digest, file = TRUE, algo = "xxhash64"))

# storing/retrieving R objects: http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/

openCacheDB <- function(file = getCacheFile()) DBI::dbConnect(RSQLite::SQLite(), file)
closeCacheDB <- function(db) DBI::dbDisconnect(db)
openCacheDBScope <- withr::local_(function(x, file = getCacheFile()) openCacheDB(file), function(x) closeCacheDB(x))

loadCacheData <- function(category, hashes, dbArg = NULL)
{
    if (getCacheMode() == "save" || getCacheMode() == "none")
        return(NULL)

    if (is.null(dbArg))
        db <- openCacheDBScope()
    else
        db <- dbArg

    ret <- NULL

    if (nrow(DBI::dbGetQuery(db, sprintf("SELECT 1 FROM sqlite_master WHERE type='table' AND name='%s'", category))) > 0)
    {
        if (length(hashes) == 1) # select only one?
        {
            df <- DBI::dbGetQuery(db, sprintf("SELECT data FROM %s WHERE hash='%s'", category, hashes))

            if (nrow(df) > 0)
                ret <- lapply(df$data, function(x) unserialize(fst::decompress_fst(x)))[[1]]
        }
        else
        {
            df <- DBI::dbGetQuery(db, sprintf("SELECT hash,data FROM %s WHERE hash IN (%s)", category,
                                              paste0(sprintf("'%s'", hashes), collapse = ",")))

            if (nrow(df) > 0)
            {
                ret <- lapply(df$data, function(x) unserialize(fst::decompress_fst(x)))
                if (length(ret) > 0)
                    names(ret) <- df$hash
            }
        }
    }

    ret <- recursiveApplyDT(ret, setalloccol, sapply, simplify = FALSE)
    
    return(ret)
}

saveCacheData <- function(category, data, hash, dbArg = NULL)
{
    if (getCacheMode() == "load" || getCacheMode() == "none")
        return(NULL)

    if (is.null(dbArg))
        db <- openCacheDBScope()
    else
        db <- dbArg

    DBI::dbExecute(db, sprintf("CREATE TABLE IF NOT EXISTS %s (hash TEXT UNIQUE, data BLOB)", category))

    df <- data.frame(d = I(list(fst::compress_fst(serialize(data, NULL, xdr = FALSE)))))

    # From https://stackoverflow.com/a/7353236: update if already exists, otherwise insert
    DBI::dbExecute(db, sprintf("INSERT OR IGNORE INTO %s VALUES ('%s', :d)", category, hash), params=df)
    DBI::dbExecute(db, sprintf("UPDATE %s SET data=(:d) WHERE changes()=0 AND hash='%s'", category, hash), params=df)

    # remove first row (from https://www.experts-exchange.com/questions/24926777/Delete-first-row-of-table.html) if
    # too many rows
    if (DBI::dbGetQuery(db, sprintf("SELECT Count(*) FROM %s", category))[[1]] > getMaxCacheEntries())
        DBI::dbExecute(db, sprintf("DELETE FROM %s WHERE ROWID in (SELECT min(ROWID) FROM %s)", category, category))
}

loadCacheSet <- function(category, setHash, dbArg = NULL)
{
    if (getCacheMode() == "save" || getCacheMode() == "none")
        return(NULL)

    dataHashes <- loadCacheData(category, setHash, dbArg)

    if (is.null(dataHashes))
        return(NULL)

    if (length(dataHashes) == 0)
        return(list())

    return(loadCacheData(category, dataHashes, dbArg))
}

saveCacheSet <- function(category, dataHashes, setHash, dbArg = NULL)
{
    if (getCacheMode() == "load" || getCacheMode() == "none")
        return(NULL)

    # check for duplicate hashes & remove them.
    dataHashes <- unique(dataHashes)

    # stopifnot(!any(duplicated(dataHashes)))
    saveCacheData(category, dataHashes, setHash, dbArg)
}

# UNDONE: include in general cache doc?

#' Clearing of cached data
#'
#' Remove (part of) the cache database used to store (intermediate) processing
#' results.
#'
#' This function will either remove one or more tables within the cache
#' \code{sqlite} database or simply wipe the whole cache file. Removing tables
#' will \code{VACUUM} the database, which may take some time for large cache
#' files.
#'
#' @param what This argument describes what should be done. When \code{what =
#'   NULL} this function will list which tables are present along with an
#'   indication of their size (database rows). If \code{what = "all"} then the
#'   complete file will be removed. Otherwise, \code{what} should be a character
#'   string (a regular expression) that is used to match the table names that
#'   should be removed.
#' @param file The cache file. If \code{NULL} then the value of the
#'   \code{patRoon.cache.fileName} option is used.
#'
#' @export
clearCache <- function(what = NULL, file = NULL)
{
    checkmate::assertString(what, na.ok = FALSE, null.ok = TRUE)
    
    if (!is.null(file))
        checkmate::assertFile(file, "r")
    else
        file <- getCacheFile()
    
    if (!file.exists(file))
        printf("No cache file found, nothing to do ...\n")
    else if (!is.null(what) && what == "all")
    {
        cat("Clearing ALL caches...\n")
        if (unlink(file) != 0)
        {
            gc() # might be some orphaned connection open
            if (unlink(file) != 0)
                warning("Could not clear cache file!")
        }
    }
    else
    {
        db <- openCacheDBScope(file = file)
        tables <- DBI::dbListTables(db)

        if (length(tables) == 0)
            printf("Cache file is empty, nothing to do ...\n")
        else if (is.null(what) || !nzchar(what))
        {
            tableRows <- unlist(sapply(tables, function(tab) DBI::dbGetQuery(db, sprintf("SELECT Count(*) FROM %s", tab))))
            printf("Please specify which cache you want to remove. Available are:\n%s",
                   paste0(sprintf("- %s (%d rows)\n", tables, tableRows), collapse = ""))
            printf("- all (removes complete cache database)\n")
        }
        else
        {
            matchedTables <- grep(what, tables, value = TRUE)
            if (length(matchedTables) == 0)
                printf("No cache found that matches given pattern. Currently stored caches: %s\n", paste0(tables, collapse = ", "))
            else
            {
                for (tab in matchedTables)
                    DBI::dbExecute(db, sprintf("DROP TABLE IF EXISTS %s", tab))
                DBI::dbExecute(db, "VACUUM")
                printf("Removed caches: %s\n", paste0(matchedTables, collapse = ", "))
            }
        }

    }
    invisible(NULL)
}
