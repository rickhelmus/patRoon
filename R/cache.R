getCacheMode <- function() getOption("patRoon.cache.mode")
getCacheFile <- function() getOption("patRoon.cache.fileName")
getMaxCacheEntries <- function() 100000 # UNDONE

makeHash <- function(...)
{
    args <- list(...)

    # strip DT self refs as they sometimes mess up hashing
    args <- recursiveApplyDT(args, function(dt) stripDTRef(copy(dt)), sapply, simplify = FALSE)

    return(digest::digest(args, algo = "xxhash64"))
}

makeFileHash <- function(...) digest::digest(sapply(list(...), digest::digest, file = TRUE, algo = "xxhash64"))

# storing/retrieving R objects: http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/

openCacheDB <- function() dbConnect(SQLite(), getCacheFile())
closeCacheDB <- function(db) dbDisconnect(db)

loadCacheData <- function(category, hashes, dbArg = NULL)
{
    if (getCacheMode() == "save" || getCacheMode() == "none")
        return(NULL)

    if (is.null(dbArg))
        db <- openCacheDB()
    else
        db <- dbArg

    ret <- NULL

    if (nrow(dbGetQuery(db, sprintf("SELECT 1 FROM sqlite_master WHERE type='table' AND name='%s'", category))) > 0)
    {
        if (length(hashes) == 1) # select only one?
        {
            df <- dbGetQuery(db, sprintf("SELECT data FROM %s WHERE hash='%s'", category, hashes))

            if (nrow(df) > 0)
                ret <- lapply(df$data, function(x) unserialize(decompress_fst(x)))[[1]]
        }
        else
        {
            df <- dbGetQuery(db, sprintf("SELECT hash,data FROM %s WHERE hash IN (%s)", category,
                                         paste0(sprintf("'%s'", hashes), collapse = ",")))

            if (nrow(df) > 0)
            {
                ret <- lapply(df$data, function(x) unserialize(decompress_fst(x)))
                if (length(ret) > 0)
                    names(ret) <- df$hash
            }
        }
    }

    if (is.null(dbArg))
        closeCacheDB(db)

    return(ret)
}

saveCacheData <- function(category, data, hash, dbArg = NULL)
{
    if (getCacheMode() == "load" || getCacheMode() == "none")
        return(NULL)

    if (is.null(dbArg))
        db <- openCacheDB()
    else
        db <- dbArg

    dbExecute(db, sprintf("CREATE TABLE IF NOT EXISTS %s (hash TEXT UNIQUE, data BLOB)", category))

    df <- data.frame(d = I(list(compress_fst(serialize(data, NULL, xdr = FALSE)))))

    # From https://stackoverflow.com/a/7353236: update if already exists, otherwise insert
    dbExecute(db, sprintf("INSERT OR IGNORE INTO %s VALUES ('%s', :d)", category, hash), params=df)
    dbExecute(db, sprintf("UPDATE %s SET data=(:d) WHERE changes()=0 AND hash='%s'", category, hash), params=df)

    # remove first row (from https://www.experts-exchange.com/questions/24926777/Delete-first-row-of-table.html) if
    # too many rows
    if (dbGetQuery(db, sprintf("SELECT Count(*) FROM %s", category))[[1]] > getMaxCacheEntries())
        dbExecute(db, sprintf("DELETE FROM %s WHERE ROWID in (SELECT min(ROWID) FROM %s)", category, category))

    if (is.null(dbArg))
        closeCacheDB(db)
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
#' files. (intermediate) processing results.
#'
#' @param what This argument describes what should be done. When \code{what =
#'   NULL} this function will list which tables are present along with an
#'   indication of their size (database rows). If \code{what = "all"} then the
#'   complete file will be removed. Finally, \code{what} can be a character
#'   vector with table names that should be removed.
#'
#' @export
clearCache <- function(what = NULL)
{
    checkmate::assertCharacter(what, min.len = 1, null.ok = TRUE)
    
    if (!file.exists(getCacheFile()))
        printf("No cache file found, nothing to do ...\n")
    else if ("all" %in% what)
    {
        cat("Clearing ALL caches...\n")
        if (unlink(getCacheFile()) != 0)
        {
            gc() # might be some orphaned connection open
            if (unlink(getCacheFile()) != 0)
                warning("Could not clear cache file!")
        }
    }
    else
    {
        db <- openCacheDB()
        tables <- dbListTables(db)

        if (length(tables) == 0)
            printf("Cache file is empty, nothing to do ...\n")
        else if (is.null(what) || !nzchar(what))
        {
            tableRows <- unlist(sapply(tables, function(tab) dbGetQuery(db, sprintf("SELECT Count(*) FROM %s", tab))))
            printf("Please specify which cache you want to remove. Available are:\n%s",
                   paste0(sprintf("- %s (%d rows)\n", tables, tableRows), collapse = ""))
            printf("- all (removes complete cache database)\n")
        }
        else
        {
            nexist <- what[!what %in% tables]
            if (length(nexist) > 0)
                warning(sprintf("Non-existing cache: %s\n", paste0(nexist, collapse = ", ")))

            for (categ in what[what %in% tables])
                dbExecute(db, sprintf("DROP TABLE IF EXISTS %s", categ))

            if (any(what %in% tables))
                dbExecute(db, "VACUUM")
        }

        closeCacheDB(db)
    }
    invisible(NULL)
}
