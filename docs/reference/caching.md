# Utilities for caching of workflow data.

Several utility functions for caching workflow data. The most important
function is `clearCache`; other functions are primarily for internal
use.

## Usage

``` r
makeHash(..., checkDT = TRUE)

makeFileHash(..., length = Inf)

loadCacheData(category, hashes, dbArg = NULL, simplify = TRUE, fixDTs = TRUE)

saveCacheData(category, data, hash, dbArg = NULL)

clearCache(what = NULL, file = NULL, vacuum = TRUE)
```

## Arguments

- ...:

  Arguments/objects to be used for hashing.

- checkDT:

  `logical`, set to `TRUE` with (a `list` with) `data.table`s to ensure
  reproducible hashing. Otherwise can be set to `FALSE` to improve
  performance.

- length:

  Maximum file length to hash. Passed to `digest`.

- category:

  The category of the object to be cached.

- hashes:

  A `character` with one more hashes (*e.g.* obtained with `makeHash`)
  of the objects to be loaded.

- dbArg:

  Alternative connection to database. Default is `NULL` and uses the
  cache options as defined by patRoon.cache.fileName. Mainly used
  internally to improve performance.

- simplify:

  If `TRUE` and `length(hashes)==1` then the returned data is returned
  directly, otherwise the data is in a `list`.

- fixDTs:

  Should be `TRUE` if cached data consists of (nested) `data.table`s.
  Otherwise can be `FALSE` to speed up loading.

- data:

  The object to be cached.

- hash:

  The hash string of the object to be cached (*e.g.* obtained with
  `makeHash`).

- what:

  This argument describes what should be done. When `what = NULL` this
  function will list which tables are present along with an indication
  of their size (database rows). If `what = "all"` then the complete
  file will be removed. Otherwise, `what` should be a character string
  (a regular expression) that is used to match the table names that
  should be removed.

- file:

  The cache file. If `NULL` then the value of the
  `patRoon.cache.fileName` option is used.

- vacuum:

  If `TRUE` then the `VACUUM` operation will be run on the cache
  database to reduce the file size. Setting this to `FALSE` might be
  handy to avoid long processing times on large cache databases.

## Details

`makeHash` Make a hash string of given arguments.

`makeFileHash` Generates a hash from the contents of one or more files.

`loadCacheData` Loads cached data from a database.

`saveCacheData` caches data in a database.

`clearCache` will either remove one or more tables within the cache
`sqlite` database or simply wipe the whole cache file. Removing tables
will `VACUUM` the database (unless `vacuum=FALSE`), which may take some
time for large cache files.
