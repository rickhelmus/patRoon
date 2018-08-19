# From R extensions manual

# libs

files <- Sys.glob(paste0("*", SHLIB_EXT))
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
if(file.exists("symbols.rds"))
    file.copy("symbols.rds", dest, overwrite = TRUE)

# bins

execs <- "GenForm"
if(WINDOWS)
    execs <- paste0(execs, ".exe")

if (any(file.exists(execs)))
{
    dest <- file.path(R_PACKAGE_DIR, paste0('bin', R_ARCH))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    file.copy(execs, dest, overwrite = TRUE)
}
