install.packages(c("pkgdown", "bookdown", "DiagrammeR"))

pkgdown::clean_site()
pkgdown::build_site(examples = FALSE)

# make book after pkgdown, otherwise it complains
# NOTE: set clean_envir to FALSE so that 'out' variable below is recognized
out <- file.path(normalizePath("docs", mustWork = FALSE), "handbook_bd")
withr::with_dir("vignettes/handbook/", bookdown::render_book("_index.Rmd",
                                                             output_dir = out,
                                                             clean_envir = FALSE))
