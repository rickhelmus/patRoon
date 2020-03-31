install.packages(c("pkgdown", "bookdown", "DiagrammeR"))

pkgdown::clean_site()
pkgdown::build_site(examples = FALSE)

# make book after pkgdown, otherwise it complains
withr::with_dir("vignettes/handbook/", bookdown::render_book("_index.Rmd",
                                                             output_dir = file.path(normalizePath("docs"), "handbook_bd")))
