docs <- normalizePath("docs")
pkgdown::clean_site()
pkgdown::build_site(examples = FALSE)

# make book after pkgdown, otherwise it complains
withr::with_dir("vignettes/handbook/", bookdown::render_book("_index.Rmd", output_dir = file.path(docs, "handbook_bd")))
