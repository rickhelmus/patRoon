install.packages(c("pkgdown", "bookdown", "DiagrammeR", "rsvg", "tinytex", "webshot"))
remotes::install_github("rich-iannone/DiagrammeRsvg")

tinytex::install_tinytex()
webshot::install_phantomjs()
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", "/home/patRoon/bin"))

pkgdown::clean_site()
pkgdown::build_site(examples = FALSE)

# make book after pkgdown, otherwise it complains
# NOTE: set clean_envir to FALSE so that 'out' variable below is recognized
out <- file.path(normalizePath("docs", mustWork = FALSE), "handbook_bd")
withr::with_dir("vignettes/handbook/", bookdown::render_book("_index.Rmd",
                                                             output_dir = out,
                                                             clean_envir = FALSE))

# PDF versions
withr::with_dir("vignettes/handbook/", bookdown::render_book("_index.Rmd",
                                                             "bookdown::pdf_book",
                                                             output_dir = out,
                                                             clean_envir = FALSE))
rmarkdown::render("vignettes/tutorial.Rmd", "rmarkdown::pdf_document",
                  output_dir = file.path("docs", "articles"))
refp <- file.path("docs/reference")
devtools::build_manual(path = refp)
file.rename(file.path(refp, paste0(desc::desc_get_field("Package"), "_", desc::desc_get_version(), ".pdf")),
            file.path(refp, "patRoon.pdf"))