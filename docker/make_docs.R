# HACK HACK HACK: recently bookdown started to ignore files starting with
# underscores. This mechanism is also used for pkgdown to ignore files, which we
# need as it shouldn't include the bookdown sub-Rmd files in the website. For
# now we temporarily re-name the files to start with underscores and rename back
# once finished building the website...
handbookSubRmdPath <- file.path("vignettes", "handbook")
origRmds <- list.files(handbookSubRmdPath, pattern = "\\.Rmd")
disabledRmds <- paste0("_", origRmds)
origRmds <- file.path(handbookSubRmdPath, origRmds)
disabledRmds <- file.path(handbookSubRmdPath, disabledRmds)
file.rename(origRmds, disabledRmds)

install.packages(c("pkgdown", "bookdown", "DiagrammeR", "rsvg", "webshot"))
remotes::install_github("rich-iannone/DiagrammeRsvg")

webshot::install_phantomjs()

pkgdown::clean_site()
pkgdown::build_site(examples = FALSE)

# HACK HACK HACK: re-enable bookdown files again...
file.rename(disabledRmds, origRmds)

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