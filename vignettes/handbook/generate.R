withr::with_dir("vignettes/handbook/", bookdown::render_book("_index.Rmd", output_dir = "docs/handbook_bd"))
withr::with_dir("vignettes", rmarkdown::render("handbook.Rmd", output_format = "all"))
withr::with_dir("vignettes", rmarkdown::render("tutorial.Rmd", output_format = "all"))
