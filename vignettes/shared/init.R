# workaround for plotting DiagrammeR plots in PDFs, see:
# https://github.com/rich-iannone/DiagrammeR/issues/133
# NOTE: the final size is determined by the out.width chunk option for PDFs
plotGV <- function(gvCode, ...)
{
    gv <- DiagrammeR::grViz(gvCode, ...)
    if (knitr::is_html_output())
        return(gv) # nothing else to do for HTML
    outf <- tempfile(fileext = ".png")
    utils::capture.output(rsvg::rsvg_png(charToRaw(DiagrammeRsvg::export_svg(gv)), file = outf))
    return(knitr::include_graphics(outf))
}

# only change out.width for PDF
PDFOutWidth <- function(w) if (knitr::is_html_output()) "100%" else w

makeTable <- function(tab)
{
    if (knitr::is_html_output())
    {
        k <- knitr::kable(tab, format = "html")
        k <- kableExtra::kable_styling(k, font_size = 11)
        k <- kableExtra::scroll_box(k, extra_css = "overflow-y: auto; height: 350px;")
    }
    else
    {
        k <- knitr::kable(tab, format = "latex", longtable = TRUE, booktabs = TRUE)
        # k <- kableExtra::kable_styling(k, latex_options = "repeat_header", font_size = 7)
        # UNDONE: repeat_header doesn't work anymore --> make bug report
        k <- kableExtra::kable_styling(k, font_size = 7)
    }
    return(k)
}

# otherwise Linux will get into memory troubles...
knitr::knit_meta("latex_dependency", clean = TRUE)

knitr::opts_chunk$set(
    collapse = FALSE,
    comment = "#>",
    warning = FALSE,
    message = FALSE
)

options(patRoon.cache.fileName = "~/tutorial.sqlite")
options(patRoon.progress.opts = list(style = 1, width = 80))
options(width=1500) # big nr to avoid wrapping text in code chunks
options(datatable.print.nrows = 15)

# pkgdown doesn't seem to load RProfile, which may become annoying when MetFrag
# paths and other patRoon settings are set there...
if (file.exists("~/.Rprofile"))
    source("~/.Rprofile")

library(patRoon)
