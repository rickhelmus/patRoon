compScorings <- read.csv(system.file("data-raw", "compounds-scorings.csv", package = "patRoon"), stringsAsFactors = FALSE)

# adapted from GenForm/msaddion.cpp
adductsGF <- data.table::fread(system.file("data-raw", "adducts-genform.csv", package = "patRoon"))

# from http://ipb-halle.github.io/MetFrag/projects/metfragcl/ and Constants.java in MetFragLib
adductsMF <- data.table::fread(system.file("data-raw", "adducts-metfrag.csv", package = "patRoon"))

usethis::use_data(compScorings, adductsGF, adductsMF, internal = TRUE, overwrite = TRUE)
