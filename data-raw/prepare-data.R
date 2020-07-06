compScorings <- read.csv(system.file("data-raw", "compounds-scorings.csv", package = "patRoon"), stringsAsFactors = FALSE)
defIDLevelRules <- read.csv(system.file("data-raw", "id-level-rules.csv", package = "patRoon"), stringsAsFactors = FALSE)
usethis::use_data(compScorings, defIDLevelRules, internal = TRUE, overwrite = TRUE)
