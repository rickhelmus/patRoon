compScorings <- read.csv(system.file("data-raw", "compounds-scorings.csv", package = "patRoon"), stringsAsFactors = FALSE)

# adapted from GenForm/msaddion.cpp
adductsGF <- data.table::fread(system.file("data-raw", "adducts-genform.csv", package = "patRoon"))

# from http://ipb-halle.github.io/MetFrag/projects/metfragcl/ and Constants.java in MetFragLib
adductsMF <- data.table::fread(system.file("data-raw", "adducts-metfrag.csv", package = "patRoon"))

TPsLogicTransformations <- data.table::fread(system.file("data-raw", "TP-logic.csv", package = "patRoon"))

# obtain PubChem transformations and prepare data
transFile <- tempfile(fileext = ".tsv")
utils::download.file("https://git-r3lab.uni.lu/eci/pubchem/-/raw/master/annotations/tps/norman_transformations_tsv.txt",
                     destfile = transFile)

PubChemTransformations <- data.table::fread(transFile, skip = 2)

# first row contains DB var types, no need for those
PubChemTransformations <- PubChemTransformations[-1]

# for now, only keep predecessor --> successor rows
PubChemTransformations <- PubChemTransformations[`#cid` == predecessorcid]

# no need for first CID column anymore
PubChemTransformations <- PubChemTransformations[, -"#cid"]

# convert column names to generic patRoon format
setnames(PubChemTransformations,
         c("predecessor", "predecessorcid", "successor", "successorcid"),
         c("parent_name", "parent_CID", "TP_name", "TP_CID"))

# merge duplicates rows: these may occur eg when different sources exist
# NOTE: this doesn't take into account different CIDs, hence, only the first will remain
nameCols <- c("parent_name", "TP_name")
otherCols <- c("transformation", "evidencedoi", "evidenceref", "sourcecomment", "sourcecommentfull",
               "datasetdoi", "datasetref", "enzyme", "biosystem")
PubChemTransformations[, (otherCols) := lapply(.SD, function(x) paste0(unique(x), collapse = "/")),
                       by = nameCols, .SDcols = otherCols]
PubChemTransformations <- unique(PubChemTransformations, by = nameCols)

# use CIDs to get SMILES with webchem
allCIDs <- union(PubChemTransformations$parent_CID, PubChemTransformations$TP_CID)
PCInfos <- as.data.table(webchem::pc_prop(allCIDs, c("CanonicalSMILES", "XLogP")))

# and add SMILES/XLogPs to parents/TPs...
PubChemTransformations[, c("parent_SMILES", "parent_LogP") :=
                           PCInfos[match(parent_CID, CID), .(CanonicalSMILES, XLogP)]]
PubChemTransformations[, c("TP_SMILES", "TP_LogP") := PCInfos[match(TP_CID, CID), .(CanonicalSMILES, XLogP)]]

# add other chem properties
PubChemTransformations[, c("parent_formula", "TP_formula", "parent_InChI", "TP_InChI",
                           "parent_InChIKey", "TP_InChIKey") :=
                           .(patRoon:::convertToFormulaBabel(parent_SMILES, "smi", mustWork = TRUE),
                             patRoon:::convertToFormulaBabel(TP_SMILES, "smi", mustWork = TRUE),
                             patRoon:::babelConvert(parent_SMILES, "smi", "inchi", mustWork = TRUE),
                             patRoon:::babelConvert(TP_SMILES, "smi", "inchi", mustWork = TRUE),
                             patRoon:::babelConvert(parent_SMILES, "smi", "inchikey", mustWork = TRUE),
                             patRoon:::babelConvert(TP_SMILES, "smi", "inchikey", mustWork = TRUE))]
PubChemTransformations[, c("parent_neutralMass", "TP_neutralMass") :=
                           .(sapply(parent_formula, patRoon:::getFormulaMass),
                             sapply(TP_formula, patRoon:::getFormulaMass))]

usethis::use_data(compScorings, adductsGF, adductsMF, TPsLogicTransformations, PubChemTransformations, internal = TRUE,
                  overwrite = TRUE)
