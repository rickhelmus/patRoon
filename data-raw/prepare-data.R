# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

compScorings <- read.csv(system.file("data-raw", "compounds-scorings.csv", package = "patRoon"), stringsAsFactors = FALSE)

# adapted from GenForm/msaddion.cpp
adductsGF <- data.table::fread(system.file("data-raw", "adducts-genform.csv", package = "patRoon"))

# from http://ipb-halle.github.io/MetFrag/projects/metfragcl/ and Constants.java in MetFragLib
adductsMF <- data.table::fread(system.file("data-raw", "adducts-metfrag.csv", package = "patRoon"))

TPsLogicTransformations <- data.table::fread(system.file("data-raw", "TP-logic.csv", package = "patRoon"))

# obtain PubChem transformations and prepare data
transFile <- tempfile(fileext = ".tsv")
utils::download.file("https://zenodo.org/records/10998412/files/PubChem_all_transformations_wExtraInfo.csv?download=1",
                     destfile = transFile)

PubChemTransformations <- data.table::fread(transFile)

# convert column names to generic patRoon format
setnames(PubChemTransformations,
         c("predecessor", "predecessorcid", "successor", "successorcid", "predecessorEM", "successorEM", "XlogPDiff"),
         c("parent_name", "parent_CID", "TP_name", "TP_CID", "parent_neutralMass", "TP_neutralMass", "LogPDiff"))

# merge duplicates rows: these may occur eg when different sources exist
# NOTE: this doesn't take into account different CIDs, hence, only the first will remain
nameCols <- c("parent_name", "TP_name")
otherCols <- c("transformation", "evidencedoi", "evidenceref", "sourcecomment", "sourcecommentfull",
               "datasetdoi", "datasetref", "enzyme", "biosystem")
PubChemTransformations[, (otherCols) := lapply(.SD, function(x) paste0(unique(x), collapse = "/")),
                       by = nameCols, .SDcols = otherCols]
PubChemTransformations <- unique(PubChemTransformations, by = nameCols)

splitIn2Lists <- function(vals) rbindlist(lapply(vals, function(x) as.list(strsplit(x, ">>")[[1]])))

# get SMILES/IKs
PubChemTransformations[, c("parent_SMILES", "TP_SMILES") := splitIn2Lists(ReactionSMILES)]
PubChemTransformations[, c("parent_InChIKey", "TP_InChIKey") := splitIn2Lists(IK_to_IK)]

# add other chem properties
PubChemTransformations[, c("parent_formula", "TP_formula", "parent_InChI", "TP_InChI") :=
                              .(patRoon:::babelConvert(parent_SMILES, "smi", "formula", mustWork = TRUE),
                                patRoon:::babelConvert(TP_SMILES, "smi", "formula", mustWork = TRUE),
                                patRoon:::babelConvert(parent_SMILES, "smi", "inchi", mustWork = TRUE),
                                patRoon:::babelConvert(TP_SMILES, "smi", "inchi", mustWork = TRUE))]

PubChemTransformations[is.na(LogPDiff), LogPDiff := 0]

# clear unneeded columns
PubChemTransformations <- PubChemTransformations[, -c("CID_to_CID", "IK_to_IK", "IKFB_to_IKFB",
                                                      "ReactionSMILES", "DescReactionSMILES", "MassDiff")]

# UNDONE: for now remove some TPs w/out names
PubChemTransformations <- PubChemTransformations[!is.na(TP_name)]

usethis::use_data(compScorings, adductsGF, adductsMF, TPsLogicTransformations, PubChemTransformations, internal = TRUE,
                  overwrite = TRUE)
