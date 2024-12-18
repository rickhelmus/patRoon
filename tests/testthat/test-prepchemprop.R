# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

context("prepare chem properties")

# row 1: invalid formula
# row 2: charged molecule
tab <- data.table(SMILES = c("C1=CC=C2C(=C1)C=CC3=CC=CC=C3N2C(=O)N", "CCCC(=O)[O-]"),
                  formula = c("CH6", "C4H7O2"))
tabPrep <- prepareChemTable(tab, prefCalcChemProps = TRUE, neutralChemProps = FALSE)
tabPrepNoPref <- prepareChemTable(tab, prefCalcChemProps = FALSE, neutralChemProps = FALSE)
tabPrepNeut <- prepareChemTable(tab, prefCalcChemProps = TRUE, neutralChemProps = TRUE)

test_that("prepareChemTable works", {
    checkmate::expect_data_table(tabPrep, any.missing = FALSE, nrows = nrow(tab))
    checkmate::expect_names(names(tabPrep), must.include = c("SMILES", "InChI", "InChIKey", "formula", "neutralMass"))
    expect_equal(tabPrep$formula[1], "C15H12N2O")
    expect_equal(tabPrepNoPref$formula[1], tab$formula[1])
    checkmate::expect_names(names(tabPrepNeut), must.include = "molNeutralized")
    expect_equal(tabPrepNeut$SMILES[2], "CCCC(=O)O")
    expect_equal(tabPrepNeut$formula[2], "C4H8O2")
    expect_gt(tabPrepNeut$neutralMass[2], tabPrep$neutralMass[2])
})
