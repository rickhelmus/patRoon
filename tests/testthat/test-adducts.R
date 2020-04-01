context("adducts")

test_that("adduct class and utilities", {
    expect_known_show(adduct("H"), testFile("adduct-show", text = TRUE))
    
    expect_equal(adduct("H"), as.adduct("[M+H]+"))
    expect_equal(adduct("H"), as.adduct("[M+H]+"))
    expect_equal(adduct(sub = "H", charge = -1), as.adduct("[M-H]-"))
    expect_equal(adduct(add = "Na", sub = "H2", charge = -1), as.adduct("[M-H2+Na]-"))
    expect_equal(adduct(add = "H2", charge = 2), as.adduct("[M+H2]2+"))
    expect_equal(adduct(add = "H", molMult = 3), as.adduct("[3M+H]+"))
    expect_equal(adduct(add = "COOH"), as.adduct("[M+CHO2]+"))

    expect_equal(adduct(add = "H"), as.adduct("[M+H]+", format = "metfrag"))
    expect_equal(adduct(add = "H"), as.adduct(1, format = "metfrag", isPositive = TRUE))
    expect_equal(adduct(charge = -1), as.adduct(0, format = "metfrag", isPositive = FALSE))

    expect_true(all(sapply(GenFormAdducts()$adduct, function(a) is(as.adduct(a, format = "genform"), "adduct"))))
    expect_true(all(sapply(MetFragAdducts()$adduct_type, function(a) is(as.adduct(a, format = "metfrag"), "adduct"))))
    expect_true(all(mapply(MetFragAdducts()$adduct_mode, MetFragAdducts()$charge,
                           FUN = function(a, ch) is(as.adduct(a, format = "metfrag", isPositive = ch>0), "adduct"))))

    expect_equal(as.character(adduct("H")), "[M+H]+")
    expect_equal(as.character(adduct(add = "Na", sub = "H2", charge = -1)), "[M+Na-H2]-")
    expect_equal(as.character(adduct(add = "Na2", charge = 2)), "[M+Na2]2+")

    expect_equal(calculateNeutralFormula(calculateIonFormula("C2H4O", "[M+H]+"), "[M+H]+"), "C2H4O")
})
