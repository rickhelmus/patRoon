# Old setting conversion

    Code
      readProjectSettings(file.path(getTestDataPath(), "newProject-23.yml"),
      defaultTestDir)
    Output
      $general
      $general$destination
      [1] "test_temp/test-np/default"
      
      $general$scriptFile
      [1] "process-mod.R"
      
      $general$createRStudioProj
      [1] TRUE
      
      $general$IMS
      $general$IMS$mode
      [1] "none"
      
      $general$IMS$limits
      [1] "bruker"
      
      $general$IMS$CCSMethod
      [1] "none"
      
      
      $general$ionization
      [1] "both"
      
      
      $analyses
      $analyses$generateAnaInfo
      [1] "none"
      
      $analyses$analysisTableFileType
      [1] "CSV"
      
      $analyses$analysisTableFileCSV
      [1] "analyses.csv"
      
      $analyses$analysisTableFileCSVPos
      [1] "analyses-pos-mod.csv"
      
      $analyses$analysisTableFileCSVNeg
      [1] "analyses-neg-mod.csv"
      
      $analyses$analysisTableFileR
      [1] "analyses.R"
      
      $analyses$analysisTableFileRPos
      [1] "analyses-pos.R"
      
      $analyses$analysisTableFileRNeg
      [1] "analyses-neg.R"
      
      $analyses$genAnaInfoDynRaw
      [1] ""
      
      $analyses$genAnaInfoDynRawPos
      [1] ""
      
      $analyses$genAnaInfoDynRawNeg
      [1] ""
      
      $analyses$genAnaInfoDynCentroid
      [1] ""
      
      $analyses$genAnaInfoDynCentroidPos
      [1] ""
      
      $analyses$genAnaInfoDynCentroidNeg
      [1] ""
      
      $analyses$genAnaInfoDynIMS
      [1] ""
      
      $analyses$genAnaInfoDynIMSPos
      [1] ""
      
      $analyses$genAnaInfoDynIMSNeg
      [1] ""
      
      $analyses$genAnaInfoDynProfile
      [1] ""
      
      $analyses$genAnaInfoDynProfilePos
      [1] ""
      
      $analyses$genAnaInfoDynProfileNeg
      [1] ""
      
      
      $preTreatment
      $preTreatment$steps
         algorithm        from            to
            <char>      <char>        <char>
      1:      pwiz raw+agilent centroid+mzML
      
      $preTreatment$output
      $preTreatment$output$centroid
      [1] "converted/centroid"
      
      $preTreatment$output$profile
      [1] "converted/profile"
      
      $preTreatment$output$ims
      [1] "converted/ims"
      
      
      $preTreatment$brukerCalib
      $preTreatment$brukerCalib$enabled
      [1] TRUE
      
      $preTreatment$brukerCalib$method
      [1] ""
      
      $preTreatment$brukerCalib$methodPos
      [1] "da-pos"
      
      $preTreatment$brukerCalib$methodNeg
      [1] "da-neg"
      
      
      
      $features
      $features$featAlgo
      [1] "XCMS"
      
      $features$piekParams
      $features$piekParams$filter
      [1] "none"
      
      $features$piekParams$filterIMS
      [1] "none"
      
      $features$piekParams$suspects
      $features$piekParams$suspects$single
      [1] ""
      
      $features$piekParams$suspects$sets
      $features$piekParams$suspects$sets$pos
      [1] ""
      
      $features$piekParams$suspects$sets$neg
      [1] ""
      
      
      
      $features$piekParams$peaksAlgo
      [1] "piek"
      
      
      $features$fGroupsAlgo
      [1] "KPIC2"
      
      $features$IMSPeaksMob
      [1] "piek"
      
      $features$IMSPeaksChrom
      [1] "piek"
      
      $features$suspects
      $features$suspects$single
      [1] ""
      
      $features$suspects$sets
      $features$suspects$sets$pos
      [1] "susps-pos"
      
      $features$suspects$sets$neg
      [1] "susps-neg"
      
      
      
      $features$exSuspList
      [1] FALSE
      
      $features$IMSSuspCCSPred
      [1] "none"
      
      $features$fGroupsAdv
      $features$fGroupsAdv$preIntThr
      [1] 100
      
      $features$fGroupsAdv$intThr
      [1] 10000
      
      $features$fGroupsAdv$repAbundance
      [1] 1
      
      $features$fGroupsAdv$maxRepRSD
      [1] 0.75
      
      $features$fGroupsAdv$blankThr
      [1] 5
      
      $features$fGroupsAdv$removeBlanks
      [1] TRUE
      
      $features$fGroupsAdv$featNorm
      [1] "istd"
      
      $features$fGroupsAdv$groupNorm
      [1] TRUE
      
      $features$fGroupsAdv$retention
      [1]   10 1000
      
      $features$fGroupsAdv$mz
      [1]     100 3000000
      
      $features$fGroupsAdv$ISTDLists
      $features$fGroupsAdv$ISTDLists$single
      [1] ""
      
      $features$fGroupsAdv$ISTDLists$sets
      $features$fGroupsAdv$ISTDLists$sets$pos
      [1] "istd-pos"
      
      $features$fGroupsAdv$ISTDLists$sets$neg
      [1] "istd-neg"
      
      
      
      
      $features$retention
      [1]   10 1000
      
      $features$mz
      [1]     100 3000000
      
      
      $annotations
      $annotations$componAlgo
      [1] "RAMClustR"
      
      $annotations$selectIons
      [1] FALSE
      
      $annotations$formulasAlgo
      [1] "GenForm"
      
      $annotations$compoundsAlgo
      [1] "MetFrag"
      
      $annotations$MSLibraryFormat
      [1] "msp"
      
      $annotations$MSLibraryPath
      [1] ""
      
      $annotations$estIDConf
      [1] "formulas"  "compounds"
      
      $annotations$compCCSPred
      [1] "none"
      
      
      $TPs
      $TPs$TPsAlgo
      [1] "biotransformer"
      
      $TPs$TPGenInput
      [1] "screening"
      
      $TPs$TPSuspectList
      [1] ""
      
      $TPs$TPDoMFDB
      [1] FALSE
      
      
      $report
      $report$reportGen
      [1] "HTML"   "legacy"
      
      $report$reportLegacy
      [1] "CSV" "PDF"
      
      

