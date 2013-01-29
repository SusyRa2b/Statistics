void makeFitTrees(){


  gROOT->LoadMacro("reducedTree.C+");


  float lumiForMCweight = 19.399; //in fb-1
  // this needs to be hardcoded below in the definition of makeTree
  


  //TTbar (3x MG, 1x POWHEG, 1x MC@NLO)
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596ra2b_v66.root","files20fb/TT_FullLept.root", false, false, -1.); 
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606ra2b_v66.root","files20fb/TT_SemiLept.root", false, false, -1.);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613ra2b_v66.root","files20fb/TT_FullHad.root", false, false, -1.);

  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1558ra2b_v66.root", "files20fb/TT-powheg.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TT_8TeV-mcatnlo_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1599ra2b_v66.root", "files20fb/TT-MCatNLO.root", false, false, -1);

  //Single Top (6)

  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1533ra2b_v66.root","files20fb/Tbar-s.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1562ra2b_v66.root","files20fb/Tbar-t.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1534ra2b_v66.root","files20fb/Tbar-tW.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1530ra2b_v66.root","files20fb/T-s.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1531ra2b_v66.root","files20fb/T-t.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1532ra2b_v66.root","files20fb/T-tW.root", false, false, -1);
  
  //WJets (3)
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJetsToLNu_HT-250To300_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1611ra2b_v66.root","files20fb/WJets-250to300.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJetsToLNu_HT-300To400_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1610ra2b_v66.root","files20fb/WJets-300to400.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1587ra2b_v66.root","files20fb/WJets-400.root", false, false, -1);

  //DY (2) -- use extensions
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1526ra2b_v66.root","files20fb/DY-200to400_orig.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1608ra2b_v66.root","files20fb/DY-200to400_ext.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1535ra2b_v66.root","files20fb/DY-400_orig.root", false, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1595ra2b_v66.root","files20fb/DY-400_ext.root", false, false, -1);

  //QCD (9)
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1513ra2b_v66.root","files20fb/QCD-120to170.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1603ra2b_v66.root","files20fb/QCD-170to300.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1561ra2b_v66.root","files20fb/QCD-170to300-smallSample.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1609ra2b_v66.root","files20fb/QCD-300to470.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1560ra2b_v66.root","files20fb/QCD-300to470-smallSample.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1515ra2b_v66.root","files20fb/QCD-470to600.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1516ra2b_v66.root","files20fb/QCD-600to800.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1559ra2b_v66.root","files20fb/QCD-800to1000.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577ra2b_v66.root","files20fb/QCD-1000to1400.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1578ra2b_v66.root","files20fb/QCD-1400to1800.root", true, false, -1);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1585ra2b_v66.root","files20fb/QCD-1800.root", true, false, -1);

  //Zinv (3) - use extensions
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525ra2b_v66.root", "files20fb/Zinv-100to200_orig.root", false, false, -1); 
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607ra2b_v66.root", "files20fb/Zinv-100to200_ext.root", false, false, -1); 
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524ra2b_v66.root", "files20fb/Zinv-200to400_orig.root", false, false, -1); 
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594ra2b_v66.root", "files20fb/Zinv-200to400_ext.root", false, false, -1); 
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523ra2b_v66.root", "files20fb/Zinv-400_orig.root", false, false, -1); 
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602ra2b_v66.root", "files20fb/Zinv-400_ext.root", false, false, -1); 

  //Diboson (3)
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1552ra2b_v66.root", "files20fb/WZ.root", false, false, -1); 
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1563ra2b_v66.root", "files20fb/WW.root", false, false, -1); 
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1551ra2b_v66.root", "files20fb/ZZ.root", false, false, -1); 
 
 // ttV (2)
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1605ra2b_v66.root", "files20fb/ttW.root", false, false, -1. );
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1604ra2b_v66.root", "files20fb/ttZ.root", false, false, -1. );
 
  //data 
  // for data, the second option allows for possible blinding. The weight is set by the last parameter and should always be 1.
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012D-PromptReco-v1_AOD_UCSB1623ra2b_v67.root", "files20fb/HTMHT_2012D1.root", false, false, 1.);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012D-PromptReco-v1_AOD_UCSB1624ra2b_v67.root", "files20fb/MET_2012D1.root", false, false, 1.);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012D-PromptReco-v1_AOD_UCSB1625ra2b_v67.root", "files20fb/JetHT_2012D1.root", false, false, 1.);

  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012D-PromptReco-v1_AOD_UCSB1636ra2b_v67.root", "files20fb/HTMHT_2012D2.root", false, false, 1.);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012D-PromptReco-v1_AOD_UCSB1638ra2b_v67.root", "files20fb/MET_2012D2.root", false, false, 1.);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012D-PromptReco-v1_AOD_UCSB1639ra2b_v67.root", "files20fb/JetHT_2012D2.root", false, false, 1.);

  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012A-13Jul2012-v1_AOD_UCSB1492ra2b_v66.root", "files20fb/MET_2012A.root", false, false, 1.0);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012B-13Jul2012-v1_AOD_UCSB1506ra2b_v66.root", "files20fb/MET_2012B.root", false, false, 1.0);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012C-24Aug2012-v1_AOD_UCSB1505ra2b_v66.root", "files20fb/MET_2012C_rr.root", false, false, 1.0);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012C-PromptReco-v2_AOD_UCSB1508ra2b_v66.root", "files20fb/MET_2012C_pr.root", false, false, 1.0);

  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012B-13Jul2012-v1_AOD_UCSB1501ra2b_v66.root", "files20fb/HTMHT_2012B.root", false, false, 1.0);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012C-24Aug2012-v1_AOD_UCSB1504ra2b_v66.root", "files20fb/HTMHT_2012C_rr.root", false, false, 1.0);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012C-PromptReco-v2_AOD_UCSB1502ra2b_v66.root", "files20fb/HTMHT_2012C_pr.root", false, false, 1.0);

  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HT_Run2012A-13Jul2012-v1_AOD_UCSB1491ra2b_v66.root", "files20fb/HT_2012A.root", false, false, 1.0);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012B-13Jul2012-v1_AOD_UCSB1518ra2b_v66.root", "files20fb/JetHT_2012B.root", false, false, 1.0);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012C-24Aug2012-v2_AOD_UCSB1510ra2b_v66.root", "files20fb/JetHT_2012C_rr.root", false, false, 1.0);
  makeTree("v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012C-PromptReco-v2_AOD_UCSB1511ra2b_v66.root", "files20fb/JetHT_2012C_pr.root", false, false, 1.0);


  // will need signal MC, too, I guess.
// do 1-50 for T1tttt
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.1.root", "files20fb/T1tttt_1.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.2.root", "files20fb/T1tttt_2.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.3.root", "files20fb/T1tttt_3.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.4.root", "files20fb/T1tttt_4.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.5.root", "files20fb/T1tttt_5.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.6.root", "files20fb/T1tttt_6.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.7.root", "files20fb/T1tttt_7.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.8.root", "files20fb/T1tttt_8.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.9.root", "files20fb/T1tttt_9.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.10.root", "files20fb/T1tttt_10.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.11.root", "files20fb/T1tttt_11.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.12.root", "files20fb/T1tttt_12.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.13.root", "files20fb/T1tttt_13.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.14.root", "files20fb/T1tttt_14.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.15.root", "files20fb/T1tttt_15.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.16.root", "files20fb/T1tttt_16.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.17.root", "files20fb/T1tttt_17.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.18.root", "files20fb/T1tttt_18.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.19.root", "files20fb/T1tttt_19.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.20.root", "files20fb/T1tttt_20.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.21.root", "files20fb/T1tttt_21.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.22.root", "files20fb/T1tttt_22.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.23.root", "files20fb/T1tttt_23.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.24.root", "files20fb/T1tttt_24.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.25.root", "files20fb/T1tttt_25.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.26.root", "files20fb/T1tttt_26.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.27.root", "files20fb/T1tttt_27.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.28.root", "files20fb/T1tttt_28.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.29.root", "files20fb/T1tttt_29.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.30.root", "files20fb/T1tttt_30.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.31.root", "files20fb/T1tttt_31.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.32.root", "files20fb/T1tttt_32.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.33.root", "files20fb/T1tttt_33.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.34.root", "files20fb/T1tttt_34.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.35.root", "files20fb/T1tttt_35.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.36.root", "files20fb/T1tttt_36.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.37.root", "files20fb/T1tttt_37.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.38.root", "files20fb/T1tttt_38.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.39.root", "files20fb/T1tttt_39.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.40.root", "files20fb/T1tttt_40.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.41.root", "files20fb/T1tttt_41.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.42.root", "files20fb/T1tttt_42.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.43.root", "files20fb/T1tttt_43.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.44.root", "files20fb/T1tttt_44.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.45.root", "files20fb/T1tttt_45.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.46.root", "files20fb/T1tttt_46.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.47.root", "files20fb/T1tttt_47.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.48.root", "files20fb/T1tttt_48.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.49.root", "files20fb/T1tttt_49.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v3_AODSIM_UCSB1601reshuf_v66.50.root", "files20fb/T1tttt_50.root", false, false, 1.0);


  // T2tt and T1bbbb done by Pawandeep as single files
  makeTree("v66_10/RT/reducedTree_T2tt_final.BtagEff05.UpdateEle.root", "files20fb/T2tt.root", false, false, 1.0);
  makeTree("v66_10/RT/reducedTree_T1bbbb_final.BtagEff05.UpateEle.root", "files20fb/T1bbbb.root", false, false, 1.0);

  /////////////////////////////////////
  /////////////////////////////////////
  /////////////////////////////////////


}


// If weight argument is <0, weight3 from reducedTree is used.  Lumi is based on lumiForMCweight
void makeTree(TString inputfile, TString outputfile, bool doQCD, bool blinding, float weight){
  TFile *f = new TFile(inputfile);
  TTree *t;
  f->GetObject("reducedTree",t);
  reducedTree rt(t);
  char* output4loop = outputfile;
  cout << " new output file name = " << output4loop << endl;
  
  if(weight<0) {
    float lumiForMCweight = 19.399;
    double lumi = lumiForMCweight*1000; //number of inverse pb.
    cout << " weight3 from first entry = " << rt.getWeight3() << " and lumi = " << lumi << "/pb" << endl;
    weight = lumi * rt.getWeight3();
    cout << " weight for sample = " << weight << endl;
  }
  
  rt.Loop(output4loop, doQCD, blinding);
 
  TFile f2(outputfile,"update");  
  f2.cd();
  TTree *ptree;
  f2.GetObject("tree",ptree);
  ptree->SetWeight(weight);
  ptree->AutoSave();
  f2.cd();
  f2.Close();
}
