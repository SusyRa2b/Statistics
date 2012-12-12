void makeFitTrees(){


  gROOT->LoadMacro("reducedTree.C+");

  // Set lumi.  Note: if passing -1, you have to set the lumi in makeTree() instead.
  float lumiForMCweight = 17.61; //in fb-1
  

  //Ben running on v66_7 reducedTrees
  
  //TTbar (3,1,1)
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/TT_FullLept.root", false, false, -1); 
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/TT_SemiLept.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1613ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/TT_FullHad.root", false, false,-1);
 makeTree("/cu4/ra2b/reducedTrees/v66_7/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTbarJetsPowheg.root", "/cu4/kreis/tinyTrees/v66_7_tt8/TT-powheg.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1489ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/TT.root", false, false, -1);

  //Single Top (6)
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1533ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/Tbar-s.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1562ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/Tbar-t.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1534ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/Tbar-tW.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1530ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/T-s.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1531ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/T-t.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1532ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/T-tW.root", false, false, -1);
  
  //WJets (3)
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJetsToLNu_HT-250To300_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1611ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/WJets-250to300.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJetsToLNu_HT-300To400_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1610ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/WJets-300to400.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1587ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/WJets-400.root", false, false, -1);

  //DY (2) -- use extensions
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1526ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/DY-200to400_orig.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1608ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/DY-200to400_ext.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1535ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/DY-400_orig.root", false, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1595ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/DY-400_ext.root", false, false, -1);

  //QCD Pythia (9)
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1513ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-120to170.root", true, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1603ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-170to300.root", true, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_v3_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1609ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-300to470.root", true, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1515ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-470to600.root", true, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1516ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-600to800.root", true, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1516ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-800to1000.root", true, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1577ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-1000to1400.root", true, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1578ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-1400to1800.root", true, false, -1);
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1585ra2b_v66.root","/cu4/kreis/tinyTrees/v66_7_tt8/QCD-1800.root", true, false, -1);

  //Zinv (3) - use extensions
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/Zinv-100to200_orig.root", false, false, -1); 
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/Zinv-100to200_ext.root", false, false, -1); 
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/Zinv-200to400_orig.root", false, false, -1); 
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/Zinv-200to400_ext.root", false, false, -1); 
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/Zinv-400_orig.root", false, false, -1); 
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/Zinv-400_ext.root", false, false, -1); 

  //Diboson (3)
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1552ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/WZ.root", false, false, -1); 
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1563ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/WW.root", false, false, -1); 
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1551ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/ZZ.root", false, false, -1); 

  //Rare (2)
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1605ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/ttW.root", false, false, -1 );
  makeTree("/cu4/ra2b/reducedTrees/v66_7/ORIGINALS/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1604ra2b_v66.root", "/cu4/kreis/tinyTrees/v66_7_tt8/ttZ.root", false, false, -1 );

  //T1bbbb (1) -- pass 1 for weight
  makeTree("/cu4/ra2b/reducedTrees/v66_7/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T1bbbb.root", "/cu4/kreis/tinyTrees/v66_7_tt8/T1bbbb.root", false, false, 1.0);
  
  //Not doing data, madgraph QCD, or T1tttt right now.


  /////////////////////////////////////
  /////////////////////////////////////
  /////////////////////////////////////


   // MC without JER corrections	 cfA v66		     																										 
  //makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1489ra2b_v66.root", "files12p0fb_trkVeto/TT.root", false, false, (lumiForMCweight/12.03)*0.40657);

  //makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1558ra2b_v66.root", "files12p0fb_trkVeto/TT-powheg.root", false, false, (lumiForMCweight/12.03)*0.12987 );

  //makeTree("reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596ra2b_v66.root","files12p0fb_trkVeto/TT_FullLept.root", false, false, (lumiForMCweight/12.03)*0.02599); 
  //makeTree("reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606ra2b_v66.root","files12p0fb_trkVeto/TT_SemiLept.root", false, false, (lumiForMCweight/12.03)*0.04907);
  //makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1586ra2b_v66.root","files12p0fb_trkVeto/TT_FullHad.root", false, false, (lumiForMCweight/12.03)*0.11887);


/*

  
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJets_HT250To300.root", "files12p0fb_trkVeto/WJets-250to300.root", false, false, (lumiForMCweight/12.03)*0.42171);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJets_HT300To400.root", "files12p0fb_trkVeto/WJets-300to400.root", false, false,(lumiForMCweight/12.03)*0.32349 );
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WJets_HT400ToInf.root", "files12p0fb_trkVeto/WJets-400.root", false, false, (lumiForMCweight/12.03)*0.07283);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1530ra2b_v66.root", "files12p0fb_trkVeto/T-s.root", false, false, (lumiForMCweight/12.03)*0.17539 );
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1533ra2b_v66.root", "files12p0fb_trkVeto/Tbar-s.root", false, false,(lumiForMCweight/12.03)*0.15126 );
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1531ra2b_v66.root", "files12p0fb_trkVeto/T-t.root", false, false, (lumiForMCweight/12.03)*28.536);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1562ra2b_v66.root", "files12p0fb_trkVeto/Tbar-t.root", false, false, (lumiForMCweight/12.03)*0.19199);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1532ra2b_v66.root", "files12p0fb_trkVeto/T-tW.root", false, false, (lumiForMCweight/12.03)*0.26832 );  
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1534ra2b_v66.root", "files12p0fb_trkVeto/Tbar-tW.root",false, false, (lumiForMCweight/12.03)*0.27061 );

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD120.root", "files12p0fb_trkVeto/QCD-120to170.root", true, false, (lumiForMCweight/12.03)*326.67);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD170.root", "files12p0fb_trkVeto/QCD-170to300.root", true, false, (lumiForMCweight/12.03)*70.632 ); 
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD300.root", "files12p0fb_trkVeto/QCD-300to470.root", true, false, (lumiForMCweight/12.03)*3.5712);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD470.root", "files12p0fb_trkVeto/QCD-470to600.root", true, false, (lumiForMCweight/12.03)*0.34293);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD600.root", "files12p0fb_trkVeto/QCD-600to800.root", true, false, (lumiForMCweight/12.03)*0.081326 ); 
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD800.root", "files12p0fb_trkVeto/QCD-800to1000.root", true, false, (lumiForMCweight/12.03)*0.010681 ); 

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD1000.root", "files12p0fb_trkVeto/QCD-1000to1400.root", true, false, (lumiForMCweight/12.03)*0.0045193 ); 
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD1400.root", "files12p0fb_trkVeto/QCD-1400to1800.root", true, false, (lumiForMCweight/12.03)*0.00020163);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.QCD1800.root", "files12p0fb_trkVeto/QCD-1800.root", true, false, (lumiForMCweight/12.03)*0.00002251);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Zinvisible_HT100To200.root", "files12p0fb_trkVeto/Zinv-100to200.root", false, false, (lumiForMCweight/12.03)*0.51970);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Zinvisible_HT200To400.root", "files12p0fb_trkVeto/Zinv-200to400.root", false, false, (lumiForMCweight/12.03)*0.11747);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.Zinvisible_HT400ToInf.root", "files12p0fb_trkVeto/Zinv-400.root", false, false, (lumiForMCweight/12.03)*0.07498);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJets_HT200To400.root", "files12p0fb_trkVeto/DY-200to400.root", false, false, (lumiForMCweight/12.03)*0.07437 );
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZJets_HT400ToInf.root", "files12p0fb_trkVeto/DY-400.root", false, false, (lumiForMCweight/12.03)*0.02372);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WZ.root", "files12p0fb_trkVeto/WZ.root", false, false, (lumiForMCweight/12.03)*0.038856 ); 
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.WW.root", "files12p0fb_trkVeto/WW.root", false, false, (lumiForMCweight/12.03)*0.066162 ); 
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.ZZ.root", "files12p0fb_trkVeto/ZZ.root", false, false, (lumiForMCweight/12.03)*0.006378 );

  makeTree("reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1605ra2b_v66.root", "files12p0fb_trkVeto/ttW.root", false, false, (lumiForMCweight/12.03)*0.01319 );
  makeTree("reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1604ra2b_v66.root", "files12p0fb_trkVeto/ttZ.root", false, false, (lumiForMCweight/12.03)*0.00985 );

*/

  //makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.T1bbbb.root", "files12p0fb_trkVeto/T1bbbb.root", false, false, 1.0);

/*
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012A-13Jul2012-v1_AOD_UCSB1492ra2b_v66.root", "files12p0fb_trkVeto/MET_2012A.root", false, true, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012B-13Jul2012-v1_AOD_UCSB1506ra2b_v66.root", "files12p0fb_trkVeto/MET_2012B.root", false, true, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012C-24Aug2012-v1_AOD_UCSB1505ra2b_v66.root", "files12p0fb_trkVeto/MET_2012C_rr.root", false, true, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012C-PromptReco-v2_AOD_UCSB1508ra2b_v66.root", "files12p0fb_trkVeto/MET_2012C_pr.root", false, true, 1.0);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012B-13Jul2012-v1_AOD_UCSB1501_v66.root", "files12p0fb_trkVeto/HTMHT_2012B.root", false, true, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012C-24Aug2012-v1_AOD_UCSB1504_v66.root", "files12p0fb_trkVeto/HTMHT_2012C_rr.root", false, true, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012C-PromptReco-v2_AOD_UCSB1502ra2b_v66.root", "files12p0fb_trkVeto/HTMHT_2012C_pr.root", false, true, 1.0);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HT_Run2012A-13Jul2012-v1_AOD_UCSB1491_v66.root", "files12p0fb_trkVeto/HT_2012A.root", false, true, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012B-13Jul2012-v1_AOD_UCSB1518ra2b_v66.root", "files12p0fb_trkVeto/JetHT_2012B.root", false, true, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012C-24Aug2012-v2_AOD_UCSB1510ra2b_v66.root", "files12p0fb_trkVeto/JetHT_2012C_rr.root", false, true, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012C-PromptReco-v2_AOD_UCSB1511ra2b_v66.root", "files12p0fb_trkVeto/JetHT_2012C_pr.root", false, true, 1.0);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012A-13Jul2012-v1_AOD_UCSB1492ra2b_v66.root", "files12p0fb_trkVeto/MET_2012A.root", false, false, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012B-13Jul2012-v1_AOD_UCSB1506ra2b_v66.root", "files12p0fb_trkVeto/MET_2012B.root", false, false, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012C-24Aug2012-v1_AOD_UCSB1505ra2b_v66.root", "files12p0fb_trkVeto/MET_2012C_rr.root", false, false, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.MET_Run2012C-PromptReco-v2_AOD_UCSB1508ra2b_v66.root", "files12p0fb_trkVeto/MET_2012C_pr.root", false, false, 1.0);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012B-13Jul2012-v1_AOD_UCSB1501ra2b_v66.root", "files12p0fb_trkVeto/HTMHT_2012B.root", false, false, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012C-24Aug2012-v1_AOD_UCSB1504ra2b_v66.root", "files12p0fb_trkVeto/HTMHT_2012C_rr.root", false, false, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HTMHT_Run2012C-PromptReco-v2_AOD_UCSB1502ra2b_v66.root", "files12p0fb_trkVeto/HTMHT_2012C_pr.root", false, false, 1.0);

  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.HT_Run2012A-13Jul2012-v1_AOD_UCSB1491ra2b_v66.root", "files12p0fb_trkVeto/HT_2012A.root", false, false, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012B-13Jul2012-v1_AOD_UCSB1518ra2b_v66.root", "files12p0fb_trkVeto/JetHT_2012B.root", false, false, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012C-24Aug2012-v2_AOD_UCSB1510ra2b_v66.root", "files12p0fb_trkVeto/JetHT_2012C_rr.root", false, false, 1.0);
  makeTree("v66_reducedTrees_wIso/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.JetHT_Run2012C-PromptReco-v2_AOD_UCSB1511ra2b_v66.root", "files12p0fb_trkVeto/JetHT_2012C_pr.root", false, false, 1.0);
*/
 


  //makeTree("reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff04_HLTEff0.SMS-T1tttt_Mgluino-350to2000_mLSP-0to1650_8TeV-Pythia6Z_partial.root", "files12p0fb_lastHCP/T1tttt_partial.root", false, false, 1.0);

}


// If weight argument is <0, weight3 from reducedTree is used.  Lumi is hardcoded (see "double lumi = ").  
void makeTree(TString inputfile, TString outputfile, bool doQCD, bool blinding, float weight){
  TFile *f = new TFile(inputfile);
  TTree *t;
  f->GetObject("reducedTree",t);
  reducedTree rt(t);
  char* output4loop = outputfile;
  cout << " new output file name = " << output4loop << endl;
  
  if(weight<0) {
    double lumi = 17610.0; //number of inverse pb.
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
