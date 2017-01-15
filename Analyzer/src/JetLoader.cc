#include "../include/JetLoader.hh"
#include <cmath>
#include <iostream> 
#include <sstream>

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree) { 
  fJets  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("AK4Puppi",       &fJets);
  fJetBr = iTree->GetBranch("AK4Puppi");

  fJetsCHS  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("AK4CHS",       &fJetsCHS);
  fJetBrCHS = iTree->GetBranch("AK4CHS");

  fN = 4;
  fNV = 3; // max number of V jets to consider for dR anti-matching
  fNVars = 3; // pt, eta, phi
  fNOtherVars = 10; // Mass, b-tag, qgid, dR, dPhi, pt_cen, pt_JESUp, pt_JESDown, pt_JERUp, pt_JERDown
  r = new TRandom3(1988);
}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fJetBr;
  delete fJetsCHS;
  delete fJetBrCHS;
}
void JetLoader::reset() { 
  fNJetsPt30           = 0;
  fNFwdPt30            = 0;
  fNBTagsLPt30         = 0;
  fNBTagsMPt30         = 0;
  fNBTagsTPt30         = 0;
  fLooseJets.clear();
  fGoodJets.clear();
  selectedJets8.clear();
  selectedJets15.clear();  
  for(int i0 = 0; i0 < int(fNJetsPt30dR08.size()); i0++) {
    fNJetsPt30dR08[i0] = -999;
    fNBTagsLPt50dR08[i0] = -999;
    fNBTagsMPt50dR08[i0] = -999;
    fNBTagsTPt50dR08[i0] = -999;
    fNBTagsLPt100dR08[i0] = -999;
    fNBTagsMPt100dR08[i0] = -999;
    fNBTagsTPt100dR08[i0] = -999;
    fNBTagsLPt150dR08[i0] = -999;
    fNBTagsMPt150dR08[i0] = -999;
    fNBTagsTPt150dR08[i0] = -999;
  }  
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
}
void JetLoader::setupTree(TTree *iTree, std::string iJetLabel) { 
  reset();
  fTree = iTree;
  std::stringstream pSNPt30,pSfwdPt30,pSbPt30,pSbLPt30,pSbMPt30,pSbTPt30;
  pSNPt30     << "n" << iJetLabel << "sPt30";
  pSfwdPt30   << "n" << iJetLabel << "sfwdPt30";
  pSbPt30     << "n" << iJetLabel << "sbtagPt30";
  pSbLPt30    << "n" << iJetLabel << "sbtagLPt30";
  pSbMPt30    << "n" << iJetLabel << "sbtagMPt30";
  pSbTPt30    << "n" << iJetLabel << "sbtagTPt30";


  fTree->Branch(pSNPt30.str().c_str()           ,&fNJetsPt30           ,(pSNPt30.str()+"/I").c_str());  // jet multiplicity
  fTree->Branch(pSfwdPt30.str().c_str()         ,&fNFwdPt30            ,(pSfwdPt30.str()+"/I").c_str());
  fTree->Branch(pSbLPt30.str().c_str()          ,&fNBTagsLPt30         ,(pSbLPt30.str()+"/I").c_str()); // b tags
  fTree->Branch(pSbMPt30.str().c_str()          ,&fNBTagsMPt30         ,(pSbMPt30.str()+"/I").c_str());
  fTree->Branch(pSbTPt30.str().c_str()          ,&fNBTagsTPt30         ,(pSbTPt30.str()+"/I").c_str());


  fNJetsPt30dR08.clear();
  fNBTagsLPt50dR08.clear();
  fNBTagsMPt50dR08.clear();
  fNBTagsTPt50dR08.clear();
  fNBTagsLPt100dR08.clear();
  fNBTagsMPt100dR08.clear();
  fNBTagsTPt100dR08.clear();
  fNBTagsLPt150dR08.clear();
  fNBTagsMPt150dR08.clear();
  fNBTagsTPt150dR08.clear();
  for(int i0 = 0; i0 < fNV; i0++) {
    fNJetsPt30dR08.push_back(-999);
    fNBTagsLPt50dR08.push_back(-999);
    fNBTagsMPt50dR08.push_back(-999);
    fNBTagsTPt50dR08.push_back(-999);
    fNBTagsLPt100dR08.push_back(-999);
    fNBTagsMPt100dR08.push_back(-999);
    fNBTagsTPt100dR08.push_back(-999);
    fNBTagsLPt150dR08.push_back(-999);
    fNBTagsMPt150dR08.push_back(-999);
    fNBTagsTPt150dR08.push_back(-999);
  }  
  for(int i0 = 0; i0 < fNV; i0++) {
    std::stringstream pSNPt30dR08,pSbLPt50dR08,pSbMPt50dR08,pSbTPt50dR08,pSbLPt100dR08,pSbMPt100dR08,pSbTPt100dR08,pSbLPt150dR08,pSbMPt150dR08,pSbTPt150dR08; 
    pSNPt30dR08 << "n" << iJetLabel << "sPt30dR08_" << i0;
    pSbLPt50dR08<< "n" << iJetLabel << "sLPt50dR08_" << i0;
    pSbMPt50dR08<< "n" << iJetLabel << "sMPt50dR08_" << i0;
    pSbTPt50dR08<< "n" << iJetLabel << "sTPt50dR08_" << i0;
    pSbLPt100dR08<< "n" << iJetLabel << "sLPt100dR08_" << i0;
    pSbMPt100dR08<< "n" << iJetLabel << "sMPt100dR08_" << i0;
    pSbTPt100dR08<< "n" << iJetLabel << "sTPt100dR08_ " << i0;
    pSbLPt150dR08<< "n" << iJetLabel << "sLPt150dR08_" << i0;
    pSbMPt150dR08<< "n" << iJetLabel << "sMPt150dR08_" << i0;
    pSbTPt150dR08<< "n" << iJetLabel << "sTPt150dR08_" << i0;
    fTree->Branch(pSNPt30dR08.str().c_str()       ,&fNJetsPt30dR08[i0]       ,(pSNPt30dR08.str()+"/I").c_str());
    fTree->Branch(pSbLPt50dR08.str().c_str()      ,&fNBTagsLPt50dR08[i0]     ,(pSbLPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt50dR08.str().c_str()      ,&fNBTagsMPt50dR08[i0]     ,(pSbMPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt50dR08.str().c_str()      ,&fNBTagsTPt50dR08[i0]     ,(pSbTPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbLPt100dR08.str().c_str()      ,&fNBTagsLPt100dR08[i0]     ,(pSbLPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt100dR08.str().c_str()      ,&fNBTagsMPt100dR08[i0]     ,(pSbMPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt100dR08.str().c_str()      ,&fNBTagsTPt100dR08[i0]     ,(pSbTPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbLPt150dR08.str().c_str()      ,&fNBTagsLPt150dR08[i0]     ,(pSbLPt150dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt150dR08.str().c_str()      ,&fNBTagsMPt150dR08[i0]     ,(pSbMPt150dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt150dR08.str().c_str()      ,&fNBTagsTPt150dR08[i0]     ,(pSbTPt150dR08.str()+"/I").c_str());
  }

  for(int i0 = 0; i0 < fN*(10)+4; i0++) {double pVar = 0; fVars.push_back(pVar);}           
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars);                                            // from MonoXUtils.cc => fN=4 j*_pt,j*_eta,j*_phi for j1,j2,j3,j4 (3*4=12)
  addOthers  (iJetLabel.c_str(),iTree,fN,fVars);                                            // Mass, b-tag, qgid, dR, dPhi for j1,j2,j3,j4 (5*4=20)
}
void JetLoader::load(int iEvent) { 
  fJets   ->Clear();
  fJetBr ->GetEntry(iEvent);
}
void JetLoader::selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets, double iRho){
  reset(); 
  int lCountPt30 = 0, lNFwdPt30 = 0, lNBTagLPt30 = 0,lNBTagMPt30 = 0, lNBTagTPt30 = 0;
  
  for (int i1 = 0; i1 < int(iVJets.size()); i1++) {
    fNJetsPt30dR08[i1] = 0;
    fNBTagsLPt50dR08[i1] = 0;
    fNBTagsMPt50dR08[i1] = 0;
    fNBTagsTPt50dR08[i1] = 0;
    fNBTagsLPt100dR08[i1] = 0;
    fNBTagsMPt100dR08[i1] = 0;
    fNBTagsTPt100dR08[i1] = 0;
    fNBTagsLPt150dR08[i1] = 0;
    fNBTagsMPt150dR08[i1] = 0;
    fNBTagsTPt150dR08[i1] = 0;
  }
    
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(passVeto(pJet->eta,pJet->phi,0.4,iElectrons))                      continue;
    if(passVeto(pJet->eta,pJet->phi,0.4,iMuons))                          continue;
    if(passVeto(pJet->eta,pJet->phi,0.4,iPhotons))                        continue;
    if(pJet->pt        <=  30)                                            continue;
    if(fabs(pJet->eta) > 2.5 && fabs(pJet->eta) < 4.5) lNFwdPt30++;
    if(fabs(pJet->eta) >= 2.5)                                            continue;
    if(!passJetLooseSel(pJet))                                            continue;
    lCountPt30++;
    addJet(pJet,fLooseJets);

    TLorentzVector vPJet; vPJet.SetPtEtaPhiM(pJet->pt, pJet->eta, pJet->phi, pJet->mass);
    fGoodJets.push_back(pJet);
    
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVL){
      lNBTagLPt30++;
    }    
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVM){ 
      lNBTagMPt30++;
    }
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVT){ 
      lNBTagTPt30++;
    }
    
    // jet and b-tag multiplicity
    for (int i1 = 0; i1 < int(iVJets.size()); i1++) {
      
      if(iVJets[i1].Pt()>350 && vPJet.DeltaR(iVJets[i1])>0.8) {
	fNJetsPt30dR08[i1]++;
	if(pJet->pt>50 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVL) fNBTagsLPt50dR08[i1]++;
	if(pJet->pt>100 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVL) fNBTagsLPt100dR08[i1]++;
	if(pJet->pt>150 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVL) fNBTagsLPt150dR08[i1]++;
	
	if(pJet->pt>50 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVM) fNBTagsMPt50dR08[i1]++;
	if(pJet->pt>100 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVM) fNBTagsMPt100dR08[i1]++;
	if(pJet->pt>150 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVM) fNBTagsMPt150dR08[i1]++;
	
	if(pJet->pt>50 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVT) fNBTagsTPt50dR08[i1]++;
	if(pJet->pt>100 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVT) fNBTagsTPt100dR08[i1]++;
	if(pJet->pt>150 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVT) fNBTagsTPt150dR08[i1]++;
      }
    }
  }
  addVJet(fLooseJets,selectedJets);
  fNJetsPt30           = lCountPt30;
  fNFwdPt30            = lNFwdPt30;
  fNBTagsLPt30         = lNBTagLPt30;
  fNBTagsMPt30         = lNBTagMPt30;
  fNBTagsTPt30         = lNBTagTPt30;

  fillJetCorr(fN,fLooseJets,fVars,iRho);
  fillOthers(fN,fLooseJets,fVars,iVJets,iRho);
}
void JetLoader::fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho){ 
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) { 
    iVals[i0*3+0] = iObjects[i0]->pt;
    iVals[i0*3+1] = iObjects[i0]->eta;
    iVals[i0*3+2] = iObjects[i0]->phi;
  }
}
void JetLoader::addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals) { 
  for(int i0 = 0; i0 < iN; i0++) { 
    int lBase = iN*fNVars+i0*fNOtherVars;
    std::stringstream pSMass,pSCSV,pSQGID,pSdR,pSdP,pScen,pSjesUp,pSjesDown,pSjerUp,pSjerDown;
    pSMass  << iHeader << i0 << "_mass";
    pSCSV   << iHeader << i0 << "_csv";
    pSQGID  << iHeader << i0 << "_qgid";
    pSdR    << iHeader << i0 << "_dR08";
    pSdP    << iHeader << i0 << "_dPhi08";
    pScen << iHeader << i0 << "_pt_cen";
    pSjesUp << iHeader << i0 << "_pt_JESUp";
    pSjesDown << iHeader << i0 << "_pt_JESDown";
    pSjerUp << iHeader << i0 << "_pt_JERUp";
    pSjerDown << iHeader << i0 << "_pt_JERDown";
    
    iTree->Branch(pSMass .str().c_str(),&iVals[lBase+0],(pSMass .str()+"/D").c_str());
    iTree->Branch(pSCSV .str().c_str() ,&iVals[lBase+1],(pSCSV  .str()+"/D").c_str());
    iTree->Branch(pSQGID.str().c_str() ,&iVals[lBase+2],(pSQGID .str()+"/D").c_str());
    iTree->Branch(pSdR  .str().c_str() ,&iVals[lBase+3],(pSdR   .str()+"/D").c_str());
    iTree->Branch(pSdP  .str().c_str() ,&iVals[lBase+4],(pSdP   .str()+"/D").c_str());
    iTree->Branch(pScen  .str().c_str() ,&iVals[lBase+5],(pScen   .str()+"/D").c_str());
    iTree->Branch(pSjesUp  .str().c_str() ,&iVals[lBase+6],(pSjesUp   .str()+"/D").c_str());
    iTree->Branch(pSjesDown.str().c_str() ,&iVals[lBase+7],(pSjesDown .str()+"/D").c_str());
    iTree->Branch(pSjerUp  .str().c_str() ,&iVals[lBase+8],(pSjerUp   .str()+"/D").c_str());
    iTree->Branch(pSjerDown.str().c_str() ,&iVals[lBase+9],(pSjerDown .str()+"/D").c_str());
  }
}
void JetLoader::fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, std::vector<TLorentzVector> iVJets, double iRho){ 
  int lBase = fNVars*fN;
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) {
    TLorentzVector vPJet; vPJet.SetPtEtaPhiM(iObjects[i0]->pt, iObjects[i0]->eta, iObjects[i0]->phi, iObjects[i0]->mass); 
    iVals[lBase+i0*fNOtherVars+0] = iObjects[i0]->mass;
    iVals[lBase+i0*fNOtherVars+1] = iObjects[i0]->csv;
    iVals[lBase+i0*fNOtherVars+2] = iObjects[i0]->qgid;
    if(iVJets.size()>0) {
      iVals[lBase+i0*fNOtherVars+3] = vPJet.DeltaR(iVJets[0]);
      iVals[lBase+i0*fNOtherVars+4] = vPJet.DeltaPhi(iVJets[0]);
    }
    
    double x1 = r->Gaus();
    double x2 = r->Gaus();
    double x3 = r->Gaus();
    double sf = getJerSF(iObjects[i0]->eta,0);
    double sfUp = getJerSF(iObjects[i0]->eta,1);
    double sfDown = getJerSF(iObjects[i0]->eta,-1);
    double sigma_MC = 0.07; // assume 7% for now
    double jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    double jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x2;
    double jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x3;

    double unc = iObjects[i0]->unc;
    double jetCorrPt = iObjects[i0]->pt;
    double jetCorrPtSmear = (iObjects[i0]->pt)*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;
    
    //double jetE = vPJet.E()  
    //double JEC = JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
    //					   iRho, iObjects[i0]->area, 
    //					   runNum,
    //					   JetCorrectorIOV,JetCorrector);
    //double unc = getJecUnc( jetCorrPt, iObjects[i0]->eta, runNum ); //use run=999 as default
    
    iVals[lBase+i0*fNOtherVars+5] = jetCorrPtSmear;
    iVals[lBase+i0*fNOtherVars+6] = jetPtJESUp;
    iVals[lBase+i0*fNOtherVars+7] = jetPtJESDown;
    iVals[lBase+i0*fNOtherVars+8] = jetPtJERUp;    
    iVals[lBase+i0*fNOtherVars+9] = jetPtJERDown;
      
  }
}
//2016 Prompt Reco
void JetLoader::loadJECs(bool isData) {
    std::cout << "JetLoader: loading jet energy correction constants" << std::endl;
    // initialize
    std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    if (isData) {      
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_MC/Spring16_25nsV6_MC_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }

}
void JetLoader::loadJECs_Rereco(bool isData) {
    // initialize
    std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
 
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016BCDV2_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016BCDV2_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016BCDV2_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016BCDV2_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersBCD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Spring16_23Sep2016BCDV2_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);
      SimpleJetResolution* JetResolutionCalculatorBCD = new SimpleJetResolution(*JetResolutionParametersBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetResolutionParameters.push_back(JetResolutionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorBCD);
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016EFV2_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016EFV2_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016EFV2_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016EFV2_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersEF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      std::string jecUncPathEF = jecPathname+"/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016EFV2_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(jecUncPathEF);
      SimpleJetResolution* JetResolutionCalculatorEF = new SimpleJetResolution(*JetResolutionParametersEF);

      correctionParameters.push_back(correctionParametersEF);
      JetResolutionParameters.push_back(JetResolutionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorEF);
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278801 ));

      //IOV: 2016G
      std::vector<JetCorrectorParameters> correctionParametersG = std::vector<JetCorrectorParameters> ();
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016GV2_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016GV2_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016GV2_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016GV2_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersG = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorG = new FactorizedJetCorrector(correctionParametersG);
      std::string jecUncPathG = jecPathname+"/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016GV2_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncG = new JetCorrectionUncertainty(jecUncPathG);
      SimpleJetResolution* JetResolutionCalculatorG = new SimpleJetResolution(*JetResolutionParametersG);

      correctionParameters.push_back(correctionParametersG);
      JetResolutionParameters.push_back(JetResolutionParametersG);
      JetCorrector.push_back( JetCorrectorG );
      JetResolutionCalculator.push_back(JetResolutionCalculatorG);
      jecUnc.push_back(jecUncG);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 280385 ));

      //IOV: 2016H
      std::vector<JetCorrectorParameters> correctionParametersH = std::vector<JetCorrectorParameters> ();
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016HV2_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016HV2_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016HV2_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016HV2_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersH = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorH = new FactorizedJetCorrector(correctionParametersH);
      std::string jecUncPathH = jecPathname+"/Spring16_23Sep2016_V2_DATA/Spring16_23Sep2016HV2_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncH = new JetCorrectionUncertainty(jecUncPathH);
      SimpleJetResolution* JetResolutionCalculatorH = new SimpleJetResolution(*JetResolutionParametersH);

      correctionParameters.push_back(correctionParametersH);
      JetResolutionParameters.push_back(JetResolutionParametersH);
      JetCorrector.push_back( JetCorrectorH );
      JetResolutionCalculator.push_back(JetResolutionCalculatorH);
      jecUnc.push_back(jecUncH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 280919, 99999999 ));

    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_MC/Spring16_23Sep2016V2_MC_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_MC/Spring16_23Sep2016V2_MC_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V2_MC/Spring16_23Sep2016V2_MC_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Spring16_23Sep2016_V2_MC/Spring16_23Sep2016V2_MC_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }
  
}

void JetLoader::loadCMSSWPath() {
    char* cmsswPathChar = getenv("CMSSW_BASE");
    if (cmsswPathChar == NULL) {
        std::cout << "Warning in JetLoader::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
        cmsswPath = "";
    }
    cmsswPath = std::string(cmsswPathChar);
}


// Retrieve jet energy uncertainty as a function of pt and eta
double JetLoader::getJecUnc( float pt, float eta , int run) {

  int foundIndex = -1;
  for (uint i=0; i<JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    std::cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
    foundIndex = 0;
  }

  jecUnc[foundIndex]->setJetPt(pt);
  jecUnc[foundIndex]->setJetEta(eta);
  return jecUnc[foundIndex]->getUncertainty(true);
}


// Retrieve jet energy resolution scale factor as a function of eta
double JetLoader::getJerSF( float eta, int nsigma) {
  double sf[3] = {1.0, 1.0, 1.0};
  if (fabs(eta) < 0.5) {
    sf[0] = 1.109; 
    sf[1] = 1.101; 
    sf[2] = 1.117;
  }
  else if (fabs(eta) >= 0.5 && fabs(eta) < 0.8) {
    sf[0] = 1.138;
    sf[1] = 1.125;
    sf[2] = 1.151;
  }
  else if (fabs(eta) >= 0.8 && fabs(eta) < 1.1) {
    sf[0] = 1.114;
    sf[1] = 1.101; 
    sf[2] = 1.127;
  }
  else if (fabs(eta) >= 1.1 && fabs(eta) < 1.3) {
    sf[0] = 1.123;
    sf[1] = 1.099;
    sf[2] = 1.147;
  }
  else if (fabs(eta) >= 1.3 && fabs(eta) < 1.7) {
    sf[0] = 1.084;
    sf[1] = 1.073;
    sf[2] = 1.095;
  }
  else if (fabs(eta) >= 1.7 && fabs(eta) < 1.9) {
    sf[0] = 1.082;
    sf[1] = 1.047;
    sf[2] = 1.117;
  }
  else if (fabs(eta) >= 1.9 && fabs(eta) < 2.1) {
    sf[0] = 1.140;
    sf[1] = 1.093;
    sf[2] = 1.187;
  }
  else if (fabs(eta) >= 2.1 && fabs(eta) < 2.3) {
    sf[0] = 1.067;
    sf[1] = 1.014;
    sf[2] = 1.120;
  }
  else if (fabs(eta) >= 2.3 && fabs(eta) < 2.5) {
    sf[0] = 1.177;
    sf[1] = 1.136;
    sf[2] = 1.218;
  }
  else if (fabs(eta) >= 2.5 && fabs(eta) < 2.8) {
    sf[0] = 1.364;
    sf[1] = 1.325;
    sf[2] = 1.403;
  }
  else if (fabs(eta) >= 2.8 && fabs(eta) < 3.0) {
    sf[0] = 1.857;
    sf[1] = 1.786;
    sf[2] = 1.928;
  }
  else if (fabs(eta) >= 3.0 && fabs(eta) < 3.2) {
    sf[0] = 1.328;
    sf[1] = 1.306;
    sf[2] = 1.350;
  }
  else if (fabs(eta) >= 3.2 && fabs(eta) < 4.7) {
    sf[0] = 1.160;
    sf[1] = 1.131;
    sf[2] = 1.189;
  }

  if (nsigma == 0) {
    return sf[0];
  }
  else if (nsigma == -1) {
    return sf[1];
  }
  else if (nsigma == 1) {
    return sf[2];
  }
  // if nsigma ! = 0, -1, 1
  std::cout << "Warning: nsigma = " << nsigma << " is not -1, 0, 1" << std::endl;
  return 1;
}

