#include "../include/JetLoader.hh"
#include <cmath>
#include <iostream> 
#include <sstream>

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree, bool iData) { 
  fJets  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("AK4Puppi",       &fJets);
  fJetBr = iTree->GetBranch("AK4Puppi");

  fJetsCHS  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("AK4CHS",       &fJetsCHS);
  fJetBrCHS = iTree->GetBranch("AK4CHS");

  fN = 4;
  fNV = 3; // max number of V jets to consider for dR anti-matching
  fNVars = 3; // pt, eta, phi
  fNOtherVars = 10; // Mass, b-tag, qgid, dR, dPhi, pt_old, pt_JESUp, pt_JESDown, pt_JERUp, pt_JERDown
  isData = iData;
  loadJECs_Rereco(isData);    
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
  fNJetsPt30jesUp      = 0;
  fNJetsPt30jesDown    = 0;
  fNJetsPt30jerUp      = 0;
  fNJetsPt30jerDown    = 0;
  fNFwdPt30            = 0;
  fNBTagsLPt30         = 0;
  fNBTagsMPt30         = 0;
  fNBTagsTPt30         = 0;
  MetXCorrjesUp        = 0;
  MetYCorrjesUp        = 0;
  MetXCorrjesDown      = 0;
  MetYCorrjesDown      = 0;
  MetXCorrjerUp        = 0;
  MetYCorrjerUp        = 0;
  MetXCorrjerDown      = 0;
  MetYCorrjerDown      = 0;
  fLooseJets.clear();
  fGoodJets.clear();
  selectedJets8.clear();
  selectedJets15.clear();
  x1List.clear();
  x2List.clear();
  x3List.clear();    
  for(int i0 = 0; i0 < int(fNJetsPt30dR08.size()); i0++) {
    fNJetsPt30dR08[i0] = -999;
    fNJetsPt30dR08jesUp[i0] = -999;
    fNJetsPt30dR08jesDown[i0] = -999;
    fNJetsPt30dR08jerUp[i0] = -999;
    fNJetsPt30dR08jerDown[i0] = -999;
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
  std::stringstream pSNPt30jesUp,pSNPt30jesDown,pSNPt30jerUp,pSNPt30jerDown;
  std::stringstream pSMetXCorrjesUp,pSMetYCorrjesUp;
  std::stringstream pSMetXCorrjesDown,pSMetYCorrjesDown;
  std::stringstream pSMetXCorrjerUp,pSMetYCorrjerUp;
  std::stringstream pSMetXCorrjerDown,pSMetYCorrjerDown;
  
  pSNPt30     << "n" << iJetLabel << "sPt30";
  pSNPt30jesUp << "n" << iJetLabel << "sPt30jesUp";
  pSNPt30jesDown << "n" << iJetLabel << "sPt30jesDown";
  pSNPt30jerUp << "n" << iJetLabel << "sPt30jerUp";
  pSNPt30jerDown << "n" << iJetLabel << "sPt30jerDown";
  pSfwdPt30   << "n" << iJetLabel << "sfwdPt30";
  pSbPt30     << "n" << iJetLabel << "sbtagPt30";
  pSbLPt30    << "n" << iJetLabel << "sbtagLPt30";
  pSbMPt30    << "n" << iJetLabel << "sbtagMPt30";
  pSbTPt30    << "n" << iJetLabel << "sbtagTPt30";
  
  pSMetXCorrjesUp     << "MetXCorrjesUp";
  pSMetYCorrjesUp     << "MetYCorrjesUp";
  pSMetXCorrjesDown     << "MetXCorrjesDown";
  pSMetYCorrjesDown     << "MetYCorrjesDown";
  pSMetXCorrjerUp     << "MetXCorrjerUp";
  pSMetYCorrjerUp     << "MetYCorrjerUp";
  pSMetXCorrjerDown     << "MetXCorrjerDown";
  pSMetYCorrjerDown     << "MetYCorrjerDown";

  fTree->Branch(pSNPt30.str().c_str()           ,&fNJetsPt30           ,(pSNPt30.str()+"/I").c_str());  // jet multiplicity
  fTree->Branch(pSNPt30jesUp.str().c_str()           ,&fNJetsPt30jesUp           ,(pSNPt30jesUp.str()+"/I").c_str());  
  fTree->Branch(pSNPt30jesDown.str().c_str()           ,&fNJetsPt30jesDown           ,(pSNPt30jesDown.str()+"/I").c_str());
  fTree->Branch(pSNPt30jerUp.str().c_str()           ,&fNJetsPt30jerUp           ,(pSNPt30jerUp.str()+"/I").c_str());  
  fTree->Branch(pSNPt30jerDown.str().c_str()           ,&fNJetsPt30jerDown           ,(pSNPt30jerDown.str()+"/I").c_str());
  fTree->Branch(pSfwdPt30.str().c_str()         ,&fNFwdPt30            ,(pSfwdPt30.str()+"/I").c_str());
  fTree->Branch(pSbLPt30.str().c_str()          ,&fNBTagsLPt30         ,(pSbLPt30.str()+"/I").c_str()); // b tags
  fTree->Branch(pSbMPt30.str().c_str()          ,&fNBTagsMPt30         ,(pSbMPt30.str()+"/I").c_str());
  fTree->Branch(pSbTPt30.str().c_str()          ,&fNBTagsTPt30         ,(pSbTPt30.str()+"/I").c_str());
  
  fTree->Branch(pSMetXCorrjesUp.str().c_str()           ,&MetXCorrjesUp           ,(pSMetXCorrjesUp.str()+"/D").c_str());  // Met corrections
  fTree->Branch(pSMetYCorrjesUp.str().c_str()           ,&MetYCorrjesUp           ,(pSMetYCorrjesUp.str()+"/D").c_str());  
  fTree->Branch(pSMetXCorrjesDown.str().c_str()           ,&MetXCorrjesDown           ,(pSMetXCorrjesDown.str()+"/D").c_str());  
  fTree->Branch(pSMetYCorrjesDown.str().c_str()           ,&MetYCorrjesDown           ,(pSMetYCorrjesDown.str()+"/D").c_str());  
  fTree->Branch(pSMetXCorrjerUp.str().c_str()           ,&MetXCorrjerUp           ,(pSMetXCorrjerUp.str()+"/D").c_str());  
  fTree->Branch(pSMetYCorrjerUp.str().c_str()           ,&MetYCorrjerUp           ,(pSMetYCorrjerUp.str()+"/D").c_str());  
  fTree->Branch(pSMetXCorrjerDown.str().c_str()           ,&MetXCorrjerDown           ,(pSMetXCorrjerDown.str()+"/D").c_str());  
  fTree->Branch(pSMetYCorrjerDown.str().c_str()           ,&MetYCorrjerDown           ,(pSMetYCorrjerDown.str()+"/D").c_str());  



  fNJetsPt30dR08.clear();
  fNJetsPt30dR08jesUp.clear();
  fNJetsPt30dR08jesDown.clear();
  fNJetsPt30dR08jerUp.clear();
  fNJetsPt30dR08jerDown.clear();  
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
    fNJetsPt30dR08jesUp.push_back(-999);
    fNJetsPt30dR08jesDown.push_back(-999);
    fNJetsPt30dR08jerUp.push_back(-999);
    fNJetsPt30dR08jerDown.push_back(-999);
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
    std::stringstream pSNPt30dR08jesUp,pSNPt30dR08jesDown,pSNPt30dR08jerUp,pSNPt30dR08jerDown;
    pSNPt30dR08 << "n" << iJetLabel << "sPt30dR08_" << i0;
    pSNPt30dR08jesUp << "n" << iJetLabel << "sPt30dR08jesUp_" << i0;
    pSNPt30dR08jesDown << "n" << iJetLabel << "sPt30dR08jesDown_" << i0;
    pSNPt30dR08jerUp << "n" << iJetLabel << "sPt30dR08jerUp_" << i0;
    pSNPt30dR08jerDown << "n" << iJetLabel << "sPt30dR08jerDown_" << i0;
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
    fTree->Branch(pSNPt30dR08jesUp.str().c_str()       ,&fNJetsPt30dR08jesUp[i0]       ,(pSNPt30dR08jesUp.str()+"/I").c_str());
    fTree->Branch(pSNPt30dR08jesDown.str().c_str()       ,&fNJetsPt30dR08jesDown[i0]       ,(pSNPt30dR08jesDown.str()+"/I").c_str());
    fTree->Branch(pSNPt30dR08jerUp.str().c_str()       ,&fNJetsPt30dR08jerUp[i0]       ,(pSNPt30dR08jerUp.str()+"/I").c_str());
    fTree->Branch(pSNPt30dR08jerDown.str().c_str()       ,&fNJetsPt30dR08jerDown[i0]       ,(pSNPt30dR08jerDown.str()+"/I").c_str());
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
void JetLoader::selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets, double iRho, unsigned int runNum){
  reset(); 
  int lCountPt30 = 0, lNFwdPt30 = 0, lNBTagLPt30 = 0,lNBTagMPt30 = 0, lNBTagTPt30 = 0;
  int lCountPt30jesUp = 0, lCountPt30jesDown = 0, lCountPt30jerUp = 0, lCountPt30jerDown = 0;
  
  for (int i1 = 0; i1 < int(iVJets.size()); i1++) {
    fNJetsPt30dR08[i1] = 0;
    fNJetsPt30dR08jesUp[i1] = 0;
    fNJetsPt30dR08jesDown[i1] = 0;
    fNJetsPt30dR08jerUp[i1] = 0;
    fNJetsPt30dR08jerDown[i1] = 0;
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
    
    double JEC_old = (pJet->pt)/(pJet->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(pJet->ptRaw, pJet->eta, pJet->phi, (pJet->mass)/JEC_old);
    double jetE = vPJet.E();
    double JEC = JetEnergyCorrectionFactor(pJet->ptRaw, pJet->eta, pJet->phi, jetE, 
    					   iRho, pJet->area, 
    					   runNum,
   					   JetCorrectionsIOV,JetCorrector);
    double jetCorrPt = JEC*(pJet->ptRaw);
    double jetCorrE = JEC*(vPJet.E());
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, pJet->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; // max 44.30 for Spring16_25nsV6_MC JER (CHANGE ONCE UPDATED)
    float sigma_MC = resolution.getResolution(parameters);
    float sf = resolution_sf.getScaleFactor(parameters);
    float sfUp = resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = resolution_sf.getScaleFactor(parameters, Variation::DOWN);
    double x1 = r->Gaus();
    double x2 = r->Gaus();
    double x3 = r->Gaus();
    double jetEnergySmearFactor = 1.0; 
    double jetEnergySmearFactorUp = 1.0; 
    double jetEnergySmearFactorDown = 1.0;    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
      jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x2;
      jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x3;
    }    
    double unc = getJecUnc( jetCorrPt, pJet->eta, runNum ); //use run=999 as default
    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;
    
    double jetCorrESmear = jetCorrE*jetEnergySmearFactor*(1+unc);    
    double jetEJESUp = jetCorrE*jetEnergySmearFactor*(1+unc);
    double jetEJESDown = jetCorrE*jetEnergySmearFactor/(1+unc);
    double jetEJERUp = jetCorrE*jetEnergySmearFactorUp;
    double jetEJERDown = jetCorrE*jetEnergySmearFactorDown;
    
    TLorentzVector thisJet;  thisJet.SetPtEtaPhiE(jetCorrPtSmear, pJet->eta, pJet->phi, jetCorrESmear);
    TLorentzVector thisJetJESUp;  thisJetJESUp.SetPtEtaPhiE(jetPtJESUp, pJet->eta, pJet->phi, jetEJESUp);
    TLorentzVector thisJetJESDown; thisJetJESDown.SetPtEtaPhiE(jetPtJESDown, pJet->eta, pJet->phi, jetEJESDown);
    TLorentzVector thisJetJERUp;  thisJetJERUp.SetPtEtaPhiE(jetPtJERUp,  pJet->eta, pJet->phi, jetEJERUp);
    TLorentzVector thisJetJERDown; thisJetJERDown.SetPtEtaPhiE(jetPtJERDown,  pJet->eta, pJet->phi, jetEJERDown);
    
    //Propagate uncertainties to the MET
    if (jetPtJESUp > 20) {
      MetXCorrjesUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
      MetYCorrjesUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
    }
    if (jetPtJESDown > 20) {
      MetXCorrjesDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
      MetYCorrjesDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
    }
    if (jetPtJERUp > 20) {
      MetXCorrjerUp += -1 * (thisJetJERUp.Px() - thisJet.Px());
      MetYCorrjerUp += -1 * (thisJetJERUp.Py() - thisJet.Py());
    }
    if (jetPtJERDown > 20) {
      MetXCorrjerDown += -1 * (thisJetJERDown.Px() - thisJet.Px());
      MetYCorrjerDown += -1 * (thisJetJERDown.Py() - thisJet.Py());
    }    

    if (jetCorrPtSmear  > 30) lCountPt30++;
    if (jetPtJESUp  > 30) lCountPt30jesUp++;
    if (jetPtJESDown  > 30) lCountPt30jesDown++;
    if (jetPtJERUp  > 30) lCountPt30jerUp++;
    if (jetPtJERDown  > 30) lCountPt30jerDown++;
    
    // jet and b-tag multiplicity
    
    vPJet.SetPtEtaPhiE(jetCorrPtSmear, pJet->eta, pJet->phi, jetCorrESmear);
    for (int i1 = 0; i1 < int(iVJets.size()); i1++) {      
      if(iVJets[i1].Pt()>350 && vPJet.DeltaR(iVJets[i1])>0.8) {
	if (jetCorrPtSmear  > 30) fNJetsPt30dR08[i1]++;
	if (jetPtJESUp  > 30) fNJetsPt30dR08jesUp[i1]++;
	if (jetPtJESDown  > 30) fNJetsPt30dR08jesDown[i1]++;
	if (jetPtJESUp  > 30) fNJetsPt30dR08jesUp[i1]++;
	if (jetPtJESDown  > 30) fNJetsPt30dR08jesDown[i1]++;
      }
    }
	
    
    if(jetCorrPtSmear  <=  30)                                            continue;    
    if(fabs(pJet->eta) > 2.5 && fabs(pJet->eta) < 4.5) lNFwdPt30++;
    if(fabs(pJet->eta) >= 2.5)                                            continue;
    if(!passJetLooseSel(pJet))                                            continue;
    x1List.push_back(x1);
    x2List.push_back(x2);
    x3List.push_back(x3);
    addJet(pJet,fLooseJets);
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
  fNJetsPt30jesUp      = lCountPt30jesUp;
  fNJetsPt30jesDown    = lCountPt30jesDown;
  fNJetsPt30jerUp      = lCountPt30jerUp;
  fNJetsPt30jerDown    = lCountPt30jerDown;
  fNFwdPt30            = lNFwdPt30;
  fNBTagsLPt30         = lNBTagLPt30;
  fNBTagsMPt30         = lNBTagMPt30;
  fNBTagsTPt30         = lNBTagTPt30;

  fillJetCorr(fN,fLooseJets,fVars,iRho,runNum);
  fillOthers(fN,fLooseJets,fVars,iVJets,iRho, runNum);
}
void JetLoader::fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum){ 
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  
  TLorentzVector vPJet;
  
  
  for(int i0 = 0; i0 < lMin; i0++) {
    
    double JEC_old = (iObjects[i0]->pt)/(iObjects[i0]->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, (iObjects[i0]->mass)/JEC_old);
    double jetE = vPJet.E();
    double JEC = JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
    					   iRho, iObjects[i0]->area, 
    					   runNum,
   					   JetCorrectionsIOV,JetCorrector);    
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    double unc = getJecUnc( jetCorrPt, iObjects[i0]->eta, runNum ); //use run=999 as default
    double x1 = x1List[i0];
    double jetEnergySmearFactor = 1.0;    
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; // max 44.30 for Spring16_25nsV6_MC JER (CHANGE ONCE UPDATED)
    float sigma_MC = resolution.getResolution(parameters);
    float sf = resolution_sf.getScaleFactor(parameters);    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    }    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    iVals[i0*3+0] = jetCorrPtSmear;
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
    pScen << iHeader << i0 << "_pt_old";
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
void JetLoader::fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, std::vector<TLorentzVector> iVJets, double iRho, unsigned int runNum){ 
  int lBase = fNVars*fN;
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) {
            
    double JEC_old = (iObjects[i0]->pt)/(iObjects[i0]->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, (iObjects[i0]->mass)/JEC_old);
    double jetE = vPJet.E();

    double JEC = JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
    					   iRho, iObjects[i0]->area, 
    					   runNum,
   					   JetCorrectionsIOV,JetCorrector);
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    
    double unc_old = iObjects[i0]->unc;
    double unc = getJecUnc( jetCorrPt, iObjects[i0]->eta, runNum ); //use run=999 as default
    
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}};
    float sigma_MC = resolution.getResolution(parameters);
    float sf = resolution_sf.getScaleFactor(parameters);
    float sfUp = resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = resolution_sf.getScaleFactor(parameters, Variation::DOWN);

    
    double x1 = x1List[i0];
    double x2 = x2List[i0];
    double x3 = x3List[i0];
    double jetEnergySmearFactor = 1.0; 
    double jetEnergySmearFactorUp = 1.0; 
    double jetEnergySmearFactorDown = 1.0;    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
      jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x2;
      jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x3;
    }    
    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;
    /*
    if (true){
      std::cout << "Jet" << std::endl;
      std::cout << "i0 =" << i0 << std::endl;
      std::cout << "ptraw = " << iObjects[i0]->ptRaw << std::endl;
      std::cout << "eta = " << iObjects[i0]->eta << std::endl;
      std::cout << "rho = " << iRho << std::endl;
      std::cout << "runNum = " << runNum << std::endl;
      std::cout << "sf = " << sf << ", " <<  sfUp << ", " << sfDown << std::endl;
      std::cout << "sigma_MC = " << sigma_MC << std::endl;    
      std::cout << "runNum = " << runNum << std::endl;
      std::cout << "unc_old = " << unc_old << std::endl;
      std::cout << "unc = " << unc << std::endl;
      std::cout << "JEC_old = " << JEC_old << std::endl;
      std::cout << "JEC = " << JEC << std::endl;
      std::cout << "x1 = " << x1 << std::endl;
      std::cout << "x2 = " << x2 << std::endl;
      std::cout << "x3 = " << x3 << std::endl;
      std::cout << "ptcorr = " << jetCorrPt << std::endl;
      std::cout << "ptcorrsmear = " << jetCorrPtSmear << std::endl;
      std::cout << "jesup = " << jetPtJESUp << std::endl;
      std::cout << "jesdown = " << jetPtJESDown << std::endl;
      std::cout << "jerup = " << jetPtJERUp << std::endl;
      std::cout << "jerdown = " << jetPtJERDown << std::endl;
    }
    */
    iVals[lBase+i0*fNOtherVars+0] = JEC*jetEnergySmearFactor*(iObjects[i0]->mass);
    iVals[lBase+i0*fNOtherVars+1] = iObjects[i0]->csv;
    iVals[lBase+i0*fNOtherVars+2] = iObjects[i0]->qgid;
    if(iVJets.size()>0) {
      iVals[lBase+i0*fNOtherVars+3] = vPJet.DeltaR(iVJets[0]);
      iVals[lBase+i0*fNOtherVars+4] = vPJet.DeltaPhi(iVJets[0]);
    }    
    iVals[lBase+i0*fNOtherVars+5] = iObjects[i0]->pt;
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
    loadCMSSWPath();
    std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    resolution = JME::JetResolution(Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_PtResolution_AK4PFPuppi.txt",jecPathname.c_str()));
    resolution_sf = JME::JetResolutionScaleFactor(Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_SF_AK4PFPuppi.txt",jecPathname.c_str()));
    
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
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
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

      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_MC/Spring16_25nsV6_MC_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }

}
void JetLoader::loadJECs_Rereco(bool isData) {
    std::cout << "JetLoader: loading Rereco jet energy correction constants" << std::endl;
    // initialize
    loadCMSSWPath();
    std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    resolution = JME::JetResolution(Form("%s/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt",jecPathname.c_str()));
    resolution_sf = JME::JetResolutionScaleFactor(Form("%s/Spring16_25nsV10_MC/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt",jecPathname.c_str()));
 
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016EF
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      std::string jecUncPathEF = jecPathname+"/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(jecUncPathEF);

      correctionParameters.push_back(correctionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278801 ));

      //IOV: 2016G
      std::vector<JetCorrectorParameters> correctionParametersG = std::vector<JetCorrectorParameters> ();
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorG = new FactorizedJetCorrector(correctionParametersG);
      std::string jecUncPathG = jecPathname+"/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncG = new JetCorrectionUncertainty(jecUncPathG);

      correctionParameters.push_back(correctionParametersG);
      JetCorrector.push_back( JetCorrectorG );
      jecUnc.push_back(jecUncG);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 280385 ));

      //IOV: 2016H
      std::vector<JetCorrectorParameters> correctionParametersH = std::vector<JetCorrectorParameters> ();
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorH = new FactorizedJetCorrector(correctionParametersH);
      std::string jecUncPathH = jecPathname+"/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncH = new JetCorrectionUncertainty(jecUncPathH);

      correctionParameters.push_back(correctionParametersH);
      JetCorrector.push_back( JetCorrectorH );
      jecUnc.push_back(jecUncH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 280919, 99999999 ));

    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFPuppi.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
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


//Jet Energy Corrections
double JetLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
						 double rho, double jetArea,
						 int run,
						 std::vector<std::pair<int,int> > JetCorrectionsIOV,
						 std::vector<FactorizedJetCorrector*> jetcorrector,
						 int jetCorrectionLevel,
						 bool printDebug) {

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

  if (!jetcorrector[foundIndex]) {
    std::cout << "WWARNING: Jet corrector pointer is null. Returning JEC = 0. \n";
    return 0;
  }

  jetcorrector[foundIndex]->setJetEta(jetEta);
  jetcorrector[foundIndex]->setJetPt(jetRawPt);
  jetcorrector[foundIndex]->setJetPhi(jetPhi);
  jetcorrector[foundIndex]->setJetE(jetE);
  jetcorrector[foundIndex]->setRho(rho);
  jetcorrector[foundIndex]->setJetA(jetArea);

  std::vector<float> corrections;
  corrections = jetcorrector[foundIndex]->getSubCorrections();

  if (printDebug) std::cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jetRawPt << " " << jetEta << " " << jetPhi << "\n";

  double cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {

    //only correct up to the required level. if -1, then do all correction levels
    if (jetCorrectionLevel >= 0 && int(j) > jetCorrectionLevel) continue;

    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) std::cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) std::cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";
  
  return cumulativeCorrection;

}

//Jet Energy Corrections
double JetLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
						 double rho, double jetArea,
						 FactorizedJetCorrector *jetcorrector,
						 int jetCorrectionLevel,
						 bool printDebug) {
  if (!jetcorrector) {
    std::cout << "WWARNING: Jet corrector pointer is null. Returning JEC = 0. \n";
    return 0;
  }

  jetcorrector->setJetEta(jetEta);
  jetcorrector->setJetPt(jetRawPt);
  jetcorrector->setJetPhi(jetPhi);
  jetcorrector->setJetE(jetE);
  jetcorrector->setRho(rho);
  jetcorrector->setJetA(jetArea);

  std::vector<float> corrections;
  corrections = jetcorrector->getSubCorrections();

  if (printDebug) std::cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jetRawPt << " " << jetEta << " " << jetPhi << "\n";

  double cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {

    //only correct up to the required level. if -1, then do all correction levels
    if (jetCorrectionLevel >= 0 && int(j) > jetCorrectionLevel) continue;

    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) std::cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) std::cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";
  
  return cumulativeCorrection;

}

