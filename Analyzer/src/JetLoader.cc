#include "../include/JetLoader.hh"
#include <cmath>
#include <iostream> 
#include <sstream>

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree, bool iData, std::string iLabel) {
  fJets  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("AK4Puppi",       &fJets);
  fJetBr = iTree->GetBranch("AK4Puppi");

  fN = 6; // number of jets to save
  fNV = 4; // max number of V jets to consider for dR anti-matching
  fNVars = 3; // pt, eta, phi
  fNOtherVars = 20; // Mass, b-tag, qgid, dR, dPhi, pt_old, pt_JESUp, pt_JESDown, pt_JERUp, pt_JERDown + 8 DeepCSV vars + flavors

  fJEC = new JECLoader(iData,iLabel,"AK4PFPuppi");
  r = new TRandom3(1988);
  isData=iData;
  fYear=iLabel;

  // modify B-tag wp by year
  if(iLabel==label2016){
    DEEPCSVL = fdeepCSVL2016; 
    DEEPCSVM = fdeepCSVM2016;
    DEEPCSVT = fdeepCSVT2016;
  }
  else if(iLabel==label2017){
    DEEPCSVL = fdeepCSVL2017;
    DEEPCSVM = fdeepCSVM2017;
    DEEPCSVT = fdeepCSVT2017;
  }
  else if(iLabel==label2018){
    DEEPCSVL = fdeepCSVL2018;
    DEEPCSVM = fdeepCSVM2018;
    DEEPCSVT = fdeepCSVT2018;
  }
  else {
    DEEPCSVL = fdeepCSVL2017;
    DEEPCSVM = fdeepCSVM2017;
    DEEPCSVT = fdeepCSVT2017;
  }

}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fJetBr;
  delete fJEC;
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
  fTightJets.clear();
  fGoodJets.clear();
  selectedJets8.clear();
  selectedJets15.clear();
  x1List.clear();
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
  
  pSNPt30        << "n" << iJetLabel << "sPt30";
  pSNPt30jesUp   << "n" << iJetLabel << "sPt30jesUp";
  pSNPt30jesDown << "n" << iJetLabel << "sPt30jesDown";
  pSNPt30jerUp   << "n" << iJetLabel << "sPt30jerUp";
  pSNPt30jerDown << "n" << iJetLabel << "sPt30jerDown";
  pSfwdPt30      << "n" << iJetLabel << "sfwdPt30";
  pSbPt30        << "n" << iJetLabel << "sbtagPt30";
  pSbLPt30       << "n" << iJetLabel << "sbtagLPt30";
  pSbMPt30       << "n" << iJetLabel << "sbtagMPt30";
  pSbTPt30       << "n" << iJetLabel << "sbtagTPt30";
  
  pSMetXCorrjesUp     << "MetXCorrjesUp";
  pSMetYCorrjesUp     << "MetYCorrjesUp";
  pSMetXCorrjesDown   << "MetXCorrjesDown";
  pSMetYCorrjesDown   << "MetYCorrjesDown";
  pSMetXCorrjerUp     << "MetXCorrjerUp";
  pSMetYCorrjerUp     << "MetYCorrjerUp";
  pSMetXCorrjerDown   << "MetXCorrjerDown";
  pSMetYCorrjerDown   << "MetYCorrjerDown";

  fTree->Branch(pSNPt30.str().c_str()           ,&fNJetsPt30           ,(pSNPt30.str()+"/I").c_str());  // jet multiplicity
  fTree->Branch(pSNPt30jesUp.str().c_str()      ,&fNJetsPt30jesUp      ,(pSNPt30jesUp.str()+"/I").c_str());  
  fTree->Branch(pSNPt30jesDown.str().c_str()    ,&fNJetsPt30jesDown    ,(pSNPt30jesDown.str()+"/I").c_str());
  fTree->Branch(pSNPt30jerUp.str().c_str()      ,&fNJetsPt30jerUp      ,(pSNPt30jerUp.str()+"/I").c_str());  
  fTree->Branch(pSNPt30jerDown.str().c_str()    ,&fNJetsPt30jerDown    ,(pSNPt30jerDown.str()+"/I").c_str());
  fTree->Branch(pSfwdPt30.str().c_str()         ,&fNFwdPt30            ,(pSfwdPt30.str()+"/I").c_str());
  fTree->Branch(pSbLPt30.str().c_str()          ,&fNBTagsLPt30         ,(pSbLPt30.str()+"/I").c_str()); // b tags
  fTree->Branch(pSbMPt30.str().c_str()          ,&fNBTagsMPt30         ,(pSbMPt30.str()+"/I").c_str());
  fTree->Branch(pSbTPt30.str().c_str()          ,&fNBTagsTPt30         ,(pSbTPt30.str()+"/I").c_str());
  
  fTree->Branch(pSMetXCorrjesUp.str().c_str()   ,&MetXCorrjesUp        ,(pSMetXCorrjesUp.str()+"/D").c_str());  // Met corrections
  fTree->Branch(pSMetYCorrjesUp.str().c_str()   ,&MetYCorrjesUp        ,(pSMetYCorrjesUp.str()+"/D").c_str());  
  fTree->Branch(pSMetXCorrjesDown.str().c_str() ,&MetXCorrjesDown      ,(pSMetXCorrjesDown.str()+"/D").c_str());  
  fTree->Branch(pSMetYCorrjesDown.str().c_str() ,&MetYCorrjesDown      ,(pSMetYCorrjesDown.str()+"/D").c_str());  
  fTree->Branch(pSMetXCorrjerUp.str().c_str()   ,&MetXCorrjerUp        ,(pSMetXCorrjerUp.str()+"/D").c_str());  
  fTree->Branch(pSMetYCorrjerUp.str().c_str()   ,&MetYCorrjerUp        ,(pSMetYCorrjerUp.str()+"/D").c_str());  
  fTree->Branch(pSMetXCorrjerDown.str().c_str() ,&MetXCorrjerDown      ,(pSMetXCorrjerDown.str()+"/D").c_str());  
  fTree->Branch(pSMetYCorrjerDown.str().c_str() ,&MetYCorrjerDown      ,(pSMetYCorrjerDown.str()+"/D").c_str());  

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
    pSbTPt100dR08<< "n" << iJetLabel << "sTPt100dR08_" << i0;
    pSbLPt150dR08<< "n" << iJetLabel << "sLPt150dR08_" << i0;
    pSbMPt150dR08<< "n" << iJetLabel << "sMPt150dR08_" << i0;
    pSbTPt150dR08<< "n" << iJetLabel << "sTPt150dR08_" << i0;
    fTree->Branch(pSNPt30dR08.str().c_str()        ,&fNJetsPt30dR08[i0]        ,(pSNPt30dR08.str()+"/I").c_str());
    fTree->Branch(pSNPt30dR08jesUp.str().c_str()   ,&fNJetsPt30dR08jesUp[i0]   ,(pSNPt30dR08jesUp.str()+"/I").c_str());
    fTree->Branch(pSNPt30dR08jesDown.str().c_str() ,&fNJetsPt30dR08jesDown[i0] ,(pSNPt30dR08jesDown.str()+"/I").c_str());
    fTree->Branch(pSNPt30dR08jerUp.str().c_str()   ,&fNJetsPt30dR08jerUp[i0]   ,(pSNPt30dR08jerUp.str()+"/I").c_str());
    fTree->Branch(pSNPt30dR08jerDown.str().c_str() ,&fNJetsPt30dR08jerDown[i0] ,(pSNPt30dR08jerDown.str()+"/I").c_str());
    fTree->Branch(pSbLPt50dR08.str().c_str()       ,&fNBTagsLPt50dR08[i0]      ,(pSbLPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt50dR08.str().c_str()       ,&fNBTagsMPt50dR08[i0]      ,(pSbMPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt50dR08.str().c_str()       ,&fNBTagsTPt50dR08[i0]      ,(pSbTPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbLPt100dR08.str().c_str()      ,&fNBTagsLPt100dR08[i0]     ,(pSbLPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt100dR08.str().c_str()      ,&fNBTagsMPt100dR08[i0]     ,(pSbMPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt100dR08.str().c_str()      ,&fNBTagsTPt100dR08[i0]     ,(pSbTPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbLPt150dR08.str().c_str()      ,&fNBTagsLPt150dR08[i0]     ,(pSbLPt150dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt150dR08.str().c_str()      ,&fNBTagsMPt150dR08[i0]     ,(pSbMPt150dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt150dR08.str().c_str()      ,&fNBTagsTPt150dR08[i0]     ,(pSbTPt150dR08.str()+"/I").c_str());
  }

  for(int i0 = 0; i0 < fN*(fNOtherVars)+fN*fNVars; i0++) {double pVar = 0; fVars.push_back(pVar);}           
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars); // fNVars*fN first vars are j*_pt,j*_eta,j*_phi for fN = j1,j2,j3,j4 (3*4=12)
  addOthers  (iJetLabel.c_str(),iTree,fN,fVars); // next fNOtherVars*fN vars (Mass, b-tag, qgid, dR, dPhi...) for j1,j2,j3,j4 (18*4=72)
}
void JetLoader::addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals) {
  for(int i0 = 0; i0 < iN; i0++) {
    int lBase = iN*fNVars+i0*fNOtherVars;
    std::stringstream pSMass,pSCSV,pSQGID,pSdR,pSdP,pScen,pSjesUp,pSjesDown,pSjerUp,pSjerDown;
    std::stringstream pSDeepCSVb,pSDeepCSVc,pSDeepCSVl,pSDeepCSVbb;
    std::stringstream pSDeepCMVAb,pSDeepCMVAc,pSDeepCMVAl,pSDeepCMVAbb;
    std::stringstream pSpartonF,pShadronF;
    pSMass      << iHeader << i0 << "_mass";
    pSCSV       << iHeader << i0 << "_csv";
    pSQGID      << iHeader << i0 << "_qgid";
    pSdR        << iHeader << i0 << "_dR08";
    pSdP        << iHeader << i0 << "_dPhi08";
    pScen       << iHeader << i0 << "_pt_old";
    pSjesUp     << iHeader << i0 << "_pt_JESUp";
    pSjesDown   << iHeader << i0 << "_pt_JESDown";
    pSjerUp     << iHeader << i0 << "_pt_JERUp";
    pSjerDown   << iHeader << i0 << "_pt_JERDown";
    pSDeepCSVb  << iHeader << i0 << "_deepcsvb";
    pSDeepCSVc  << iHeader << i0 << "_deepcsvc";
    pSDeepCSVl  << iHeader << i0 << "_deepcsvl";
    pSDeepCSVbb << iHeader << i0 << "_deepcsvbb";
    pSDeepCMVAb << iHeader << i0 << "_deepcmvab";
    pSDeepCMVAc << iHeader << i0 << "_deepcmvac";
    pSDeepCMVAl << iHeader << i0 << "_deepcmval";
    pSDeepCMVAbb<< iHeader << i0 << "_deepcmvabb";
    pSpartonF   << iHeader << i0 << "_partonFlavor";
    pShadronF   << iHeader << i0 << "_hadronFlavor";
    iTree->Branch(pSMass     .str().c_str() ,&iVals[lBase+0], (pSMass       .str()+"/D").c_str());
    iTree->Branch(pSCSV      .str().c_str() ,&iVals[lBase+1], (pSCSV        .str()+"/D").c_str());
    iTree->Branch(pSQGID     .str().c_str() ,&iVals[lBase+2], (pSQGID       .str()+"/D").c_str());
    iTree->Branch(pSdR       .str().c_str() ,&iVals[lBase+3], (pSdR         .str()+"/D").c_str());
    iTree->Branch(pSdP       .str().c_str() ,&iVals[lBase+4], (pSdP         .str()+"/D").c_str());
    iTree->Branch(pScen      .str().c_str() ,&iVals[lBase+5], (pScen        .str()+"/D").c_str());
    iTree->Branch(pSjesUp    .str().c_str() ,&iVals[lBase+6], (pSjesUp      .str()+"/D").c_str());
    iTree->Branch(pSjesDown  .str().c_str() ,&iVals[lBase+7], (pSjesDown    .str()+"/D").c_str());
    iTree->Branch(pSjerUp    .str().c_str() ,&iVals[lBase+8], (pSjerUp      .str()+"/D").c_str());
    iTree->Branch(pSjerDown  .str().c_str() ,&iVals[lBase+9], (pSjerDown    .str()+"/D").c_str());

    iTree->Branch(pSDeepCSVb  .str().c_str() ,&iVals[lBase+10],(pSDeepCSVb  .str()+"/D").c_str());
    iTree->Branch(pSDeepCSVc  .str().c_str() ,&iVals[lBase+11],(pSDeepCSVc  .str()+"/D").c_str());
    iTree->Branch(pSDeepCSVl  .str().c_str() ,&iVals[lBase+12],(pSDeepCSVl  .str()+"/D").c_str());
    iTree->Branch(pSDeepCSVbb .str().c_str() ,&iVals[lBase+13],(pSDeepCSVbb .str()+"/D").c_str());
    iTree->Branch(pSDeepCMVAb .str().c_str() ,&iVals[lBase+14],(pSDeepCMVAb .str()+"/D").c_str());
    iTree->Branch(pSDeepCMVAc .str().c_str() ,&iVals[lBase+15],(pSDeepCMVAc .str()+"/D").c_str());
    iTree->Branch(pSDeepCMVAl .str().c_str() ,&iVals[lBase+16],(pSDeepCMVAl .str()+"/D").c_str());
    iTree->Branch(pSDeepCMVAbb.str().c_str() ,&iVals[lBase+17],(pSDeepCMVAbb.str()+"/D").c_str());

    iTree->Branch(pSpartonF   .str().c_str() ,&iVals[lBase+18],(pSpartonF   .str()+"/I").c_str());
    iTree->Branch(pShadronF   .str().c_str() ,&iVals[lBase+19],(pShadronF   .str()+"/I").c_str());
  }
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
    double JEC = fJEC->JetEnergyCorrectionFactor(pJet->ptRaw, pJet->eta, pJet->phi, jetE, 
						 iRho, pJet->area, 
						 runNum);
    double jetCorrPt = JEC*(pJet->ptRaw);
    double jetCorrE = JEC*(vPJet.E());
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, pJet->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; 
    float sigma_MC = fJEC->resolution.getResolution(parameters);
    float sf = fJEC->resolution_sf.getScaleFactor(parameters);
    float sfUp = fJEC->resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = fJEC->resolution_sf.getScaleFactor(parameters, Variation::DOWN);
    double x1 = r->Gaus();
    double jetEnergySmearFactor = 1.0; 
    double jetEnergySmearFactorUp = 1.0; 
    double jetEnergySmearFactorDown = 1.0;    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
      if(sfUp < 1) jetEnergySmearFactorUp = jetEnergySmearFactor;
      else jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x1;
      if(sfDown < 1) jetEnergySmearFactorDown = jetEnergySmearFactor;
      else jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x1;
    }    
    double unc = fJEC->getJecUnc( jetCorrPt, pJet->eta, runNum ); //use run=999 as default
    
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

    if (jetCorrPtSmear > 30) lCountPt30++;
    if (jetPtJESUp     > 30) lCountPt30jesUp++;
    if (jetPtJESDown   > 30) lCountPt30jesDown++;
    if (jetPtJERUp     > 30) lCountPt30jerUp++;
    if (jetPtJERDown   > 30) lCountPt30jerDown++;
    
    // jet and b-tag multiplicity    
    vPJet.SetPtEtaPhiE(jetCorrPtSmear, pJet->eta, pJet->phi, jetCorrESmear);
    for (int i1 = 0; i1 < int(iVJets.size()); i1++) {      
      if(iVJets[i1].Pt()>200 && vPJet.DeltaR(iVJets[i1])>0.8) {
	if (jetCorrPtSmear  > 30) fNJetsPt30dR08[i1]++;
	if (jetPtJESUp      > 30) fNJetsPt30dR08jesUp[i1]++;
	if (jetPtJESDown    > 30) fNJetsPt30dR08jesDown[i1]++;
	if (jetPtJERUp      > 30) fNJetsPt30dR08jerUp[i1]++;
	if (jetPtJERDown    > 30) fNJetsPt30dR08jerDown[i1]++;
      }
    }
	
    if(jetCorrPtSmear  <=  30) continue;
    if(fabs(pJet->eta) > 2.5 && fabs(pJet->eta) < 4.5) lNFwdPt30++;
    //if(fabs(pJet->eta) >= 2.5) continue;
    if(!passJetTightSel(pJet,fYear)) continue;
    x1List.push_back(x1);
    addJet(pJet,fTightJets);
    fGoodJets.push_back(pJet);
    
    if(fabs(pJet->eta) < 2.5){
      if((pJet->deepcsvbb + pJet->deepcsvb) > DEEPCSVL) lNBTagLPt30++;
      if((pJet->deepcsvbb + pJet->deepcsvb) > DEEPCSVM) lNBTagMPt30++;
      if((pJet->deepcsvbb + pJet->deepcsvb) > DEEPCSVT) lNBTagTPt30++;
    }
    
    //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation
    // jet and b-tag multiplicity
    for (int i1 = 0; i1 < int(iVJets.size()); i1++) {
      if(iVJets[i1].Pt()>200 && vPJet.DeltaR(iVJets[i1])>0.8) {
	if(pJet->pt>50 && fabs(pJet->eta) < 2.5) {
	  if((pJet->deepcsvbb + pJet->deepcsvb) > DEEPCSVL) {
	    fNBTagsLPt50dR08[i1]++;
	    if(pJet->pt>100) fNBTagsLPt100dR08[i1]++;
            if(pJet->pt>150) fNBTagsLPt150dR08[i1]++;
	  }
          if((pJet->deepcsvbb + pJet->deepcsvb) > DEEPCSVM) {
            fNBTagsMPt50dR08[i1]++;
            if(pJet->pt>100) fNBTagsMPt100dR08[i1]++;
            if(pJet->pt>150) fNBTagsMPt150dR08[i1]++;
          }
          if((pJet->deepcsvbb + pJet->deepcsvb) > DEEPCSVT) {
            fNBTagsTPt50dR08[i1]++;
            if(pJet->pt>100) fNBTagsTPt100dR08[i1]++;
            if(pJet->pt>150) fNBTagsTPt150dR08[i1]++;
          }
	}
      }
    }

  } //end loop over jets
  addVJet(fTightJets,selectedJets);
  fNJetsPt30           = lCountPt30;
  fNJetsPt30jesUp      = lCountPt30jesUp;
  fNJetsPt30jesDown    = lCountPt30jesDown;
  fNJetsPt30jerUp      = lCountPt30jerUp;
  fNJetsPt30jerDown    = lCountPt30jerDown;
  fNFwdPt30            = lNFwdPt30;
  fNBTagsLPt30         = lNBTagLPt30;
  fNBTagsMPt30         = lNBTagMPt30;
  fNBTagsTPt30         = lNBTagTPt30;

  fillJetCorr(fN,fTightJets,fVars,iRho,runNum);
  fillOthers(fN,fTightJets,fVars,iVJets,iRho, runNum);
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
    double JEC = fJEC->JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
						 iRho, iObjects[i0]->area, 
						 runNum);
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    double x1 = x1List[i0];
    double jetEnergySmearFactor = 1.0;    
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; 
    float sigma_MC = fJEC->resolution.getResolution(parameters);
    float sf = fJEC->resolution_sf.getScaleFactor(parameters);    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    }    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    iVals[i0*3+0] = jetCorrPtSmear;
    iVals[i0*3+1] = iObjects[i0]->eta;
    iVals[i0*3+2] = iObjects[i0]->phi;
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

    double JEC = fJEC->JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
						 iRho, iObjects[i0]->area, 
						 runNum);
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    
    double unc = fJEC->getJecUnc( jetCorrPt, iObjects[i0]->eta, runNum ); //use run=999 as default
    
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}};
    float sigma_MC = fJEC->resolution.getResolution(parameters);
    float sf = fJEC->resolution_sf.getScaleFactor(parameters);
    float sfUp = fJEC->resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = fJEC->resolution_sf.getScaleFactor(parameters, Variation::DOWN);

    double x1 = x1List[i0];
    double jetEnergySmearFactor = 1.0; 
    double jetEnergySmearFactorUp = 1.0; 
    double jetEnergySmearFactorDown = 1.0;    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
      if(sfUp < 1) jetEnergySmearFactorUp = jetEnergySmearFactor;
      else jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x1;
      if(sfDown < 1) jetEnergySmearFactorDown = jetEnergySmearFactor;
      else jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x1;
    }    
    
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;
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

    // deepcsv
    iVals[lBase+i0*fNOtherVars+10] = iObjects[i0]->deepcsvb;
    iVals[lBase+i0*fNOtherVars+11] = iObjects[i0]->deepcsvc;
    iVals[lBase+i0*fNOtherVars+12] = iObjects[i0]->deepcsvl;
    iVals[lBase+i0*fNOtherVars+13] = iObjects[i0]->deepcsvbb;
    iVals[lBase+i0*fNOtherVars+14] = iObjects[i0]->deepcmvab;
    iVals[lBase+i0*fNOtherVars+15] = iObjects[i0]->deepcmvac;
    iVals[lBase+i0*fNOtherVars+16] = iObjects[i0]->deepcmval;
    iVals[lBase+i0*fNOtherVars+17] = iObjects[i0]->deepcmvabb;
    
    // flavour
    iVals[lBase+i0*fNOtherVars+18] = iObjects[i0]->partonFlavor;
    iVals[lBase+i0*fNOtherVars+19] = iObjects[i0]->hadronFlavor;

  }
}
