#include "../include/PerJetLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>
#include <unordered_set>

#define NARR 50

using namespace baconhep;

PerJetLoader::PerJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS,std::string iAddJetCHS,int iN, bool iData) { 
  fVJets         = new TClonesArray("baconhep::TJet");
  fVAddJets      = new TClonesArray("baconhep::TAddJet");
  fGens         = new TClonesArray("baconhep::TGenParticle");
//  fVJetsCHS      = new TClonesArray("baconhep::TJet");
//  fVAddJetsCHS   = new TClonesArray("baconhep::TAddJet");

  iTree->SetBranchAddress(iJet.c_str(),       &fVJets);
  iTree->SetBranchAddress(iAddJet.c_str(),    &fVAddJets);
  iTree->SetBranchAddress("GenParticle", &fGens);
//  iTree->SetBranchAddress(iJetCHS.c_str(),    &fVJetsCHS);
//  iTree->SetBranchAddress(iAddJetCHS.c_str(), &fVAddJetsCHS);

  fVJetBr        = iTree->GetBranch(iJet.c_str());
  fVAddJetBr     = iTree->GetBranch(iAddJet.c_str());
  fGenBr         = iTree->GetBranch("GenParticle");
//  fVJetBrCHS     = iTree->GetBranch(iJetCHS.c_str());
//  fVAddJetBrCHS  = iTree->GetBranch(iAddJetCHS.c_str());

  fN = iN;

  isData = iData;  
  loadJECs_Rereco(isData);

  r = new TRandom3(1988);

  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
}
PerJetLoader::~PerJetLoader() { 
  delete fVJets;
  delete fVJetBr;
  delete fVAddJets;
  delete fVAddJetBr;
  delete fGens;
  delete fGenBr;
 /* delete fVJetsCHS;
  delete fVJetBrCHS;
  delete fVAddJetsCHS;
  delete fVAddJetBrCHS;
  */
}
void PerJetLoader::reset() { 
  fNLooseVJets        = 0;
  fNTightVJets        = 0;  
  for(int i0 = 0; i0 < int(fisTightVJet.size()); i0++) fisTightVJet[i0] = -999;
  selectedVJets.clear();
  fLooseVJets.clear();
  x1List.clear();
  x2List.clear();
  x3List.clear();    
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
  resetZprime();  
}
void PerJetLoader::resetCHS() {
  fNLooseVJetsCHS     = 0;
  fNTightVJetsCHS     = 0;
  for(int i0 = 0; i0 < int(fisTightVJetCHS.size()); i0++) fisTightVJetCHS[i0] = -999;
  selectedVJetsCHS.clear();
  fLooseVJetsCHS.clear();
  fdoublecsvCHS.clear();
  fdoublesubCHS.clear();
  fptCHS.clear();
  fetaCHS.clear();
  fphiCHS.clear();
}
void PerJetLoader::resetDoubleB() {
  fLooseVJetsByDoubleB.clear();
  selectedVJetsByDoubleB.clear();
  //fLooseVJetsCHSByDoubleB.clear();
  //selectedVJetsCHSByDoubleB.clear();
}
void PerJetLoader::resetZprime() {
  fvSize              = -999;
  fvMatching          = -999;
  fisHadronicV        = 0;
  //fRatioPt =0;
}
void PerJetLoader::setupTree(TTree *iTree, std::string iJetLabel) { 
  reset();

  fSingletons.clear();
  fArrays.clear();

  fSingletons["mass"] = 0;
  fSingletons["csv"] = 0;
  fSingletons["CHF"] = 0;
  fSingletons["NHF"] = 0;
  fSingletons["NEMF"] = 0;
  fSingletons["tau21"] = 0;
  fSingletons["tau32"] = 0;
  fSingletons["msd"] = 0;
  fSingletons["rho"] = 0;
  fSingletons["minsubcsv"] = 0;
  fSingletons["maxsubcsv"] = 0;
  fSingletons["doublecsv"] = 0;
  fSingletons["doublesub"] = 0;
  fSingletons["ptraw"] = 0;
  fSingletons["genpt"] = 0;
  fSingletons["e2_b1"] = 0; // Correlation function inputs beta=1
  fSingletons["e3_b1"] = 0;
  fSingletons["e3_v1_b1"] = 0;
  fSingletons["e3_v2_b1"] = 0;
  fSingletons["e4_v1_b1"] = 0;
  fSingletons["e4_v2_b1"] = 0;
  fSingletons["e2_b2"] = 0; // Correlation function inputs beta=2
  fSingletons["e3_b2"] = 0;
  fSingletons["e3_v1_b2"] = 0;
  fSingletons["e3_v2_b2"] = 0;
  fSingletons["e4_v1_b2"] = 0;
  fSingletons["e4_v2_b2"] = 0;
  fSingletons["e2_sdb1"] = 0; // Correlation function inputs beta=1 soft-dropped 
  fSingletons["e3_sdb1"] = 0;
  fSingletons["e3_v1_sdb1"] = 0;
  fSingletons["e3_v2_sdb1"] = 0;
  fSingletons["e4_v1_sdb1"] = 0;
  fSingletons["e4_v2_sdb1"] = 0;
  fSingletons["e2_sdb2"] = 0; // Correlation function inputs beta=2 soft-dropped 
  fSingletons["e3_sdb2"] = 0;
  fSingletons["e3_v1_sdb2"] = 0;
  fSingletons["e3_v2_sdb2"] = 0;
  fSingletons["e4_v1_sdb2"] = 0;
  fSingletons["e4_v2_sdb2"] = 0;
  fSingletons["N2sdb1"] = 0; // 2-prong ECFs observables
  fSingletons["N2sdb2"] = 0;
  fSingletons["M2sdb1"] = 0;
  fSingletons["M2sdb2"] = 0;
  fSingletons["D2sdb1"] = 0;
  fSingletons["D2sdb2"] = 0;
  fSingletons["N2b1"] = 0;
  fSingletons["N2b2"] = 0;
  fSingletons["M2b1"] = 0;
  fSingletons["M2b2"] = 0;
  fSingletons["D2b1"] = 0;
  fSingletons["D2b2"] = 0;
  fSingletons["pt_old"] = 0;
  fSingletons["pt_JESUp"] = 0;
  fSingletons["pt_JESDown"] = 0;
  fSingletons["pt_JERUp"] = 0;
  fSingletons["pt_JERDown"] = 0;
  fSingletons["e2_sdb05"] = 0; // Correlation function inputs beta=0.5 soft-dropped 
  fSingletons["e3_sdb05"] = 0;
  fSingletons["e3_v1_sdb05"] = 0;
  fSingletons["e3_v2_sdb05"] = 0;
  fSingletons["e4_v1_sdb05"] = 0;
  fSingletons["e4_v2_sdb05"] = 0;
  fSingletons["e2_sdb4"] = 0; // Correlation function inputs beta=4 soft-dropped 
  fSingletons["e3_sdb4"] = 0;
  fSingletons["e3_v1_sdb4"] = 0;
  fSingletons["e3_v2_sdb4"] = 0;
  fSingletons["e4_v1_sdb4"] = 0;
  fSingletons["e4_v2_sdb4"] = 0;
  fSingletons["flavour"] = 0;
  fSingletons["nbHadrons"] = 0;
  fSingletons["nSV"] = 0;
  fSingletons["jetNTracks"] = 0;
  fSingletons["tau_flightDistance2dSig_1"] = 0;
  fSingletons["SubJet_csv"] = 0;
  fSingletons["z_ratio"] = 0;
  fSingletons["trackSipdSig_3"] = 0;
  fSingletons["trackSipdSig_2"] = 0;
  fSingletons["trackSipdSig_1"] = 0;
  fSingletons["trackSipdSig_0"] = 0;
  fSingletons["trackSipdSig_1_0"] = 0;
  fSingletons["trackSipdSig_0_0"] = 0;
  fSingletons["trackSipdSig_1_1"] = 0;
  fSingletons["trackSipdSig_0_1"] = 0;
  fSingletons["trackSip2dSigAboveCharm_0"] = 0;
  fSingletons["trackSip2dSigAboveBottom_0"] = 0;
  fSingletons["trackSip2dSigAboveBottom_1"] = 0;
  fSingletons["tau1_trackEtaRel_0"] = 0;
  fSingletons["tau1_trackEtaRel_1"] = 0;
  fSingletons["tau1_trackEtaRel_2"] = 0;
  fSingletons["tau0_trackEtaRel_0"] = 0;
  fSingletons["tau0_trackEtaRel_1"] = 0;
  fSingletons["tau0_trackEtaRel_2"] = 0;
  fSingletons["tau_vertexMass_0"] = 0;
  fSingletons["tau_vertexEnergyRatio_0"] = 0;
  fSingletons["tau_vertexDeltaR_0"] = 0;
  fSingletons["tau_flightDistance2dSig_0"] = 0;
  fSingletons["tau_vertexMass_1"] = 0;
  fSingletons["tau_vertexEnergyRatio_1"] = 0;
  fSingletons["nProngs"] = 0;

  fArrays["cpfPt"] = new float[NARR]; // example - add more

  fTree = iTree;

  for (auto &iter : fSingletons) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    fTree->Branch(bname.str().c_str(), &(iter.second), (bname+"/f").c_str());
  }
  for (auto &iter : fArrays) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    std::stringstream bname2;
    bname2 << bname << "[" << NARR << "]/f";
    fTree->Branch(bname.str().c_str(), &(iter.second), (bname+"/f").c_str());
  }

}

void PerJetLoader::setupTreeZprime(TTree *iTree, std::string iJetLabel) {
  resetZprime();
  std::stringstream pSiV;   pSiV << iJetLabel << "0_isHadronicV";
  std::stringstream pSVM;   pSVM << iJetLabel << "0_vMatching";
  std::stringstream pSVS;   pSVS << iJetLabel << "0_vSize";
  std::stringstream pSpF;   pSpF << iJetLabel << "0_partonFlavor";
  std::stringstream pShF;   pShF << iJetLabel << "0_hadronFlavor";
  std::stringstream pSnC;   pSnC << iJetLabel << "0_nCharged";
  std::stringstream pSnN;   pSnN << iJetLabel << "0_nNeutrals";
  std::stringstream pSnP;   pSnP << iJetLabel << "0_nParticles";
  std::stringstream pSvF;   pSvF << iJetLabel << "0_vertexFlavor"; 
  std::stringstream pSvFI;   pSvFI << iJetLabel << "0_vertexFlavorInfo";

  fTree = iTree;
  fTree->Branch(pSiV.str().c_str() ,&fisHadronicV         ,(pSiV.str()+"/I").c_str());
  fTree->Branch(pSVM.str().c_str() ,&fvMatching           ,(pSVM.str()+"/D").c_str());
  fTree->Branch(pSVS.str().c_str() ,&fvSize               ,(pSVS.str()+"/D").c_str());
  fTree->Branch(pSpF.str().c_str() ,&fpartonFlavor        ,(pSpF.str()+"/I").c_str());
  fTree->Branch(pShF.str().c_str() ,&fhadronFlavor        ,(pShF.str()+"/I").c_str());
  fTree->Branch(pSnC.str().c_str() ,&fnCharged            ,(pSnC.str()+"/I").c_str());
  fTree->Branch(pSnN.str().c_str() ,&fnNeutrals           ,(pSnN.str()+"/I").c_str());
  fTree->Branch(pSnP.str().c_str() ,&fnParticles          ,(pSnP.str()+"/I").c_str());
  fTree->Branch(pSvF.str().c_str() ,&fnVtxFlavor          ,(pSvF.str()+"/I").c_str());
  fTree->Branch(pSvFI.str().c_str() ,&fnVtxFlavInfo, (pSvFI.str()+"/I").c_str());
}

void PerJetLoader::load(int iEvent) { 
  fVJets       ->Clear();
  fVJetBr      ->GetEntry(iEvent);
  fVAddJets    ->Clear();
  fVAddJetBr   ->GetEntry(iEvent);
  fGens        ->Clear();
  fGenBr       ->GetEntry(iEvent);
}

void PerJetLoader::selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum){
  reset();  
  int lCount(0), lCountT(0);
  for  (int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    TJet *pVJet = (TJet*)((*fVJets)[i0]);    
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    
    double JEC_old = (pVJet->pt)/(pVJet->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(pVJet->ptRaw, pVJet->eta, pVJet->phi, (pVJet->mass)/JEC_old);
    double jetE = vPJet.E();
    double JEC = JetEnergyCorrectionFactor(pVJet->ptRaw, pVJet->eta, pVJet->phi, jetE, 
                 iRho, pVJet->area, 
                 runNum,
                JetCorrectionsIOV,JetCorrector);
    double jetCorrPt = JEC*(pVJet->ptRaw);
    double jetCorrE = JEC*(vPJet.E());
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, pVJet->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; // max 44.30 for Spring16_25nsV6_MC JER (CHANGE ONCE UPDATED)
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
    double unc = getJecUnc( jetCorrPt, pVJet->eta, runNum ); //use run=999 as default
    
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
    
    TLorentzVector thisJet;  thisJet.SetPtEtaPhiE(jetCorrPtSmear, pVJet->eta, pVJet->phi, jetCorrESmear);
    TLorentzVector thisJetJESUp;  thisJetJESUp.SetPtEtaPhiE(jetPtJESUp, pVJet->eta, pVJet->phi, jetEJESUp);
    TLorentzVector thisJetJESDown; thisJetJESDown.SetPtEtaPhiE(jetPtJESDown, pVJet->eta, pVJet->phi, jetEJESDown);
    TLorentzVector thisJetJERUp;  thisJetJERUp.SetPtEtaPhiE(jetPtJERUp,  pVJet->eta, pVJet->phi, jetEJERUp);
    TLorentzVector thisJetJERDown; thisJetJERDown.SetPtEtaPhiE(jetPtJERDown,  pVJet->eta, pVJet->phi, jetEJERDown);
    
    if(jetCorrPtSmear   <=  200)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJets);
    lCount++;
    x1List.push_back(x1);
    x2List.push_back(x2);
    x3List.push_back(x3);
    if(!passJetTightLepVetoSel(pVJet))                                     continue;
    lCountT++;
  }
  addVJet(fLooseVJets,selectedVJets);

  for  (int i0 = 0; i0 < int(selectedVJets.size()); i0++) { 
    if(passJetTightLepVetoSel(fLooseVJets[i0])) fisTightVJet[i0] = 1;
  }
  fNLooseVJets = lCount;
  fNTightVJets = lCountT;

  fillJetCorr( fN,fLooseVJets,fVars,iRho,runNum,fillEach);
  fillVJet(fN,fLooseVJets,fVars,dR,iRho,runNum,fillEach);
}

void PerJetLoader::fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum, bool fillEach){ 
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
    unsigned fillIdx = fillEach ? 0 : i0;
    iVals[fillIdx*3+0] = jetCorrPtSmear;
    iVals[fillIdx*3+1] = iObjects[i0]->eta;
    iVals[fillIdx*3+2] = iObjects[i0]->phi;
  }
}

void PerJetLoader::countVJetProngs(double dR) {
}

void PerJetLoader::selectVJetsByDoubleBCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum){
  // first do Puppi jets (pT > 500 GeV)
  resetDoubleB(); 
  for  (int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    TJet *pVJet = (TJet*)((*fVJets)[i0]);
    if(pVJet->pt        <=  500)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJetsByDoubleB);
    if(!passJetTightLepVetoSel(pVJet))                                     continue;
  }
  addVJet(fLooseVJetsByDoubleB,selectedVJetsByDoubleB);

  // now do CHS jets (pT > 400 GeV)
  for  (int i0 = 0; i0 < fVJetsCHS->GetEntriesFast(); i0++) {
    TJet *pVJet = (TJet*)((*fVJetsCHS)[i0]);
    if(pVJet->pt        <=  400)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJetsCHSByDoubleB);
    if(!passJetTightLepVetoSel(pVJet))                                     continue;
  }
  addVJet(fLooseVJetsCHSByDoubleB,selectedVJetsCHSByDoubleB);
  
  std::vector<int> indexCHS;
  std::vector<double> doubleBCHS;
  for (int i0 = 0; i0 < int(selectedVJetsByDoubleB.size()); i0++) {    
    int iCHSJet = -999;
    double dbCHSJet = -999;
    iCHSJet = getMatchedCHSJetIndex(selectedVJetsCHSByDoubleB,selectedVJetsByDoubleB[i0],0.8);
    if (iCHSJet > -999) {
      //selectedVJetsCHS[iCHSJet];      
      TAddJet *pAddJetCHS = getAddJetCHS(fLooseVJetsCHSByDoubleB[iCHSJet]);  
      dbCHSJet = pAddJetCHS->doublecsv;      
    }
    std::cout << "index CHS    = " << iCHSJet << std::endl;
    std::cout << "double-b CHS = " << dbCHSJet << std::endl;
    indexCHS.push_back(iCHSJet);
    doubleBCHS.push_back(dbCHSJet);
  }

  
}
void PerJetLoader::selectVJetsCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum){
  resetCHS();
  int lCount(0), lCountT(0);
  for  (int i0 = 0; i0 < fVJetsCHS->GetEntriesFast(); i0++) {
    TJet *pVJet = (TJet*)((*fVJetsCHS)[i0]);
    if(pVJet->pt        <=  150)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJetsCHS);
    lCount++;

    if(!passJetTightLepVetoSel(pVJet))                                     continue;
    lCountT++;
  }
  addVJet(fLooseVJetsCHS,selectedVJetsCHS);

  
  for  (int i0 = 0; i0 < int(selectedVJetsCHS.size()); i0++) { 
    if(passJetTightLepVetoSel(fLooseVJetsCHS[i0])) fisTightVJetCHS[i0] = 1;
  }
  fNLooseVJetsCHS = lCount;
  fNTightVJetsCHS = lCountT;
}
void PerJetLoader::fillVJet(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double dR, double iRho, unsigned int runNum, bool fillEach){ 
  int lBase = 3.*fN;
  int lMin = iObjects.size();
  int lNLabel = int(fLabels.size());
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) { 
    TAddJet *pAddJet = getAddJet(iObjects[i0]);
    //JEC    
    double x1 = x1List[i0];
    double x2 = x2List[i0];
    double x3 = x3List[i0];
    
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
    
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; // max 44.30 for Spring16_25nsV6_MC JER (CHANGE ONCE UPDATED)
    float sigma_MC = resolution.getResolution(parameters);
    float sf = resolution_sf.getScaleFactor(parameters);
    float sfUp = resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = resolution_sf.getScaleFactor(parameters, Variation::DOWN);

    
    double jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    double jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x2;
    double jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x3;
    
    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;
    /*
    if (true){
      std::cout << "VJet" << std::endl;
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
    unsigned fillIdx = fillEach ? 0 : i0;

    iVals[lBase+fillIdx*lNLabel+0]  = JEC*jetEnergySmearFactor*(iObjects[i0]->mass);
    iVals[lBase+fillIdx*lNLabel+1]  = iObjects[i0]->csv;
    iVals[lBase+fillIdx*lNLabel+2]  = iObjects[i0]->chHadFrac;
    iVals[lBase+fillIdx*lNLabel+3]  = iObjects[i0]->neuHadFrac;
    iVals[lBase+fillIdx*lNLabel+4]  = iObjects[i0]->neuEmFrac;
    iVals[lBase+fillIdx*lNLabel+5]  = (pAddJet->tau2/pAddJet->tau1);
    iVals[lBase+fillIdx*lNLabel+6]  = (pAddJet->tau3/pAddJet->tau2);
    iVals[lBase+fillIdx*lNLabel+7]  = pAddJet->mass_sd0;
    iVals[lBase+fillIdx*lNLabel+8]  = log((pAddJet->mass_sd0*pAddJet->mass_sd0)/iObjects[i0]->pt);
    iVals[lBase+fillIdx*lNLabel+9]  = TMath::Min(pAddJet->sj1_csv,pAddJet->sj2_csv);
    iVals[lBase+fillIdx*lNLabel+10] = TMath::Max(TMath::Max(pAddJet->sj1_csv,pAddJet->sj2_csv),TMath::Max(pAddJet->sj3_csv,pAddJet->sj4_csv));
    iVals[lBase+fillIdx*lNLabel+11] = pAddJet->doublecsv;
    iVals[lBase+fillIdx*lNLabel+12] = pAddJet->Double_sub;
    iVals[lBase+fillIdx*lNLabel+13] = iObjects[i0]->ptRaw;
    iVals[lBase+fillIdx*lNLabel+14] = iObjects[i0]->genpt;
    iVals[lBase+fillIdx*lNLabel+15] = pAddJet->e2_b1;
    iVals[lBase+fillIdx*lNLabel+16] = pAddJet->e3_b1;
    iVals[lBase+fillIdx*lNLabel+17] = pAddJet->e3_v1_b1;
    iVals[lBase+fillIdx*lNLabel+18] = pAddJet->e3_v2_b1;
    iVals[lBase+fillIdx*lNLabel+19] = pAddJet->e4_v1_b1;
    iVals[lBase+fillIdx*lNLabel+20] = pAddJet->e4_v2_b1;
    iVals[lBase+fillIdx*lNLabel+21] = pAddJet->e2_b2;
    iVals[lBase+fillIdx*lNLabel+22] = pAddJet->e3_b2;
    iVals[lBase+fillIdx*lNLabel+23] = pAddJet->e3_v1_b2;
    iVals[lBase+fillIdx*lNLabel+24] = pAddJet->e3_v2_b2;
    iVals[lBase+fillIdx*lNLabel+25] = pAddJet->e4_v1_b2;
    iVals[lBase+fillIdx*lNLabel+26] = pAddJet->e4_v2_b2;
    iVals[lBase+fillIdx*lNLabel+27] = pAddJet->e2_sdb1;
    iVals[lBase+fillIdx*lNLabel+28] = pAddJet->e3_sdb1;
    iVals[lBase+fillIdx*lNLabel+29] = pAddJet->e3_v1_sdb1;
    iVals[lBase+fillIdx*lNLabel+30] = pAddJet->e3_v2_sdb1;
    iVals[lBase+fillIdx*lNLabel+31] = pAddJet->e4_v1_sdb1;
    iVals[lBase+fillIdx*lNLabel+32] = pAddJet->e4_v2_sdb1;
    iVals[lBase+fillIdx*lNLabel+33] = pAddJet->e2_sdb2;
    iVals[lBase+fillIdx*lNLabel+34] = pAddJet->e3_sdb2;
    iVals[lBase+fillIdx*lNLabel+35] = pAddJet->e3_v1_sdb2;
    iVals[lBase+fillIdx*lNLabel+36] = pAddJet->e3_v2_sdb2;
    iVals[lBase+fillIdx*lNLabel+37] = pAddJet->e4_v1_sdb2;
    iVals[lBase+fillIdx*lNLabel+38] = pAddJet->e4_v2_sdb2;
    iVals[lBase+fillIdx*lNLabel+39] = pAddJet->e3_v2_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    iVals[lBase+fillIdx*lNLabel+40] = pAddJet->e3_v2_sdb2/(pAddJet->e2_sdb2*pAddJet->e2_sdb2);
    iVals[lBase+fillIdx*lNLabel+41] = pAddJet->e3_v1_sdb1/(pAddJet->e2_sdb1);
    iVals[lBase+fillIdx*lNLabel+42] = pAddJet->e3_v1_sdb2/(pAddJet->e2_sdb2);
    iVals[lBase+fillIdx*lNLabel+43] = pAddJet->e3_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    iVals[lBase+fillIdx*lNLabel+44] = pAddJet->e3_sdb2/(pAddJet->e2_sdb2*pAddJet->e2_sdb2*pAddJet->e2_sdb2);
    iVals[lBase+fillIdx*lNLabel+45] = pAddJet->e3_v2_b1/(pAddJet->e2_b1*pAddJet->e2_b1);
    iVals[lBase+fillIdx*lNLabel+46] = pAddJet->e3_v2_b2/(pAddJet->e2_b2*pAddJet->e2_b2);
    iVals[lBase+fillIdx*lNLabel+47] = pAddJet->e3_v1_b1/(pAddJet->e2_b1);
    iVals[lBase+fillIdx*lNLabel+48] = pAddJet->e3_v1_b2/(pAddJet->e2_b2);
    iVals[lBase+fillIdx*lNLabel+49] = pAddJet->e3_b1/(pAddJet->e2_b1*pAddJet->e2_b1*pAddJet->e2_b1);
    iVals[lBase+fillIdx*lNLabel+50] = pAddJet->e3_b2/(pAddJet->e2_b2*pAddJet->e2_b2*pAddJet->e2_b2);
    iVals[lBase+fillIdx*lNLabel+51] = iObjects[i0]->pt;
    iVals[lBase+fillIdx*lNLabel+52] = jetPtJESUp;
    iVals[lBase+fillIdx*lNLabel+53] = jetPtJESDown;
    iVals[lBase+fillIdx*lNLabel+54] = jetPtJERUp;
    iVals[lBase+fillIdx*lNLabel+55] = jetPtJERDown;
    iVals[lBase+fillIdx*lNLabel+56] = pAddJet->e2_sdb05;
    iVals[lBase+fillIdx*lNLabel+57] = pAddJet->e3_sdb05;
    iVals[lBase+fillIdx*lNLabel+58] = pAddJet->e3_v1_sdb05;
    iVals[lBase+fillIdx*lNLabel+59] = pAddJet->e3_v2_sdb05;
    iVals[lBase+fillIdx*lNLabel+60] = pAddJet->e4_v1_sdb05;
    iVals[lBase+fillIdx*lNLabel+61] = pAddJet->e4_v2_sdb05;
    iVals[lBase+fillIdx*lNLabel+62] = pAddJet->e2_sdb4;
    iVals[lBase+fillIdx*lNLabel+63] = pAddJet->e3_sdb4;
    iVals[lBase+fillIdx*lNLabel+64] = pAddJet->e3_v1_sdb4;
    iVals[lBase+fillIdx*lNLabel+65] = pAddJet->e3_v2_sdb4;
    iVals[lBase+fillIdx*lNLabel+66] = pAddJet->e4_v1_sdb4;
    iVals[lBase+fillIdx*lNLabel+67] = pAddJet->e4_v2_sdb4;
    iVals[lBase+fillIdx*lNLabel+68] = pAddJet->flavour;
    iVals[lBase+fillIdx*lNLabel+69] = pAddJet->nbHadrons;
    iVals[lBase+fillIdx*lNLabel+70] = pAddJet->nSV;
    iVals[lBase+fillIdx*lNLabel+71] = pAddJet->jetNTracks;
    iVals[lBase+fillIdx*lNLabel+72] = pAddJet->tau_flightDistance2dSig_1;
    iVals[lBase+fillIdx*lNLabel+73] = pAddJet->SubJet_csv;
    iVals[lBase+fillIdx*lNLabel+74] = pAddJet->z_ratio;
    iVals[lBase+fillIdx*lNLabel+75] = pAddJet->trackSipdSig_3;
    iVals[lBase+fillIdx*lNLabel+76] = pAddJet->trackSipdSig_2;
    iVals[lBase+fillIdx*lNLabel+77] = pAddJet->trackSipdSig_1;
    iVals[lBase+fillIdx*lNLabel+78] = pAddJet->trackSipdSig_0;
    iVals[lBase+fillIdx*lNLabel+79] = pAddJet->trackSipdSig_1_0;
    iVals[lBase+fillIdx*lNLabel+80] = pAddJet->trackSipdSig_0_0;
    iVals[lBase+fillIdx*lNLabel+81] = pAddJet->trackSipdSig_1_1;
    iVals[lBase+fillIdx*lNLabel+82] = pAddJet->trackSipdSig_0_1;
    iVals[lBase+fillIdx*lNLabel+83] = pAddJet->trackSip2dSigAboveCharm_0;
    iVals[lBase+fillIdx*lNLabel+84] = pAddJet->trackSip2dSigAboveBottom_0;
    iVals[lBase+fillIdx*lNLabel+85] = pAddJet->trackSip2dSigAboveBottom_1;
    iVals[lBase+fillIdx*lNLabel+86] = pAddJet->tau1_trackEtaRel_0;
    iVals[lBase+fillIdx*lNLabel+87] = pAddJet->tau1_trackEtaRel_1;
    iVals[lBase+fillIdx*lNLabel+88] = pAddJet->tau1_trackEtaRel_2;
    iVals[lBase+fillIdx*lNLabel+89] = pAddJet->tau0_trackEtaRel_0;
    iVals[lBase+fillIdx*lNLabel+90] = pAddJet->tau0_trackEtaRel_1;
    iVals[lBase+fillIdx*lNLabel+91] = pAddJet->tau0_trackEtaRel_2;
    iVals[lBase+fillIdx*lNLabel+92] = pAddJet->tau_vertexMass_0;
    iVals[lBase+fillIdx*lNLabel+93] = pAddJet->tau_vertexEnergyRatio_0;
    iVals[lBase+fillIdx*lNLabel+94] = pAddJet->tau_vertexDeltaR_0;
    iVals[lBase+fillIdx*lNLabel+95] = pAddJet->tau_flightDistance2dSig_0;
    iVals[lBase+fillIdx*lNLabel+96] = pAddJet->tau_vertexMass_1;
    iVals[lBase+fillIdx*lNLabel+97] = pAddJet->tau_vertexEnergyRatio_1;

    unsigned nG = fGens->GetEntriesFast();
    unsigned nP = 0;
    double dR2 = dR * dR;
    double threshold = 0.2 * iObjects[i0]->pt;
    std::unordered_set<TGenParticle*> partons; // avoid double-counting
    for (unsigned iG = 0; iG != nG; ++iG) {
      TGenParticle *part = (TGenParticle*)((*fGens)[iG]);
      unsigned apdgid = abs(part->pdgId);
      if (apdgid > 5 &&
          apdgid != 21 &&
          apdgid != 15 &&
          apdgid != 11 &&
          apdgid != 13) 
        continue;

      TGenParticle *parent = part;
      bool foundParent = false;
      while (parent->parent > 0) {
        parent = (TGenParticle*)((*fGens)[parent->parent]);
        if (partons.find(parent) != partons.end()) {
          foundParent = true;
          break;
        }
      }
      if (foundParent)
        continue;

      if (part->pt < threshold)
        continue;
      if (deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, part->eta, part->phi) > dR2)
        continue;
      partons.insert(part);
      ++nP;
    }  
    iVals[lBase+fillIdx*lNLabel+98] = nP;


    fpartonFlavor   = iObjects[0]->partonFlavor;
    fhadronFlavor   = iObjects[0]->hadronFlavor;
    fnCharged       = iObjects[0]->nCharged;
    fnNeutrals      = iObjects[0]->nNeutrals;
    fnParticles     = iObjects[0]->nParticles;
    fnVtxFlavor     = iObjects[0]->vtxFlavor;
    fnVtxFlavInfo   = iObjects[0]->vtxFlavInfo;

    if (fillEach)
      fTree->Fill(); // do the fill here to get all arrays

  }
}
int PerJetLoader::getMatchedCHSJetIndex(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR) {  
  TLorentzVector iJet1;
  int iJet1id(0), nmatched(0);
  float mindR = dR;
  for(int i0 = 0; i0 < int(iJets1.size()); i0++) {
    if ((iJets1[i0].DeltaR(iJet2) < mindR) && (fabs(iJets1[i0].Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))) {
      nmatched++;
      iJet1 = iJets1[i0];
      iJet1id = i0;
      mindR= iJets1[i0].DeltaR(iJet2);  
    }
  }
  if (nmatched >0 && (iJet1.DeltaR(iJet2) < dR) && (fabs(iJet1.Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))){    
    return iJet1id;
  }
  return -999;
}
void PerJetLoader::matchJet(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR, int jIndex){
  TLorentzVector iJet1;
  int iJet1id(0), nmatched(0);
  float mindR = dR;
  for(int i0 = 0; i0 < int(iJets1.size()); i0++) {
    if ((iJets1[i0].DeltaR(iJet2) < mindR) && (fabs(iJets1[i0].Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))) {
      nmatched++;
      iJet1 = iJets1[i0];
      iJet1id = i0;
      mindR= iJets1[i0].DeltaR(iJet2);  
    }
  }
  if (nmatched >0 && (iJet1.DeltaR(iJet2) < dR) && (fabs(iJet1.Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))){
    fillVJetCHS(fLooseVJetsCHS[iJet1id], jIndex);
  }
}

/*void PerJetLoader::matchJet15(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR){
  TLorentzVector iJet1;
  int nmatched(0);
  float mindR = dR;  
  for(int i0 = 0; i0 < int(iJets1.size()); i0++) {
    if ((iJets1[i0].DeltaR(iJet2) < mindR)) {
      nmatched++;
      iJet1 = iJets1[i0];
      mindR= iJets1[i0].DeltaR(iJet2);
    }
  }
  if (nmatched >0 && (iJet1.DeltaR(iJet2) < dR)){
    fRatioPt = iJet2.Pt()/iJet1.Pt();
  }
}
*/
void PerJetLoader::fillVJetCHS(TJet *iJet, int jIndex){
  TAddJet *pAddJet = getAddJetCHS(iJet);
  fdoublecsvCHS[jIndex] = double(pAddJet->doublecsv);
  fdoublesubCHS[jIndex] = double(pAddJet->Double_sub);
  fptCHS[jIndex] = double(iJet->pt);
  fetaCHS[jIndex] = double(iJet->eta);
  fphiCHS[jIndex] = double(iJet->phi);
}
TAddJet *PerJetLoader::getAddJet(TJet *iJet) { 
  int lIndex = -1;
  TAddJet *lJet = 0; 
  for(int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    if((*fVJets)[i0] == iJet) { lIndex = i0; break;}
  }
  if(lIndex == -1) return 0;
  for  (int i0 = 0; i0 < fVAddJets->GetEntriesFast(); i0++) { 
    TAddJet *pJet = (TAddJet*)((*fVAddJets)[i0]);
    if(pJet->index == fabs(lIndex)) { lJet = pJet; break;}
  }
  return lJet;
}
TAddJet *PerJetLoader::getAddJetCHS(TJet *iJet) {
  int lIndex = -1;
  TAddJet *lJet = 0;
  for(int i0 = 0; i0 < fVJetsCHS->GetEntriesFast(); i0++) {
    if((*fVJetsCHS)[i0] == iJet) { lIndex = i0; break;}
  }
  if(lIndex == -1) return 0;
  for  (int i0 = 0; i0 < fVAddJetsCHS->GetEntriesFast(); i0++) {
    TAddJet *pJet = (TAddJet*)((*fVAddJetsCHS)[i0]);
    if(pJet->index == fabs(lIndex)) { lJet = pJet; break;}
  }
  return lJet;
}

//2016 Prompt Reco
void PerJetLoader::loadJECs(bool isData) {
    std::cout << "PerJetLoader: loading jet energy correction constants" << std::endl;
    // initialize
    loadCMSSWPath();
    std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    resolution = JME::JetResolution(Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_PtResolution_AK8PFPuppi.txt",jecPathname.c_str()));
    resolution_sf = JME::JetResolutionScaleFactor(Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_SF_AK8PFPuppi.txt",jecPathname.c_str()));
    
    if (isData) {      
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));

      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_MC/Spring16_25nsV6_MC_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }

}
void PerJetLoader::loadJECs_Rereco(bool isData) {
    std::cout << "PerJetLoader: loading Rereco jet energy correction constants" << std::endl;
    // initialize
    loadCMSSWPath();
    std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    resolution = JME::JetResolution(Form("%s/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt",jecPathname.c_str()));
    resolution_sf = JME::JetResolutionScaleFactor(Form("%s/Spring16_25nsV10_MC/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",jecPathname.c_str()));
 
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      std::string jecUncPathEF = jecPathname+"/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(jecUncPathEF);

      correctionParameters.push_back(correctionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278801 ));

      //IOV: 2016G
      std::vector<JetCorrectorParameters> correctionParametersG = std::vector<JetCorrectorParameters> ();
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorG = new FactorizedJetCorrector(correctionParametersG);
      std::string jecUncPathG = jecPathname+"/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncG = new JetCorrectionUncertainty(jecUncPathG);

      correctionParameters.push_back(correctionParametersG);
      JetCorrector.push_back( JetCorrectorG );
      jecUnc.push_back(jecUncG);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 280385 ));

      //IOV: 2016H
      std::vector<JetCorrectorParameters> correctionParametersH = std::vector<JetCorrectorParameters> ();
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorH = new FactorizedJetCorrector(correctionParametersH);
      std::string jecUncPathH = jecPathname+"/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncH = new JetCorrectionUncertainty(jecUncPathH);

      correctionParameters.push_back(correctionParametersH);
      JetCorrector.push_back( JetCorrectorH );
      jecUnc.push_back(jecUncH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 280919, 99999999 ));

    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }
  
}

void PerJetLoader::loadCMSSWPath() {
    char* cmsswPathChar = getenv("CMSSW_BASE");
    if (cmsswPathChar == NULL) {
        std::cout << "Warning in PerJetLoader::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
        cmsswPath = "";
    }
    cmsswPath = std::string(cmsswPathChar);
}



// Retrieve jet energy uncertainty as a function of pt and eta
double PerJetLoader::getJecUnc( float pt, float eta , int run) {

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
double PerJetLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
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
double PerJetLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
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

