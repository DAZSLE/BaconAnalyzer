#include "../include/VJetLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>

using namespace baconhep;

// sub-jet ordering
class sjpair {
public:
  sjpair (float d, float p, float j){
    dR = (d>-1) ? d : 999;
    dPhi = (p>-1) ? p : -999;
    dPhiJRF = (j>-1) ? j : -999;
  }
  ~sjpair() { }
  float dR=999;
  float dPhi=-999;
  float dPhiJRF=-999;
};
bool compsjpairs(sjpair p1, sjpair p2){
  return p1.dR<p2.dR;
}
bool compsjpairsdPhi(sjpair p1, sjpair p2){
  return p1.dPhi>p2.dPhi;
}
bool compsjpairsdPhiJRF(sjpair p1, sjpair p2){
  return p1.dPhiJRF>p2.dPhiJRF;
}
double clean(double x, double def=-1) {
  if (!(x==x)) return def;
  else return x;
}

VJetLoader::VJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS,std::string iAddJetCHS,int iN, std::string subjetbtagScaleFactorFilename) { 
  fVJets         = new TClonesArray("baconhep::TJet");
  fVAddJets      = new TClonesArray("baconhep::TAddJet");
  fVJetsCHS      = new TClonesArray("baconhep::TJet");
  fVAddJetsCHS   = new TClonesArray("baconhep::TAddJet");

  iTree->SetBranchAddress(iJet.c_str(),       &fVJets);
  iTree->SetBranchAddress(iAddJet.c_str(),    &fVAddJets);
  iTree->SetBranchAddress(iJetCHS.c_str(),    &fVJetsCHS);
  iTree->SetBranchAddress(iAddJetCHS.c_str(), &fVAddJetsCHS);

  fVJetBr        = iTree->GetBranch(iJet.c_str());
  fVAddJetBr     = iTree->GetBranch(iAddJet.c_str());
  fVJetBrCHS     = iTree->GetBranch(iJetCHS.c_str());
  fVAddJetBrCHS  = iTree->GetBranch(iAddJetCHS.c_str());

  fN = iN;

  fSubJetCalib = new BTagCalibration("csvv2",subjetbtagScaleFactorFilename);
  fSubJetreadersL.clear(); fSubJetreadersM.clear();
  fSubJetreaders.clear();
  for(auto imtype : measurementTypes) { // fSubJetreadersL 6, M6 
    for(auto ivtype : variationTypes) {
      fSubJetreadersL.push_back(new BTagCalibrationReader(fSubJetCalib, BTagEntry::OP_LOOSE,  imtype, ivtype)); // first lt(HF) then incl(LF) and first central(0,3) then up(1,4) and then down(2,5)
      fSubJetreadersM.push_back(new BTagCalibrationReader(fSubJetCalib, BTagEntry::OP_MEDIUM, imtype, ivtype));
    }
  }
  fSubJetreaders.push_back(fSubJetreadersL); fSubJetreaders.push_back(fSubJetreadersM);
}
VJetLoader::~VJetLoader() { 
  delete fVJets;
  delete fVJetBr;
  delete fVAddJets;
  delete fVAddJetBr;
  delete fVJetsCHS;
  delete fVJetBrCHS;
  delete fVAddJetsCHS;
  delete fVAddJetBrCHS;
}
void VJetLoader::resetSubJetBTag() {
  for(unsigned int i0 = 0; i0 < fSubJetBTagVars.size(); i0++) fSubJetBTagVars[i0] = 1;
}
void VJetLoader::reset() { 
  fNLooseVJets        = 0;
  fNTightVJets        = 0;  
  for(int i0 = 0; i0 < int(fisTightVJet.size()); i0++) fisTightVJet[i0] = -999;
  selectedVJets.clear();
  fLooseVJets.clear();
  fGoodVSubJets.clear();
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
  resetSubJetBTag();
  resetZprime();
  resetMonoX();
}
void VJetLoader::resetCHS() {
  fNLooseVJetsCHS     = 0;
  fNTightVJetsCHS     = 0;
  for(int i0 = 0; i0 < int(fisTightVJetCHS.size()); i0++) fisTightVJetCHS[i0] = -999;
  selectedVJetsCHS.clear();
  fLooseVJetsCHS.clear();
  fGoodVSubJetsCHS.clear();
  fdoublecsvCHS.clear();
  fdoublesubCHS.clear();
  fptCHS.clear();
  fetaCHS.clear();
  fphiCHS.clear();
}
void VJetLoader::resetDoubleB() {
  fLooseVJetsByDoubleB.clear();
  selectedVJetsByDoubleB.clear();
  fLooseVJetsCHSByDoubleB.clear();
  selectedVJetsCHSByDoubleB.clear();
}
void VJetLoader::resetZprime() {
  fvSize              = 999;
  fvMatching          = 999;
  fisHadronicV        = 0;
  fRatioPt =0;
}
void VJetLoader::resetMonoX() {
  fVMT                = 0;
  fdR_sj0dR           = 999;
  fdPhi_sj0dPhi       = -999;
  fdPhiJRF_sj0dPhiJRF = -999;
  ftopSize            = 999;
  ftopMatching        = 999;
  fisHadronicTop      = 0;
}
void VJetLoader::setupTree(TTree *iTree, std::string iJetLabel) { 
  reset();
  fLabels.clear();
  fLabels.push_back("mass");
  fLabels.push_back("csv");
  fLabels.push_back("CHF");
  fLabels.push_back("NHF");
  fLabels.push_back("NEMF");
  fLabels.push_back("tau21");
  fLabels.push_back("tau32");
  fLabels.push_back("msd");
  fLabels.push_back("rho");
  fLabels.push_back("minsubcsv");
  fLabels.push_back("maxsubcsv");
  fLabels.push_back("doublecsv");
  fLabels.push_back("doublesub");
  fLabels.push_back("ptraw");
  fLabels.push_back("genpt");
  fLabels.push_back("e2_b1"); // Correlation function inputs beta=1
  fLabels.push_back("e3_b1");
  fLabels.push_back("e3_v1_b1");
  fLabels.push_back("e3_v2_b1");
  fLabels.push_back("e4_v1_b1");
  fLabels.push_back("e4_v2_b1");
  fLabels.push_back("e2_b2"); // Correlation function inputs beta=2
  fLabels.push_back("e3_b2");
  fLabels.push_back("e3_v1_b2");
  fLabels.push_back("e3_v2_b2");
  fLabels.push_back("e4_v1_b2");
  fLabels.push_back("e4_v2_b2");
  fLabels.push_back("e2_sdb1"); // Correlation function inputs beta=1 soft-dropped 
  fLabels.push_back("e3_sdb1");
  fLabels.push_back("e3_v1_sdb1");
  fLabels.push_back("e3_v2_sdb1");
  fLabels.push_back("e4_v1_sdb1");
  fLabels.push_back("e4_v2_sdb1");
  fLabels.push_back("e2_sdb2"); // Correlation function inputs beta=2 soft-dropped 
  fLabels.push_back("e3_sdb2");
  fLabels.push_back("e3_v1_sdb2");
  fLabels.push_back("e3_v2_sdb2");
  fLabels.push_back("e4_v1_sdb2");
  fLabels.push_back("e4_v2_sdb2");
  fLabels.push_back("N2sdb1"); // 2-prong ECFs observables
  fLabels.push_back("N2sdb2");
  fLabels.push_back("M2sdb1");
  fLabels.push_back("M2sdb2");
  fLabels.push_back("D2sdb1");
  fLabels.push_back("D2sdb2");
  fLabels.push_back("N2b1");
  fLabels.push_back("N2b2");
  fLabels.push_back("M2b1");
  fLabels.push_back("M2b2");
  fLabels.push_back("D2b1");
  fLabels.push_back("D2b2");

  std::stringstream pSNJ;   pSNJ << "n" << iJetLabel << "s";
  fTree = iTree;
  for(int i0 = 0; i0 < fN*4.;                    i0++) {double pVar = 0; fVars.push_back(pVar);} // declare array of vars
  for(int i0 = 0; i0 < fN*(int(fLabels.size())); i0++) {double pVar = 0; fVars.push_back(pVar);} 
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars);                                                 // from MonoXUtils.cc => fN =1 *_pt,*_eta,*_phi for vjet0 (3*1=3)
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars,fN*3,fLabels);
  fTree->Branch(pSNJ.str().c_str() ,&fNLooseVJets         ,(pSNJ.str()+"/I").c_str());  
  for(int i0 = 0; i0 < fN;                    i0++) fisTightVJet.push_back(-999);
  for(int i0 = 0; i0 < fN;                    i0++) {
    std::stringstream pSTJ;   pSTJ << iJetLabel << i0 << "_isTightVJet";
    fTree->Branch(pSTJ.str().c_str() ,&fisTightVJet[i0]         ,(pSTJ.str()+"/I").c_str());
  }
}
void VJetLoader::setupTreeZprime(TTree *iTree, std::string iJetLabel) {
  resetZprime();
  std::stringstream pSiV;   pSiV << iJetLabel << "0_isHadronicV";
  std::stringstream pSVM;   pSVM << iJetLabel << "0_vMatching";
  std::stringstream pSVS;   pSVS << iJetLabel << "0_vSize";
  std::stringstream pSpF;   pSpF << iJetLabel << "0_partonFlavor";
  std::stringstream pShF;   pShF << iJetLabel << "0_hadronFlavor";
  std::stringstream pSnC;   pSnC << iJetLabel << "0_nCharged";
  std::stringstream pSnN;   pSnN << iJetLabel << "0_nNeutrals";
  std::stringstream pSnP;   pSnP << iJetLabel << "0_nParticles";
  std::stringstream pSratio; pSratio << iJetLabel << "0_ratioCA15_04";

  fTree = iTree;
  fTree->Branch(pSiV.str().c_str() ,&fisHadronicV         ,(pSiV.str()+"/I").c_str());
  fTree->Branch(pSVM.str().c_str() ,&fvMatching           ,(pSVM.str()+"/D").c_str());
  fTree->Branch(pSVS.str().c_str() ,&fvSize               ,(pSVS.str()+"/D").c_str());
  fTree->Branch(pSpF.str().c_str() ,&fpartonFlavor        ,(pSpF.str()+"/I").c_str());
  fTree->Branch(pShF.str().c_str() ,&fhadronFlavor        ,(pShF.str()+"/I").c_str());
  fTree->Branch(pSnC.str().c_str() ,&fnCharged            ,(pSnC.str()+"/I").c_str());
  fTree->Branch(pSnN.str().c_str() ,&fnNeutrals           ,(pSnN.str()+"/I").c_str());
  fTree->Branch(pSnP.str().c_str() ,&fnParticles          ,(pSnP.str()+"/I").c_str());
  fTree->Branch(pSratio.str().c_str() ,&fRatioPt          ,(pSratio.str()+"/D").c_str());
}
void VJetLoader::setupTreeCHS(TTree *iTree, std::string iJetLabel) {
  resetCHS();  
  fTree = iTree;
  for(int i0 = 0; i0 < fN; i0++) {
    fdoublecsvCHS.push_back(-999);
    fdoublesubCHS.push_back(-999);
    fptCHS.push_back(-999);
    fetaCHS.push_back(-999);
    fphiCHS.push_back(-999);
    fisTightVJetCHS.push_back(-999);
  }
  for(int i0 = 0; i0 < fN; i0++) {
    std::stringstream pSdc;   pSdc << iJetLabel << i0 << "_doublecsv";
    std::stringstream pSds;   pSds << iJetLabel << i0 << "_doublesub";    
    std::stringstream pSpt;   pSpt << iJetLabel << i0 << "_pt";    
    std::stringstream pSeta;   pSeta << iJetLabel << i0 << "_eta";  
    std::stringstream pSphi;   pSphi << iJetLabel << i0 << "_phi";
    std::stringstream pSTJ;   pSTJ << iJetLabel << i0 << "_isTightVJet";    
    fTree->Branch(pSTJ.str().c_str() ,&fisTightVJetCHS[i0]         ,(pSTJ.str()+"/I").c_str());
    fTree->Branch(pSdc.str().c_str() ,&fdoublecsvCHS.at(i0)        ,(pSdc.str()+"/D").c_str());
    fTree->Branch(pSds.str().c_str() ,&fdoublesubCHS.at(i0)       ,(pSds.str()+"/D").c_str());
    fTree->Branch(pSpt.str().c_str() ,&fptCHS.at(i0)       ,(pSpt.str()+"/D").c_str());
    fTree->Branch(pSeta.str().c_str() ,&fetaCHS.at(i0)       ,(pSeta.str()+"/D").c_str());
    fTree->Branch(pSphi.str().c_str() ,&fphiCHS.at(i0)       ,(pSphi.str()+"/D").c_str());
  }
}
void VJetLoader::setupTreeMonoX(TTree *iTree, std::string iJetLabel) {
  std::stringstream pSMT;   pSMT << iJetLabel << "0_mT";
  std::stringstream pSdR;   pSdR << iJetLabel << "0_dR_sj0dR";
  std::stringstream pSdJ;   pSdJ << iJetLabel << "0_dPhiJRF_sj0dPhiJRF";
  std::stringstream pSiT;   pSiT << iJetLabel << "0_isHadronicTop";
  std::stringstream pSTM;   pSTM << iJetLabel << "0_topMatching";
  std::stringstream pSTS;   pSTS << iJetLabel << "0_topSize";
  fTree = iTree;
  fTree->Branch(pSMT.str().c_str() ,&fVMT                 ,(pSMT.str()+"/F").c_str());
  fTree->Branch(pSdR.str().c_str() ,&fdR_sj0dR            ,(pSdR.str()+"/F").c_str());
  fTree->Branch(pSdJ.str().c_str() ,&fdPhiJRF_sj0dPhiJRF  ,(pSdJ.str()+"/F").c_str());
  fTree->Branch(pSiT.str().c_str() ,&fisHadronicTop       ,(pSiT.str()+"/I").c_str());
  fTree->Branch(pSTM.str().c_str() ,&ftopMatching         ,(pSTM.str()+"/D").c_str());
  fTree->Branch(pSTS.str().c_str() ,&ftopSize             ,(pSTS.str()+"/D").c_str());
}
void VJetLoader::setupTreeSubJetBTag(TTree *iTree, std::string iJetLabel) {
  resetSubJetBTag();
  fTree = iTree;
  for(int i0 = 0; i0 < 40; i0++) {float pBTagVar = 1; fSubJetBTagVars.push_back(pBTagVar);} // declare array of 40 vars ( L0,L1,Lminus1,L2, M0,M1,Mminus1,M2) for (CENT,MISTAGUP,MISTAGDO,BTAGUP,BTAGDO)
  int i1 = 0;
  for(auto iwptype : wpTypes) {
    addSubJetBTag(iJetLabel.c_str(),iTree,iwptype,fBtagLabels,i1,fSubJetBTagVars);
    i1 += 20;
  }
}
void VJetLoader::load(int iEvent) { 
  fVJets       ->Clear();
  fVJetBr      ->GetEntry(iEvent);
  fVAddJets    ->Clear();
  fVAddJetBr   ->GetEntry(iEvent);
  fVJetsCHS    ->Clear();
  fVJetBrCHS   ->GetEntry(iEvent);
  fVAddJetsCHS ->Clear();
  fVAddJetBrCHS->GetEntry(iEvent);
}
void VJetLoader::selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho){
  reset(); 
  int lCount(0), lCountT(0);
  for  (int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    TJet *pVJet = (TJet*)((*fVJets)[i0]);
    if(pVJet->pt        <=  350)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJets);
    lCount++;

    if(!passJetTightLepVetoSel(pVJet))                                     continue;
    lCountT++;
  }
  addVJet(fLooseVJets,selectedVJets);

  for  (int i0 = 0; i0 < int(selectedVJets.size()); i0++) { 
    if(passJetTightLepVetoSel(fLooseVJets[i0])) fisTightVJet[i0] = 1;
  }
  fNLooseVJets = lCount;
  fNTightVJets = lCountT;

  fillJet( fN,fLooseVJets,fVars);
  fillVJet(fN,fLooseVJets,fVars,iRho);
}

void VJetLoader::selectVJetsByDoubleBCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho){
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
void VJetLoader::selectVJetsCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho){
  resetCHS();
  int lCount(0), lCountT(0);
  for  (int i0 = 0; i0 < fVJetsCHS->GetEntriesFast(); i0++) {
    TJet *pVJet = (TJet*)((*fVJetsCHS)[i0]);
    if(pVJet->pt        <=  300)                                           continue;
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
void VJetLoader::fillVJet(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho){ 
  int lBase = 3.*fN;
  int lMin = iObjects.size();
  int lNLabel = int(fLabels.size());
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) { 
    TAddJet *pAddJet = getAddJet(iObjects[i0]);
    iVals[lBase+i0*lNLabel+0]  = iObjects[i0]->mass;
    iVals[lBase+i0*lNLabel+1]  = iObjects[i0]->csv;
    iVals[lBase+i0*lNLabel+2]  = iObjects[i0]->chHadFrac;
    iVals[lBase+i0*lNLabel+3]  = iObjects[i0]->neuHadFrac;
    iVals[lBase+i0*lNLabel+4]  = iObjects[i0]->neuEmFrac;
    iVals[lBase+i0*lNLabel+5]  = (pAddJet->tau2/pAddJet->tau1);
    iVals[lBase+i0*lNLabel+6]  = (pAddJet->tau3/pAddJet->tau2);
    iVals[lBase+i0*lNLabel+7]  = pAddJet->mass_sd0;
    iVals[lBase+i0*lNLabel+8]  = log((pAddJet->mass_sd0*pAddJet->mass_sd0)/iObjects[i0]->pt);
    iVals[lBase+i0*lNLabel+9]  = TMath::Min(pAddJet->sj1_csv,pAddJet->sj2_csv);
    iVals[lBase+i0*lNLabel+10] = TMath::Max(TMath::Max(pAddJet->sj1_csv,pAddJet->sj2_csv),TMath::Max(pAddJet->sj3_csv,pAddJet->sj4_csv));
    iVals[lBase+i0*lNLabel+11] = pAddJet->doublecsv;
    iVals[lBase+i0*lNLabel+12] = pAddJet->Double_sub;
    iVals[lBase+i0*lNLabel+13] = iObjects[i0]->ptRaw;
    iVals[lBase+i0*lNLabel+14] = iObjects[i0]->genpt;
    iVals[lBase+i0*lNLabel+15] = pAddJet->e2_b1;
    iVals[lBase+i0*lNLabel+16] = pAddJet->e3_b1;
    iVals[lBase+i0*lNLabel+17] = pAddJet->e3_v1_b1;
    iVals[lBase+i0*lNLabel+18] = pAddJet->e3_v2_b1;
    iVals[lBase+i0*lNLabel+19] = pAddJet->e4_v1_b1;
    iVals[lBase+i0*lNLabel+20] = pAddJet->e4_v2_b1;
    iVals[lBase+i0*lNLabel+21] = pAddJet->e2_b2;
    iVals[lBase+i0*lNLabel+22] = pAddJet->e3_b2;
    iVals[lBase+i0*lNLabel+23] = pAddJet->e3_v1_b2;
    iVals[lBase+i0*lNLabel+24] = pAddJet->e3_v2_b2;
    iVals[lBase+i0*lNLabel+25] = pAddJet->e4_v1_b2;
    iVals[lBase+i0*lNLabel+26] = pAddJet->e4_v2_b2;
    iVals[lBase+i0*lNLabel+27] = pAddJet->e2_sdb1;
    iVals[lBase+i0*lNLabel+28] = pAddJet->e3_sdb1;
    iVals[lBase+i0*lNLabel+29] = pAddJet->e3_v1_sdb1;
    iVals[lBase+i0*lNLabel+30] = pAddJet->e3_v2_sdb1;
    iVals[lBase+i0*lNLabel+31] = pAddJet->e4_v1_sdb1;
    iVals[lBase+i0*lNLabel+32] = pAddJet->e4_v2_sdb1;
    iVals[lBase+i0*lNLabel+33] = pAddJet->e2_sdb2;
    iVals[lBase+i0*lNLabel+34] = pAddJet->e3_sdb2;
    iVals[lBase+i0*lNLabel+35] = pAddJet->e3_v1_sdb2;
    iVals[lBase+i0*lNLabel+36] = pAddJet->e3_v2_sdb2;
    iVals[lBase+i0*lNLabel+37] = pAddJet->e4_v1_sdb2;
    iVals[lBase+i0*lNLabel+38] = pAddJet->e4_v2_sdb2;
    iVals[lBase+i0*lNLabel+39] = pAddJet->e3_v2_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+40] = pAddJet->e3_v2_sdb2/(pAddJet->e2_sdb2*pAddJet->e2_sdb2);
    iVals[lBase+i0*lNLabel+41] = pAddJet->e3_v1_sdb1/(pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+42] = pAddJet->e3_v1_sdb2/(pAddJet->e2_sdb2);
    iVals[lBase+i0*lNLabel+43] = pAddJet->e3_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+44] = pAddJet->e3_sdb2/(pAddJet->e2_sdb2*pAddJet->e2_sdb2*pAddJet->e2_sdb2);
    iVals[lBase+i0*lNLabel+45] = pAddJet->e3_v2_b1/(pAddJet->e2_b1*pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+46] = pAddJet->e3_v2_b2/(pAddJet->e2_b2*pAddJet->e2_b2);
    iVals[lBase+i0*lNLabel+47] = pAddJet->e3_v1_b1/(pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+48] = pAddJet->e3_v1_b2/(pAddJet->e2_b2);
    iVals[lBase+i0*lNLabel+49] = pAddJet->e3_b1/(pAddJet->e2_b1*pAddJet->e2_b1*pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+50] = pAddJet->e3_b2/(pAddJet->e2_b2*pAddJet->e2_b2*pAddJet->e2_b2);

    fpartonFlavor   = iObjects[0]->partonFlavor;
    fhadronFlavor   = iObjects[0]->hadronFlavor;
    fnCharged       = iObjects[0]->nCharged;
    fnNeutrals      = iObjects[0]->nNeutrals;
    fnParticles     = iObjects[0]->nParticles;

    // SubJets
    int lNSubJets(0);
    TLorentzVector vJ;   vJ.SetPtEtaPhiM(iObjects[i0]->pt,iObjects[i0]->eta,iObjects[i0]->phi,iObjects[i0]->mass);
    TLorentzVector vSJ1; vSJ1.SetPtEtaPhiM(pAddJet->sj1_pt, pAddJet->sj1_eta, pAddJet->sj1_phi, pAddJet->sj1_m); if(pAddJet->sj1_pt>0){ lNSubJets++; fGoodVSubJets.push_back(vSJ1);}
    TLorentzVector vSJ2; vSJ2.SetPtEtaPhiM(pAddJet->sj2_pt, pAddJet->sj2_eta, pAddJet->sj2_phi, pAddJet->sj2_m); if(pAddJet->sj2_pt>0){ lNSubJets++; fGoodVSubJets.push_back(vSJ2);}
    TLorentzVector vSJ3; vSJ3.SetPtEtaPhiM(pAddJet->sj3_pt, pAddJet->sj3_eta, pAddJet->sj3_phi, pAddJet->sj3_m); if(pAddJet->sj3_pt>0){ lNSubJets++; fGoodVSubJets.push_back(vSJ3);}
    TLorentzVector vSJ4; vSJ4.SetPtEtaPhiM(pAddJet->sj4_pt, pAddJet->sj4_eta, pAddJet->sj4_phi, pAddJet->sj4_m); if(pAddJet->sj4_pt>0){ lNSubJets++; fGoodVSubJets.push_back(vSJ4);}

    float sj12dR(-1), sj13dR(-1), sj23dR(-1), sj14dR(-1), sj24dR(-1), sj34dR(-1);
    float sj12dPhi(-1), sj13dPhi(-1),sj23dPhi(-1), sj14dPhi(-1), sj24dPhi(-1), sj34dPhi(-1);
    float sj12dPhiJRF(-1), sj13dPhiJRF(-1),sj23dPhiJRF(-1), sj14dPhiJRF(-1), sj24dPhiJRF(-1), sj34dPhiJRF(-1);

    if (lNSubJets>1) {
      sj12dR = vSJ1.DeltaR(vSJ2);  sj12dPhi = dPhi(vSJ1,vSJ2,vSJ1+vSJ2);  sj12dPhiJRF = dPhi(vSJ1,vSJ2,vJ);
      if (lNSubJets>2) {
	sj13dR = vSJ1.DeltaR(vSJ3);  sj13dPhi = dPhi(vSJ1,vSJ3,vSJ1+vSJ3); sj13dPhiJRF = dPhi(vSJ1,vSJ3,vJ);  
	sj23dR = vSJ2.DeltaR(vSJ3);  sj23dPhi = dPhi(vSJ2,vSJ3,vSJ2+vSJ3); sj23dPhiJRF = dPhi(vSJ2,vSJ3,vJ);
	if (lNSubJets>3) {
	  sj14dR = vSJ1.DeltaR(vSJ4);  sj14dPhi = dPhi(vSJ1,vSJ4,vSJ1+vSJ4); sj14dPhiJRF = dPhi(vSJ1,vSJ4,vJ); 
	  sj24dR = vSJ2.DeltaR(vSJ4);  sj24dPhi = dPhi(vSJ2,vSJ4,vSJ2+vSJ4); sj24dPhiJRF = dPhi(vSJ2,vSJ4,vJ);
	  sj34dR = vSJ3.DeltaR(vSJ4);  sj34dPhi = dPhi(vSJ3,vSJ4,vSJ3+vSJ4); sj34dPhiJRF = dPhi(vSJ3,vSJ4,vJ);
	}
      }
    }

    std::vector<sjpair> sjpairs;
    sjpairs.push_back(sjpair(clean(sj12dR),clean(sj12dPhi),clean(sj12dPhiJRF)));
    sjpairs.push_back(sjpair(clean(sj13dR),clean(sj13dPhi),clean(sj13dPhiJRF)));
    sjpairs.push_back(sjpair(clean(sj23dR),clean(sj23dPhi),clean(sj23dPhiJRF)));
    sjpairs.push_back(sjpair(clean(sj14dR),clean(sj14dPhi),clean(sj14dPhiJRF)));
    sjpairs.push_back(sjpair(clean(sj24dR),clean(sj24dPhi),clean(sj24dPhiJRF)));
    sjpairs.push_back(sjpair(clean(sj34dR),clean(sj34dPhi),clean(sj34dPhiJRF)));

    float ldR_sj0dR(999), ldPhi_sj0dPhi(-999), ldPhiJRF_sj0dPhiJRF(-999);
    std::sort(sjpairs.begin(),sjpairs.end(),compsjpairs);
    ldR_sj0dR = sjpairs[0].dR;
    std::sort(sjpairs.begin(),sjpairs.end(),compsjpairsdPhi);
    ldPhi_sj0dPhi = sjpairs[0].dPhi;
    std::sort(sjpairs.begin(),sjpairs.end(),compsjpairsdPhiJRF);
    ldPhiJRF_sj0dPhiJRF = sjpairs[0].dPhiJRF;

    fdR_sj0dR           = ldR_sj0dR;
    fdPhi_sj0dPhi       = ldPhi_sj0dPhi;
    fdPhiJRF_sj0dPhiJRF = ldPhiJRF_sj0dPhiJRF;
  }
}
int VJetLoader::getMatchedCHSJetIndex(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR) {  
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
void VJetLoader::matchJet(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR, int jIndex){
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

void VJetLoader::matchJet15(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR){
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
void VJetLoader::fillVJetCHS(TJet *iJet, int jIndex){
  TAddJet *pAddJet = getAddJetCHS(iJet);
  fdoublecsvCHS[jIndex] = double(pAddJet->doublecsv);
  fdoublesubCHS[jIndex] = double(pAddJet->Double_sub);
  fptCHS[jIndex] = double(iJet->pt);
  fetaCHS[jIndex] = double(iJet->eta);
  fphiCHS[jIndex] = double(iJet->phi);
}
void VJetLoader::addSubJetBTag(std::string iHeader,TTree *iTree,std::string iLabel,std::vector<std::string> &iLabels,int iN,std::vector<float> &iVals) {
  int iBase=iN;
  for(int i0 = 0; i0 < int(iLabels.size()); i0++) {
    std::stringstream pVal0,pVal1,pValminus1,pVal2;
    pVal0       << iHeader << "btagw" << iLabel << "0"      << "_" << iLabels[i0 % iLabels.size()]; //btagwL0_CENT -- where iLabel(L,M) and iLabels(CENT,MISTAGUP,MISTAGD0,BTAGUP,BTAGD0)
    pVal1       << iHeader << "btagw" << iLabel << "1"      << "_" << iLabels[i0 % iLabels.size()]; //btagwL1_CENT
    pValminus1  << iHeader << "btagw" << iLabel << "minus1" << "_" << iLabels[i0 % iLabels.size()]; //btagwLminus1_CENT
    pVal2       << iHeader << "btagw" << iLabel << "2"      << "_" << iLabels[i0 % iLabels.size()]; //btagwL2_CENT
    iTree->Branch(pVal0      .str().c_str(),&iVals[iBase+0],(pVal0      .str()+"/F").c_str());
    iTree->Branch(pVal1      .str().c_str(),&iVals[iBase+1],(pVal1      .str()+"/F").c_str());
    iTree->Branch(pValminus1 .str().c_str(),&iVals[iBase+2],(pValminus1 .str()+"/F").c_str());
    iTree->Branch(pVal2      .str().c_str(),&iVals[iBase+3],(pVal2      .str()+"/F").c_str());
    iBase+=4;
  }
}
void VJetLoader::fillSubJetBTag(const TClonesArray* iGens, std::vector<TLorentzVector> iObjects) {
  // vSFL should contain CENT (), MISTAG(Ms), BTAG(Bs)  - 5 - CENT(vSFL.at(0)),MsUP(vSFL.at(1)),MsDO(vSFL.at(2)),BsUP(vSFL.at(3)),BsDO(vSFL.at(4))
  int iN = 0;
  for(unsigned int j0=0; j0<2; j0++){  // L, M
    std::vector<std::vector<float>> vSFL,vSFL_nominal;
    vSFL_nominal.clear(); vSFL_nominal.push_back(getSubJetSFs("nominal",iGens,iObjects, fSubJetreaders[j0].at(0), fSubJetreaders[j0].at(3)));
    vSFL.clear(); vSFL.push_back(vSFL_nominal.at(0));

    for(auto iftype :flavorTypes) {
      vSFL_nominal.push_back(getSubJetSFs(iftype,iGens,iObjects, fSubJetreaders[j0].at(0), fSubJetreaders[j0].at(3))); // 0 and 3 HF and LF respectively - flavor types: Ms,Bs
    }

    for(auto iftype :flavorTypes) {
      for(unsigned int i0=1; i0<3; i0++){
	std::vector<float> vSF0; vSF0.clear();
        for(unsigned int i1=0; i1<(vSFL.at(0)).size(); i1++) {
          if(iftype.compare("Ms")==0) vSF0.push_back( (getSubJetSFs(iftype,iGens,iObjects, fSubJetreaders[j0].at(i0), fSubJetreaders[j0].at(i0+3))).at(i1) * (vSFL_nominal.at(2).at(i1)));
          if(iftype.compare("Bs")==0) vSF0.push_back( (getSubJetSFs(iftype,iGens,iObjects, fSubJetreaders[j0].at(i0), fSubJetreaders[j0].at(i0+3))).at(i1) * (vSFL_nominal.at(1).at(i1)));
        }
        vSFL.push_back(vSF0);
      }
    }

    // Fill SubJet btag
    for(unsigned int j1=0; j1<5; j1++){
      int lBase = j1*4+iN;
      fSubJetBTagVars[lBase+0] = getSubJetBtagEventReweight(iGens, 0,  iObjects, vSFL.at(j1));
      fSubJetBTagVars[lBase+1] = getSubJetBtagEventReweight(iGens, 1,  iObjects, vSFL.at(j1));
      fSubJetBTagVars[lBase+2] = getSubJetBtagEventReweight(iGens, -1, iObjects, vSFL.at(j1));
      fSubJetBTagVars[lBase+3] = getSubJetBtagEventReweight(iGens, 2,  iObjects, vSFL.at(j1));
    }
    iN += 20;
  }

}
TAddJet *VJetLoader::getAddJet(TJet *iJet) { 
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
TAddJet *VJetLoader::getAddJetCHS(TJet *iJet) {
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
double VJetLoader::dPhi(TLorentzVector v1, TLorentzVector v2, TLorentzVector v3){
  TVector3 hVelocity = v3.BoostVector();
  TLorentzRotation Boost(hVelocity);
  TLorentzRotation tosubjetRest = Boost.Inverse();
  TVector3 v1Dir = (tosubjetRest*v1).Vect();
  TVector3 v2Dir = (tosubjetRest*v2).Vect();

  float angle_mht = atan2((v1.Py()+v2.Py()),(v1.Px()+v2.Px())) + TMath::Pi();

  Float_t b1 = TMath::Pi() - acos(cos(v1Dir.Phi()-angle_mht));
  Float_t b2 = TMath::Pi() - acos(cos(v2Dir.Phi()-angle_mht));

  return TMath::Min(b1,b2);
}
