#include "TFile.h"
#include "TMatrixD.h"
#include "../include/EvtLoader.hh"
#include <iostream>
#include <sstream>

using namespace baconhep;

EvtLoader::EvtLoader(TTree *iTree,std::string iName,std::string iHLTFile,std::string iPUWeight) { 
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = iHLTFile.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
    iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  fEvt      = new TEventInfo();
  iTree->SetBranchAddress("Info",       &fEvt);
  fEvtBr    = iTree->GetBranch("Info");
  std::cout << "trigger file " << iHLTFile << std::endl;
  fTrigger  = new TTrigger(iHLTFile);
  
  fVertices = new TClonesArray("baconhep::TVertex");
  iTree->SetBranchAddress("PV",       &fVertices);
  fVertexBr     = iTree->GetBranch("PV");
  
  TFile *lFile = TFile::Open(iPUWeight.c_str()); 
  fPUWeightHist =  (TH1F*) lFile->Get("puw");
  fPUWeightHist->SetDirectory(0);
  lFile->Close();
  fSample = (char*) (iName.c_str());

}
EvtLoader::~EvtLoader() { 
  delete  fEvt;
  delete  fEvtBr;
  delete  fVertices;
  delete  fVertexBr;
}
void EvtLoader::reset() { 
  fEvtV         = 0; 
  fITrigger     = 0; 
  fNVtx         = 0; 
  fNPU          = 0;
  fPUWeight     = 0; 
  fScale        = 1;
  fevtWeight    = 1;
  fRho          = 0;

  fkfactor      = 1;
  fkFactor_CENT = 1;
  fEwkCorr_CENT = 1;

  fMet          = 0; 
  fMetPhi       = 0; 
  fPuppEt       = 0; 
  fPuppEtPhi    = 0; 
  fCaloMet      = 0;
  fCaloMetPhi   = 0;
}
void EvtLoader::setupTree(TTree *iTree) {
  reset();
  fTree = iTree;
  fTree->Branch("runNum"          ,&fRun            ,"fRun/i");
  fTree->Branch("lumiSec"         ,&fLumi           ,"fLumi/i");
  fTree->Branch("evtNum"          ,&fEvtV           ,"fEvtV/i");
  fTree->Branch("passJson"        ,&fPassJson       ,"fPassJson/i");
  fTree->Branch("metfilter"       ,&fMetFilters     ,"fMetFilters/i");
  fTree->Branch("triggerBits"     ,&fITrigger       ,"fITrigger/i");
 
  fTree->Branch("npu"             ,&fNPU            ,"fNPU/i");
  fTree->Branch("npv"             ,&fNVtx           ,"fNVtx/i");
  fTree->Branch("rho"             ,&fRho            ,"fRho/F");
  fTree->Branch("puWeight"        ,&fPUWeight       ,"fPUWeight/F");
  fTree->Branch("scale1fb"        ,&fScale          ,"fScale/F");  
  fTree->Branch("evtWeight"       ,&fevtWeight      ,"fevtWeight/F");

  fTree->Branch("kfactor"         ,&fkfactor        ,"fkfactor/F");
  fTree->Branch("kfactorNLO"      ,&fkFactor_CENT   ,"fkFactor_CENT/F");

  fTree->Branch("pfmet"           ,&fMet            ,"fMet/F");
  fTree->Branch("pfmetphi"        ,&fMetPhi         ,"fMetPhi/F");
  fTree->Branch("puppet"          ,&fPuppEt         ,"fPuppEt/F");
  fTree->Branch("puppetphi"       ,&fPuppEtPhi      ,"fPuppEtPhi/F");

  fRun   = 0;
  fLumi  = 0;
  fPassJson = 0;
}
void EvtLoader::load(int iEvent) { 
  fVertices ->Clear();
  fEvtBr    ->GetEntry(iEvent);
  fVertexBr ->GetEntry(iEvent);
  fTrigString.clear();
  fRun   = fEvt->runNum;
  fLumi  = fEvt->lumiSec;
  fPu    = fEvt->nPUmean;
}
void EvtLoader::fillEvent(unsigned int trigBit,float lWeight, unsigned int passJson) { 
  reset();
  fRun          = fEvt->runNum;
  fLumi         = fEvt->lumiSec;
  fNPU          = fEvt->nPUmean;
  fRho          = fEvt->rhoJet;
  fPassJson     = passJson;
  fNVtx         = nVtx();
  fITrigger     = triggerBit();
  fPUWeight     = puWeight(float(fNPU)); 
  fScale        = lWeight;
  fevtWeight    = 1;
  fMetFilters   = fEvt->metFilterFailBits;
  fEvtV         = fEvt->evtNum;
  fkfactor      = 1;
  fMet          = fEvt->pfMETC;
  fMetPhi       = fEvt->pfMETCphi;
  fCaloMet      = fEvt->caloMET;
  fCaloMetPhi   = fEvt->caloMETphi;
  fPuppEt       = fEvt->puppETC;
  fPuppEtPhi    = fEvt->puppETCphi;
  return;
}
// Met Filter
bool EvtLoader::passFilter() { 
  return (fEvt->metFilterFailBits == 0);
}
// Vtx
int EvtLoader::nVtx() {
  unsigned int lNVertex = 0;
  for(int i0 = 0; i0 < fVertices->GetEntries(); i0++) {
    TVertex *pVertex = (TVertex*) ((*fVertices)[i0]);
    if(fabs(pVertex->z) > 24) continue;
    if(pVertex->ndof    <  4) continue;
    float pX = pVertex->x;
    float pY = pVertex->y;
    float pRho  = sqrt(pX*pX+pY*pY);
    if(pRho             > 2.) continue;
    lNVertex++;
  }
  return lNVertex;
}
bool EvtLoader::PV(){
  return fEvt->hasGoodPV;
}
void EvtLoader::addTrigger(std::string iName) { 
  fTrigString.push_back(iName);
}
bool EvtLoader::passTrigger() {
  bool lPass = false;
  for(unsigned int i0 = 0; i0 < fTrigString.size(); i0++) { 
    if(fTrigger->pass(fTrigString[i0],fEvt->triggerBits)) lPass = true;
  }
  return lPass;
}
bool EvtLoader::passTrigger(std::string iTrigger) {
  return fTrigger->pass(iTrigger,fEvt->triggerBits);
}
unsigned int EvtLoader::triggerBit() {
  unsigned int lBit = 0;
  //std::cout << "trig size " <<fTrigString.size() << std::endl;
  for(unsigned int i0 = 0; i0 < fTrigString.size(); i0++) { 
    if(fTrigger->pass(fTrigString[i0],fEvt->triggerBits))  lBit |= 1 << i0;
  }
  return lBit;
}
// puWeight
float EvtLoader::puWeight(float iNPU) { 
  float lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(iNPU)));
  if(iNPU > 50) lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(50)));
  if(iNPU <  1) lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(0)));
  return lNPVW;
}
// mT
float  EvtLoader::mT(float iMet,float iMetPhi,TLorentzVector &iVec) { 
  float lDPhi = fabs(iMetPhi-iVec.Phi());
  if(fabs(lDPhi) > TMath::Pi()*2.-lDPhi) lDPhi = TMath::Pi()*2.-lDPhi;
  float lMt = sqrt(2.0*(iVec.Pt()*iMet*(1.0-cos(lDPhi))));
  return lMt;
}
void  EvtLoader::fillmT(float iMet, float iMetPhi,float iFMet, float iFMetPhi, std::vector<TLorentzVector> &lCorr, float &fmT) {
  if(lCorr.size()>0){
    TLorentzVector lVecCorr;
    for(unsigned int i0 =0; i0 < lCorr.size(); i0++) lVecCorr += lCorr[i0];
    fmT = (lVecCorr.Pt()>0) ?  mT(iMet,     iMetPhi,     lVecCorr): -999;
    if(iFMet>0) fmT = (lVecCorr.Pt()>0) ? mT(iFMet,     iFMetPhi,     lVecCorr): -999;
  }
}
// kFactor and EWK
void EvtLoader::computeCorr(float iPt,std::string iHist0,std::string iHist1,std::string iHist2,std::string iNLO,std::string ikfactor){
  std::string isuffix = "";
  if(iNLO.find("G")!=std::string::npos) isuffix ="_G";

  TFile *lFile = TFile::Open(ikfactor.c_str());
  fHist0 =  (TH1F*) lFile->Get(iHist0.c_str()); // NLO
  fHist0->SetDirectory(0);
  fHist1 =  (TH1F*) lFile->Get(iHist1.c_str()); // LO
  fHist1->SetDirectory(0);
  fHist2 =  (TH1F*) lFile->Get(iHist2.c_str()); // NLO*EWK
  fHist2->SetDirectory(0);
  lFile->Close();

  fHist2->Divide(fHist1);
  fHist0->Divide(fHist1);

  fkFactor_CENT = Float_t(fHist0->GetBinContent(fHist0->FindBin(iPt)));
  //if(iPt > 700) fkFactor_CENT = Float_t(fHist0->GetBinContent(fHist0->FindBin(700)));
  if(iPt < 100) fkFactor_CENT = Float_t(fHist0->GetBinContent(fHist0->FindBin(100)));
  
  fEwkCorr_CENT = Float_t(fHist2->GetBinContent(fHist2->FindBin(iPt)));
  //if(iPt > 700) fEwkCorr_CENT = Float_t(fHist2->GetBinContent(fHist2->FindBin(700)));
  if(iPt < 100) fEwkCorr_CENT = Float_t(fHist2->GetBinContent(fHist2->FindBin(100)));

  fkfactor = fEwkCorr_CENT; //(NLO*ewk/LO)

}
