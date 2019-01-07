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

  // only valid for 2016 pu weight
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

  fkfactorQCD   = 1;
  fkfactorEWK   = 1;

  fMet          = 0; 
  fMetPhi       = 0; 
  fPuppEt       = 0; 
  fPuppEtPhi    = 0; 
}
void EvtLoader::setupTree(TTree *iTree) {
  reset();
  fTree = iTree;
  fTree->Branch("runNum"          ,&fRun            ,"fRun/i");
  fTree->Branch("lumiSec"         ,&fLumi           ,"fLumi/i");
  fTree->Branch("evtNum"          ,&fEvtV           ,"fEvtV/i");
  fTree->Branch("passJson"        ,&fPassJson       ,"fPassJson/i");
  fTree->Branch("metfilter"       ,&fMetFilters     ,"fMetFilters/i");
  fTree->Branch("triggerBits"     ,&fITrigger       ,"fITrigger/l");
 
  fTree->Branch("npu"             ,&fNPU            ,"fNPU/i");
  fTree->Branch("npv"             ,&fNVtx           ,"fNVtx/i");
  fTree->Branch("rho"             ,&fRho            ,"fRho/F");
  fTree->Branch("puWeight"        ,&fPUWeight       ,"fPUWeight/F");
  fTree->Branch("scale1fb"        ,&fScale          ,"fScale/F");  
  fTree->Branch("evtWeight"       ,&fevtWeight      ,"fevtWeight/F");

  fTree->Branch("kfactorQCD"      ,&fkfactorQCD     ,"fkfactorQCD/F");
  fTree->Branch("kfactorEWK"      ,&fkfactorEWK     ,"fkfactorEWK/F");

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
  fkfactorQCD   = 1;
  fkfactorEWK   = 1;
  fMet          = fEvt->pfMETC;
  fMetPhi       = fEvt->pfMETCphi;
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
// Trigger
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
ULong64_t EvtLoader::triggerBit() {
  ULong64_t lBit = 0;
  for(ULong64_t i0 = 0; i0 < fTrigString.size(); i0++) { 
    if(fTrigger->pass(fTrigString[i0],fEvt->triggerBits)) {
      ULong64_t lBittmp(0),lBitor(0); 
      lBittmp = 1L << i0;
      lBitor = (lBittmp |= lBit);
      lBit = lBitor;
    }
  }
  return lBit;
}
// puWeight 2016
float EvtLoader::puWeight(float iNPU) { 
  float lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(iNPU)));
  if(iNPU > 50) lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(50)));
  if(iNPU <  1) lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(0)));
  return lNPVW;
}
// kFactors
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

  fHist2->Divide(fHist1); // EWK correction = NLO*ewk/LO
  fHist0->Divide(fHist1); // Old QCD correction

  int iPtMin(200), iPtMax(1000);
  fkfactorQCD = Float_t(fHist0->GetBinContent(fHist0->FindBin(iPt)));
  if(iPt > iPtMax) fkfactorQCD = Float_t(fHist0->GetBinContent(fHist0->FindBin(iPtMax)));
  if(iPt < iPtMin) fkfactorQCD = Float_t(fHist0->GetBinContent(fHist0->FindBin(iPtMin)));
  
  fkfactorEWK = Float_t(fHist2->GetBinContent(fHist2->FindBin(iPt)));
  if(iPt > iPtMax) fkfactorEWK = Float_t(fHist2->GetBinContent(fHist2->FindBin(iPtMax)));
  if(iPt < iPtMin) fkfactorEWK = Float_t(fHist2->GetBinContent(fHist2->FindBin(iPtMin)));
}
