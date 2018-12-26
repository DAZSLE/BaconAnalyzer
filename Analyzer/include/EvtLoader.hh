#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "Utils.hh"

#include "TH1F.h"
#include "TH2D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

using namespace baconhep;

class EvtLoader { 
public:
  EvtLoader(TTree *iTree,std::string iName,
	    std::string iHLTFile="${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/HLTFile_25ns",
	    std::string iPUWeight="${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/puWeights_8X.root");
  ~EvtLoader(); 
  void reset();
  void setupTree  (TTree *iTree);
  void load (int iEvent);
  //Fillers
  void fillEvent(unsigned int trigBit, float lWeight, unsigned int passJson);
  bool passSkim();
  TLorentzVector Met(int iOption);
  //Trigger
  bool passFilter();
  void addTrigger(std::string iName);
  bool passTrigger();
  bool passTrigger(std::string iTrigger);
  unsigned int triggerBit();
  //PU
  float        puWeight(float iPU);
  int          nVtx();
  bool         PV();
  //EWK and kFactor
  void computeCorr(float iPt,
		   std::string iHist0,
		   std::string iHist1,
		   std::string iHist2,
		   std::string iNLO,
                   std::string ikfactor="${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/kfactors.root");
  //MET and mT
  void         fillRecoil(std::vector<TLorentzVector> &iVecCorr,std::vector<TLorentzVector> iPhotons);
  void         fillmT(float iMet, float iMetPhi,float iFMet, float iFMetPhi, std::vector<TLorentzVector> &lCorr, float &fmT);
  float        mT(float iMet,float iMetPhi,TLorentzVector &iVec);
  //Vars
  float fRho;
  float fMetPhi;
  float fMet;
  float fPuppEtPhi;
  float fPuppEt;
  float fCaloMet;
  float fCaloMetPhi;
  float fMetNoMu;

  float fFMetPhi;
  float fFMet;
  float fFPuppEtPhi;
  float fFPuppEt;

  unsigned int fRun;
  unsigned int fEvtV;
  unsigned int fLumi;
  unsigned int fPassJson;
  
  float fPu;

  TEventInfo   *fEvt;

  float fevtWeight;
  float fScale;

  float fkFactor_CENT;
  float fEwkCorr_CENT, fEwkCorr_UP, fEwkCorr_DO;
  float fkfactor;

  int fselectBits;
  int fNVtx;
protected: 
  TBranch      *fEvtBr;

  TClonesArray *fVertices;
  TBranch      *fVertexBr;

  TTree        *fTree;
  TTrigger     *fTrigger;
  TH1F         *fPUWeightHist;

  TH1F         *fHist0;
  TH1F         *fHist1;
  TH1F         *fHist2;

  std::vector<std::string> fTrigString;

  char*  fSample;
  unsigned int fITrigger;
  ULong64_t  fIMoreTrigger;
  unsigned int fMetFilters;
  unsigned int fNPU;
  unsigned int fNPUP;
  unsigned int fNPUM;

  float fPUWeight;
  float fMetSig;
};
