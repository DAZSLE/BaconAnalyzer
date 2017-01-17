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
	    std::string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns",
	    std::string iPUWeight="${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/puWeight_23Sept.root");
  ~EvtLoader(); 
  void reset();
  void setupTree  (TTree *iTree);
  void load (int iEvent);
  //Fillers
  void fillEvent(unsigned int trigBit, float lWeight, unsigned int passJson);
  bool passSkim();
  TLorentzVector Met(int iOption);
  //SFs
  void fillLepSF(std::vector<TLorentzVector> iLeptons, double lepPdgId, int islep0Tight, int islep1Tight);
  //Trigger
  bool passFilter();
  bool passTrigger(std::string iTrigger);
  void triggerEff(std::vector<TLorentzVector> iElectrons, std::vector<TLorentzVector> iPhotons);
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

  float fUW, fUWPhi;
  float fUA, fUAPhi;
  float fUZ, fUZPhi;

  float flep1Id, flep2Id;

  unsigned int fRun;
  unsigned int fEvtV;
  unsigned int fLumi;
  unsigned int fPassJson;
  TEventInfo   *fEvt;

  float fevtWeight;
  float fScale;

  float fkFactor_CENT;
  float fEwkCorr_CENT, fEwkCorr_UP, fEwkCorr_DO;
  float fkfactor;
  float fPDF, fPDF_UP, fPDF_DO;
  float fRenScale_UP, fRenScale_DO, fFacScale_UP, fFacScale_DO;

  double fsf_lep, fsf_lepTrack;

  double fsf_eleTrig;
  double fsf_metTrig;
  double fsf_phoTrig;
  double fLepWeight;

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
  TH1F         *fHistPDF;
  TH1F         *fHistRUP;
  TH1F         *fHistRDO;
  TH1F         *fHistFUP;
  TH1F         *fHistFDO;

  TH2D         *fhMuLoose;
  TH2D         *fhMuTight;
  TH1D         *fhMuTrack;
  TH2D         *fhEleVeto;
  TH2D         *fhEleTight;
  TH2D         *fhEleTrack;
  
  TH1D         *hEleTrigB;
  TH1D         *hEleTrigE;
  TH2D         *hEleTrigLow;
  TH1D         *hMetTrig;
  TH1D         *hPhoTrig;

  std::vector<std::vector<std::string>> fTrigString;

  char*  fSample;
  unsigned int fITrigger;
  unsigned int fMetFilters;
  unsigned int fNPU;
  unsigned int fNPUP;
  unsigned int fNPUM;

  float fPUWeight;
  float fMetSig;
};
