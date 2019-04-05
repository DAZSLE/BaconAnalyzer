#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TPSWeight.hh"

using namespace baconhep;

class GenLoader { 
public:
  GenLoader(TTree *iTree, bool isPs=false);
  ~GenLoader();
  void reset();
  void setupTree(TTree *iTree,float iXSIn = 1);
  void resetHiggs();
  void setupTreeHiggs(TTree *iTree);
  void load (int iEvent);

  // Helpers
  TGenParticle* getParticle(int i);
  bool isType(std::string boson,std::string mode);
  bool hard(int &iP);
  bool hasChild(int &iparent, bool beHard=false);
  TGenParticle* findDaughter(int iparent, int dauId);
  int findDaughterId(int iparent, int dauId);
  int findLastParent(int iparent,int iId);
  void findBoson(int iId, int lOption);

  // Matching
  int getHadronicWInTopFlavor(TGenParticle *genp,int iW,TLorentzVector jet,double dR,double &topMatching, double &topSize);
  int getHadronicVflav(TGenParticle *genp,int j,int iId, TLorentzVector jet,double dR,double &vMatching, double &vSize);
  int ismatchedJet(TLorentzVector jet0, double dR,double &matching, double &size, int iId = 6);
  int ismatchedSubJet(TLorentzVector subjet0);
  int isHWWsemilepBoson(int iH,int iId,int iIdDau,float &genSize);
  int isHDau(int iId, int iDauId, TLorentzVector jet, int iHiggs);

  // tt
  int isttbarType(int lepPdgId);
  void saveTTbarType();
  float computeTTbarCorr();

  // ps weights (only for some samples)
  void setPSWeights(TTree *iTree);
  void loadPSWeights(int iEvent);
  void fillPSWeights();

  TClonesArray  *fGens;
  TBranch       *fGenBr;
  TGenEventInfo *fGenInfo;
  TBranch       *fGenInfoBr;
  TClonesArray  *fPSWeights;
  TBranch       *fPSWeightBr;
  float fWeight;
  float fBosonPt;
  float fBosonPhi;
  float fBosonMass;
  float fBosonEta;
  int fBosonPdgId;
  int genEleFromW;
  int genMuFromW;
  int genTauFromW;
  float fTopPt;
  float fAntitopPt;
  float fTopPtWeight;

  float fgenbPt,fgenbPhi,fgenbMass,fgenbEta;
  float fgenlPt,fgenlPhi,fgenlEta;
  int fgenlId;

  std::vector<float> fgenHPt;
  std::vector<float> fgenHEta;
  std::vector<float> fgenHPhi;
  std::vector<float> fgenHMass;

  std::vector<std::vector<float>> fgenHDauPt;
  std::vector<std::vector<float>> fgenHDauEta;
  std::vector<std::vector<float>> fgenHDauPhi;
  std::vector<std::vector<float>> fgenHDauM;
  std::vector<std::vector<float>> fgenHDauId;  
  std::vector<std::vector<float>> fgenHDauDecay;

  int fNWeights = 50; // # of weights
  std::vector<float> fgenPSWeight;
protected: 
  TTree         *fTree;
};
