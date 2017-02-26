#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace baconhep;

class GenLoader { 
public:
  GenLoader(TTree *iTree);
  ~GenLoader();
  void reset();
  void setupTree(TTree *iTree,float iXSIn);
  void load (int iEvent);
  void fillGenEvent();
  //Fill specific Gen Info
  void selectBoson(int iDMu);
  float lepmatched(int iId, std::vector<TLorentzVector> vec, double dR);
  int ismatchedJet(TLorentzVector jet0, double dR,double &matching, double &size, int iId = 6);
  int ismatchedSubJet(TLorentzVector subjet0);
  //Helpers
  bool isType(std::string boson,std::string mode);
  int isttbarType(int lepPdgId);
  void saveTTbarType();
  TGenParticle *fBoson;
  float computeTTbarCorr();
  TGenParticle* findDaughter(int iparent, int dauId);
  int findDaughterId(int iparent, int dauId);
  int findLastBoson(int iparent,int iId);
  int isHadronicWInTop(TGenParticle *genp,int j,TLorentzVector jet,double dR,double &topMatching, double &topSize);
  int isHadronicTop(TGenParticle *genp,int j,TLorentzVector jet,double dR,double &topMatching, double &topSize);
  int isHadronicV(TGenParticle *genp,int j,int iId,TLorentzVector jet,double dR,double &vMatching, double &vSize);
  int isHadronicVflav(TGenParticle *genp,int j,int iId, TLorentzVector jet,double dR,double &vMatching, double &vSize, int dauId);
  void findBoson(int iId, int lOption);

  TClonesArray  *fGens;
  TBranch       *fGenBr;
  TGenEventInfo *fGenInfo;
  TBranch       *fGenInfoBr;

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
  
protected: 
  TTree         *fTree;
};
