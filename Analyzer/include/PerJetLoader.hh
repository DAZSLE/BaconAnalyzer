#ifndef PerJetLoader_H
#define PerJetLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TSVtx.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "Utils.hh"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JECLoader.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/MeasureDefinition.hh"


#include "TRandom3.h"

using namespace baconhep;

class PerJetLoader { 
public:
  PerJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS="AK8CHS",std::string iAddJetCHS="AddAK8CHS",int iN=1, bool iData=false,std::string iLabel="2017");
  ~PerJetLoader();
  double correction(TJet &iJet,double iRho);
  void reset();
  void setupTreeQbert(TTree *iTree,std::string iJetLabel);
  void load(int iEvent);
  void selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, double metphi, unsigned int runNum);
  void fillVJet(int iN,std::vector<TJet*> &iObjects,double dR, double iRho, double metphi, unsigned int runNum);
  void matchJet(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR, int jIndex);
  TAddJet *getAddJet(TJet *iJet);
  TAddJet *getAddJetCHS(TJet *iJet);

  double fvSize, fvMatching;
  int fisHadronicV;
  std::vector<int> fisTightVJet, fisTightVJetCHS;
  int fpartonFlavor, fhadronFlavor, fnbHadrons, fncHadrons, fnCharged, fnNeutrals, fnParticles, fnVtxFlavor,fnVtxFlavInfo;

  std::vector<TJet*> fLooseVJets, fLooseVJetsCHS;
  std::vector<TLorentzVector> selectedVJets, selectedVJetsCHS;
  std::vector<double> fdoublecsvCHS, fdoublesubCHS, fptCHS, fetaCHS, fphiCHS;
  std::vector<TJet*> fLooseVJetsByDoubleB, fLooseVJetsCHSByDoubleB;
  std::vector<TLorentzVector> selectedVJetsByDoubleB, selectedVJetsCHSByDoubleB;


  const double CSVL = 0.5426; // CSVv2SubJet WP 
  const double CSVM = 0.8484;
  TRandom3* r;

  //bool aktSort = true; 
  bool aktSort = false;

  TClonesArray *fVJets;
  TClonesArray *fPFs;

protected: 
  fastjet::AreaDefinition *areaDef = 0;
  fastjet::GhostedAreaSpec *activeArea = 0;
  fastjet::JetDefinition *jetDef = 0;
  
  //TClonesArray *fVJets;
  TBranch      *fVJetBr;
  TClonesArray *fVAddJets;
  TBranch      *fVAddJetBr;
  TClonesArray *fVJetsCHS;
  TBranch      *fVJetBrCHS;
  TClonesArray *fVAddJetsCHS;
  TBranch      *fVAddJetBrCHS;
  TClonesArray *fGens;
  TBranch      *fGenBr;
  //TClonesArray *fPFs;
  TBranch      *fPFBr;
  TClonesArray *fSVs;
  TBranch      *fSVBr;

  TTree        *fTree;

  JECLoader    *fJEC;

  int           fNLooseVJets, fNLooseVJetsCHS;
  int           fNTightVJets, fNTightVJetsCHS;
  int           fNProngs;

  std::map<std::string,float> fSingletons; 
  int fN_cpf, fN_ipf, fN_sv, fN_part;
	std::map<std::string,float*> fCPFArrs, fIPFArrs, fSVArrs, fPartArrs;
  std::vector<std::string> fTrigString;

  int           fN;
  
  std::string cmsswPath;
  
  // for jet energy corrections
  bool isData;
  std::string fYear;

  // Gaussian random numbers (one for each jet)
  std::vector<double> x1List;
};
#endif
