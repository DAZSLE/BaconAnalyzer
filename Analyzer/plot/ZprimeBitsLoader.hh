#ifndef ZPRIMEBITSLOADER_HH
#define ZPRIMEBITSLOADER_HH

#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting styles
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TH1D.h>                     // 1D histogram class
#include <TLorentzVector.h>           // 4-vector class
#include <vector>                     // STL vector class
#include <iostream>                   // standard I/O
#include <iomanip>                    // functions to format standard I/O
#include <fstream>                    // functions for file I/O
#include <string>                     // C++ string class
#include <cmath>                      // C++ math library
#include <cassert>

#include "CPlot.hh"                   // helper class for plots
#include "CSample.hh"                 // helper class to manage samples

using namespace std;

class ZprimeBitsLoader {
public:
  ZprimeBitsLoader(TTree *iTree=0,TString algo="PUPPI", TString jet="jet0",TString number="8");		
  ~ZprimeBitsLoader();
  bool isPho(bool isData);
  bool isHad(bool isData);
  bool selectJetAlgoAndSize(TString algo);
  bool passPreSelection(bool isData, string preselection);
  bool passBoostedZprimePreSelection();
  bool passBoostedGammaZprimeSelection();
  bool passBoostedZprimeSR(float ddtcut);
  bool passBoostedGammaZprimeSR(float ddtcut);
  bool passSelection(bool isData,string selection,string subsample,float ddt,float csv1);
  bool passBoostedZprimeBTag(float csvcut);
  double getWgt(bool isData, TString algo, double LUMI);
  double tau21DDT();
  double N2(double beta);
  double N2sd(double beta);
  double M2(double beta);
  double M2sd(double beta);
  double D2(double beta);
  double D2sd(double beta);
  
  //
  // variables to read in bacon bits
  //
  int metfilter;                                                                   // MET filter bits          
  unsigned int runNum, lumiSec, evtNum;                                            // event ID
  unsigned int triggerBits, selectBits;                                            // trigger and jet type bits
  double triggerEff;                                                               // trigger efficiency
  float evtWeight, puWeight;                                                       // pu and evt weight
  unsigned int npu, npv;                                                           // PU, PV multiplicity
  int nmu, nele, ntau, npho;                                                       // object multiplicity
  int njets;                                                                       // jet multiplicity 
  float scale1fb;                                                                  // cross section scale factor per 1/fb
  float vmetpt,vmetphi,vfakemetpt,vfakemetphi;                                     // MET and fake MET

  //
  // Jets with pT-ordering:
  //
  double           bst_jet0_pt,bst_jet0_eta,bst_jet0_phi,bst_jet0_mass;            // leading boosted jet
  double           bst_jet0_msd, bst_jet0_rho, bst_jet0_phil;                      // msd and msd/pt dependent variables
  double           bst_jet0_tau32, bst_jet0_tau21;                                 // boosted tagger variables                                                                                            
  double           bst_jet0_maxsubcsv, bst_jet0_minsubcsv;                         // subjet btag                            
  double           bst_jet0_CHF, bst_jet0_NHF, bst_jet0_NEMF;                      // boosted jet variables
  float            bst_mt;                                                         // mT
  double           bst_jet0_doublecsv;                                             // doublecsv variable for leading jet
  double           e2_b1, e3_b1, e3_v1_b1, e3_v2_b1;
  double           e2_b2, e3_b2, e3_v1_b2, e3_v2_b2;
  double           e2_sdb1, e3_sdb1, e3_v1_sdb1, e3_v2_sdb1;
  double           e2_sdb2, e3_sdb2, e3_v1_sdb2, e3_v2_sdb2;
protected:

  const unsigned int kBOOSTED8PUPPI = 2;
  const float RHO_CUT = 0.38;
  const unsigned int kMET  = 2;
  const unsigned int kHT  = 2;
  const unsigned int kSinglePhoton  = 4;

};
#endif
