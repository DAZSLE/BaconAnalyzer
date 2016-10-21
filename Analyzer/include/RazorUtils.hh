
//================================================================================================
//
// Various functions for including razor variables 
//
//________________________________________________________________________________________________
#ifndef RAZORUTILS_HH
#define RAZORUTILS_HH

// bacon object headers
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
//Root
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>
//STL
#include <vector>
#include <cassert>
#include <iostream>

#include "Hemisphere.hh"
#include "Davismt2.hh"

//=== FUNCTION DECLARATIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
double GetAlphaT(vector<TLorentzVector> jets) ;
double GetDPhiMin(vector<TLorentzVector> jets);
vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets);
std::vector< std::vector<int> > getHemispheresV2( std::vector<TLorentzVector> jets);
double computeMR(TLorentzVector hem1, TLorentzVector hem2);
double computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector met);
double calcMT2(float testMass, bool massive, std::vector<TLorentzVector> jets, TLorentzVector MET, int hemi_seed, int hemi_association);

#endif
