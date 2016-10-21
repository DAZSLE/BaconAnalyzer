//============================================================================================
// B-tagging Uncertainty functions:
//-----------------------------
//  - jet and subjet b-tagging SFs
//  - b-tag MC efficiencies for b-jets, c-jets, and light flavor(u/d/s/g) jets
//  - b-tag event weight depending on # of selected b-tags in event
//============================================================================================
#ifndef BTAGUNC_HH
#define BTAGUNC_HH

// bacon object headers
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

// B-tag calibration and SF headers
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "BTagCalibrationStandalone.h"

// ROOT headers                                                                                                                                                                                                                  
#include <TLorentzVector.h>
#include <TMath.h>
#include <TClonesArray.h>

// C++ headers                                                                                                                                             
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <cassert>

//=== FUNCTION DECLARATIONS ====================================================================================== 
std::vector<float> getSubJetSFs(std::string flavor, const TClonesArray* genParArr, std::vector<TLorentzVector> vGoodSubJet, BTagCalibrationReader *HFSJreader, BTagCalibrationReader *LFSJreader);
float getSubJetBtagEventReweight(const TClonesArray* genParArr, int NminBjets, std::vector <TLorentzVector> &vSubJet, std::vector <float> SF);

#endif
