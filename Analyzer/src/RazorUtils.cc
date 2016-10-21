#include "../include/RazorUtils.hh"
//================================================================================================
//
// Various functions for including razor variables
//
//________________________________________________________________________________________________

// bacon object headers
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#include <string>
#include <sstream>
#include <vector>

//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
double GetAlphaT(vector<TLorentzVector> jets) 
{   
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;
    double alphaT = 0;
    
    if(nJets < 2) return alphaT;
    
    int nComb = pow(2, nJets); // # possible combinations
    
    // steal from the getHemispheres method

    //step 1: store all possible partitions of the input jets
    int j_count;
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;
        
        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }
    
    //step 2: Select combination that mininize |ET1 - ET2|
    double eMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;

    for(size_t i=0; i < possibleHem1s.size(); i++)
    {
        double eTemp = fabs(possibleHem1s[i].Et() - possibleHem2s[i].Et());
        if (eMin < 0 || eTemp < eMin)
        {
            eMin = eTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }
    
    float MhtX = 0., MhtY = 0.;
    float HT = 0.; 
    for (auto& obj : jets) { HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); }

      TLorentzVector MyMHT;
      MyMHT.SetPxPyPzE(-MhtX, -MhtY, 0, sqrt(pow(MhtX,2) + pow(MhtY,2)));

    float MHT = MyMHT.Pt();

    // Calculate alphaT
    alphaT = 0.5 * (1-eMin/HT)/sqrt(1-pow(MHT/HT,2));

    return alphaT;  
}

//--------------------------------------------------------------------------------------------------
double GetDPhiMin(vector<TLorentzVector> jets)
    // This variable is used in the alphaT analysis
{
  //    int nJets = jets.size();
    double dPhiMin = -1.;
    float HT = 0.;
    float MhtX = 0.;
    float MhtY = 0.;
    // Search for min dPhi between recomputed missing HT and test jets
    for (auto& obj : jets) { HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); }
    TLorentzVector MyMHT;
    MyMHT.SetPxPyPzE(-MhtX, -MhtY, 0, sqrt(pow(MhtX,2) + pow(MhtY,2)));

    for (auto& obj : jets)
    {
    // Recompute MHT by ignoring a test jet 
        float recomputedMHTX = MhtX - obj.Px();
        float recomputedMHTY = MhtY - obj.Py();
        TLorentzVector recomputedMHT;
        recomputedMHT.SetPxPyPzE(-recomputedMHTX, -recomputedMHTY, 0, sqrt(pow(recomputedMHTX,2) + pow(recomputedMHTY,2)));
        double phiTemp = fabs(recomputedMHT.Phi() - obj.Phi());
        if (dPhiMin < 0 || phiTemp < dPhiMin)   dPhiMin = phiTemp;
    }

    return dPhiMin;
}

//--------------------------------------------------------------------------------------------------
vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets){
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;

    if(nJets < 2){ //return empty hemispheres if there are fewer than 2 jets provided
        TLorentzVector emptyHem1, emptyHem2;
        vector<TLorentzVector> emptyHemsOut;
        emptyHemsOut.push_back(emptyHem1);
        emptyHemsOut.push_back(emptyHem2);
        return emptyHemsOut;
    }

    //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
    int nComb = pow(2, nJets);

    //step 1: store all possible partitions of the input jets
    int j_count;
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;
        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }

    //step 2: choose the partition that minimizes m1^2 + m2^2
    double mMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;
    for(size_t i=0; i < possibleHem1s.size(); i++){
        double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
        if(mMin < 0 || mTemp < mMin){
            mMin = mTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }

    //return the hemispheres in decreasing order of pt
    vector<TLorentzVector> hemsOut;
    if(myHem1.Pt() > myHem2.Pt()){
        hemsOut.push_back(myHem1);
        hemsOut.push_back(myHem2);
    } else {
        hemsOut.push_back(myHem2);
        hemsOut.push_back(myHem1);
    }

    return hemsOut;
}

//--------------------------------------------------------------------------------------------------
std::vector< std::vector<int> > getHemispheresV2( std::vector<TLorentzVector> jets )
{
  //returns vector with original indices to jets
  int nJets = jets.size();
  vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations                                                              
  std::vector< std::vector<int> > index1;
  vector<TLorentzVector> possibleHem2s;
  std::vector< std::vector<int> > index2;

  if(nJets < 2){ //return empty hemispheres if there are fewer than 2 jets provided                                                           
    std::vector<int> emptyIndex1, emptyIndex2;
    std::vector< std::vector<int> > void_return;
    void_return.push_back( emptyIndex1 );
    void_return.push_back( emptyIndex2 );
    return void_return;
  }
  
  //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc                                              
  int nComb = pow(2, nJets);
  //std::cout << "njets: " << nJets << " ncomb: " << nComb << std::endl;
  //step 1: store all possible partitions of the input jets                                                                                   
  int j_count;
  for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
    //std::cout << "=iter: " << i << std::endl; 
    TLorentzVector j_temp1, j_temp2;
    std::vector<int> tmp_index1, tmp_index2;
    int itemp = i;
    j_count = nComb/2;
    int count = 0;
    while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2               
      //std::cout << "j_count: " << j_count << " itemp: " << itemp << " count: " << count << std::endl; 
      if(itemp/j_count == 1){
	j_temp1 += jets[count];
	tmp_index1.push_back( count );
      } else {
	j_temp2 += jets[count];
	tmp_index2.push_back( count );
      }
      itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count                                                            
      j_count /= 2;
      count++;
    }
    possibleHem1s.push_back(j_temp1);
    index1.push_back( tmp_index1 );
    possibleHem2s.push_back(j_temp2);
    index2.push_back( tmp_index2 );
  }
  
  //step 2: choose the partition that minimizes m1^2 + m2^2                                                                                   
  double mMin = -1;
  TLorentzVector myHem1;
  TLorentzVector myHem2;
  int partition_index = -1;
  for(size_t i=0; i < possibleHem1s.size(); i++){
    double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
    if(mMin < 0 || mTemp < mMin){
      mMin = mTemp;
      myHem1 = possibleHem1s[i];
      myHem2 = possibleHem2s[i];
      partition_index = i;
    }
  }

  //return the hemispheres in decreasing order of pt                                                                                          
  vector<TLorentzVector> hemsOut;
  std::vector< std::vector<int> > index_out;
  if(myHem1.Pt() > myHem2.Pt()){
    hemsOut.push_back(myHem1);
    hemsOut.push_back(myHem2);
    index_out.push_back( index1[partition_index] );
    index_out.push_back( index2[partition_index] );
  } else {
    hemsOut.push_back(myHem2);
    hemsOut.push_back(myHem1);
    index_out.push_back( index2[partition_index] );
    index_out.push_back( index1[partition_index] );
  }
  
  return index_out;
};

//--------------------------------------------------------------------------------------------------
double computeMR(TLorentzVector hem1, TLorentzVector hem2){
    return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

//--------------------------------------------------------------------------------------------------
double computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet){
    double mR = computeMR(hem1, hem2);
    double term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
    double term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
    double mTR = sqrt(term1 - term2);
    return (mTR / mR) * (mTR / mR);
}

//--------------------------------------------------------------------------------------------------
double calcMT2(float testMass, bool massive, vector<TLorentzVector> jets, TLorentzVector MET, int hemi_seed, int hemi_association)
{
  //computes MT2 using a test mass of testMass, with hemispheres made massless if massive is set to false
  //hemispheres are clustered by finding the grouping of input jets that minimizes the Lund distance
  
  if(jets.size() < 2) return -9999; //need at least two jets for the calculation
  vector<float> px, py, pz, E;
  for(uint i = 0; i < jets.size(); i++){
    //push 4vector components onto individual lists
    px.push_back(jets[i].Px());
    py.push_back(jets[i].Py());
    pz.push_back(jets[i].Pz());
    E.push_back(jets[i].E());
  }
  
  //form the hemispheres using the provided Hemisphere class
  Hemisphere* hemis = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
  vector<int> grouping = hemis->getGrouping();
  TLorentzVector pseudojet1(0.,0.,0.,0.);
  TLorentzVector pseudojet2(0.,0.,0.,0.);
        
  //make the hemisphere vectors
  for(uint i=0; i<px.size(); ++i){
    if(grouping[i]==1){
      pseudojet1.SetPx(pseudojet1.Px() + px[i]);
      pseudojet1.SetPy(pseudojet1.Py() + py[i]);
      pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
      pseudojet1.SetE( pseudojet1.E()  + E[i]);   
    }else if(grouping[i] == 2){
      pseudojet2.SetPx(pseudojet2.Px() + px[i]);
      pseudojet2.SetPy(pseudojet2.Py() + py[i]);
      pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
      pseudojet2.SetE( pseudojet2.E()  + E[i]);
    }
  }
  delete hemis;
  
  //now compute MT2 using the Davismt2 class
  
  //these arrays contain (mass, px, py) for the pseudojets and the MET
  double pa[3];
  double pb[3];
  double pmiss[3];
  
  pmiss[0] = 0;
  pmiss[1] = static_cast<double> (MET.Px());
  pmiss[2] = static_cast<double> (MET.Py());
  
  pa[0] = static_cast<double> (massive ? pseudojet1.M() : 0);
  pa[1] = static_cast<double> (pseudojet1.Px());
  pa[2] = static_cast<double> (pseudojet1.Py());
  
  pb[0] = static_cast<double> (massive ? pseudojet2.M() : 0);
  pb[1] = static_cast<double> (pseudojet2.Px());
  pb[2] = static_cast<double> (pseudojet2.Py());
  
  Davismt2 *mt2 = new Davismt2();
  mt2->set_momenta(pa, pb, pmiss);
  mt2->set_mn(testMass);
  Float_t MT2=mt2->get_mt2();
  delete mt2;
  return MT2;
};
