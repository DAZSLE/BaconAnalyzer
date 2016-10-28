#include "ZprimeBitsLoader.hh"  
using namespace std;

ZprimeBitsLoader::ZprimeBitsLoader(TTree *iTree,TString algo,TString jet,TString number) {
  if(iTree){
    TString met = "puppet"; if (algo!="PUPPI") met = "pfmet";
    iTree->SetBranchAddress("runNum",                            &runNum);
    iTree->SetBranchAddress("lumiSec",                           &lumiSec);
    iTree->SetBranchAddress("evtNum",                            &evtNum);
    iTree->SetBranchAddress("metfilter",                         &metfilter);
    iTree->SetBranchAddress("triggerBits",                       &triggerBits);
    iTree->SetBranchAddress("selectBits",                        &selectBits);
    iTree->SetBranchAddress("npu",                               &npu);
    iTree->SetBranchAddress("npv",                               &npv);
    iTree->SetBranchAddress("nmuLoose",                          &nmu);
    iTree->SetBranchAddress("neleLoose",                         &nele);
    iTree->SetBranchAddress("ntau",                              &ntau);
    iTree->SetBranchAddress("nphoLoose",                         &npho);
    iTree->SetBranchAddress("puWeight",                          &puWeight);
    iTree->SetBranchAddress("scale1fb",                          &scale1fb);
    iTree->SetBranchAddress("evtWeight",                         &evtWeight);
    iTree->SetBranchAddress(met,                                 &vmetpt);
    iTree->SetBranchAddress(met+"phi",                           &vmetphi);
    iTree->SetBranchAddress("nAK8Puppijets",                     &njets);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_pt",             &bst_jet0_pt);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_eta",            &bst_jet0_eta);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_phi",            &bst_jet0_phi);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_msd",            &bst_jet0_msd);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_rho",            &bst_jet0_rho);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_tau21",          &bst_jet0_tau21);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_CHF",            &bst_jet0_CHF);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_NHF",            &bst_jet0_NHF);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_NEMF",           &bst_jet0_NEMF);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_doublecsv",      &bst_jet0_doublecsv);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_minsubcsv",      &bst_jet0_minsubcsv);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_maxsubcsv",      &bst_jet0_maxsubcsv);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e2_b1",          &e2_b1);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_b1",          &e3_b1);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_v1_b1",       &e3_v1_b1);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_v2_b1",       &e3_v2_b1);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e2_b2",          &e2_b2);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_b2",          &e3_b2);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_v1_b2",       &e3_v1_b2);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_v1_b2",       &e3_v2_b2);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e2_sdb1",        &e2_sdb1);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_sdb1",        &e3_sdb1);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_v1_sdb1",     &e3_v1_sdb1);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_v2_sdb1",     &e3_v2_sdb1);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e2_sdb2",        &e2_sdb2);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_sdb2",        &e3_sdb2);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_v1_sdb2",     &e3_v1_sdb2);
    iTree->SetBranchAddress("AK8Puppijet"+jet+"_e3_v2_sdb2",     &e3_v2_sdb2);
  }
}
ZprimeBitsLoader::~ZprimeBitsLoader(){}
bool ZprimeBitsLoader::selectJetAlgoAndSize(TString algo){
  bool lPass = false;
  if((selectBits & kBOOSTED8PUPPI) && algo=="PUPPI") lPass = true;
  return lPass;
}
bool ZprimeBitsLoader::isPho(bool isData){
  bool lPass = false;
  if (nmu==0 && nele==0 && npho==1 && ntau==0){
   if(isData){
    if(triggerBits & kSinglePhoton)  lPass = true;}
   else{
    lPass = true;}
 }
 return lPass;
}

/*bool ZprimeBitsLoader::isHad(bool isData){
  bool lPass = false;
    if(nmu==0 && nele==0 && npho==0 && ntau==0 && met>175){
      if(isData){
        if((triggerBits & kMET)!=0)  lPass = true;
      }
      else{
       lPass = true;
      }
    }
   return lPass;
}*/

bool ZprimeBitsLoader::passBoostedZprimePreSelection(){
  //if((bst_jet0_msd>60) && (bst_jet0_msd<100)){ 
   return njets>0 & bst_jet0_pt>500;
   //  //else {return false;}
}

bool ZprimeBitsLoader::passPreSelection(bool isData, string preselection){
  bool lPass = false;
 // if(preselection.compare("Had")==0 && isHad(isData)) lPass = true;
  if(preselection.compare("Pho")==0 && isPho(isData)) lPass = true;
  return lPass;
}

bool ZprimeBitsLoader::passBoostedGammaZprimeSelection(){
  //if((bst_jet0_msd>60) && (bst_jet0_msd<100)){ 
 return njets>0 & bst_jet0_pt>175;
 //else {return false;}
}

bool ZprimeBitsLoader::passBoostedGammaZprimeSR(float ddtcut){
  
  return passBoostedGammaZprimeSelection() & (bst_jet0_tau21 < (-0.063*bst_jet0_rho + ddtcut));
}



bool ZprimeBitsLoader::passBoostedZprimeSR(float ddtcut){
  
  return passBoostedZprimePreSelection() & (bst_jet0_tau21 < (-0.063*bst_jet0_rho + ddtcut));
}


bool ZprimeBitsLoader::passBoostedZprimeBTag(float csvcut){
  return passBoostedZprimePreSelection() & (bst_jet0_doublecsv > csvcut);
}


bool ZprimeBitsLoader::passSelection(bool isData,string selection,string subsample, float ddt,float csv1){
  return passBoostedZprimePreSelection();
}


double ZprimeBitsLoader::getWgt(bool isData, TString algo, double LUMI){
  float wgt = 1;
  if(!isData) {     
    wgt *= LUMI*scale1fb*evtWeight;
    if (algo == "CHS") wgt *= puWeight;
  }
  return wgt;
}
double ZprimeBitsLoader::tau21DDT(){
  return bst_jet0_tau21 + 0.063*bst_jet0_rho;
}
double ZprimeBitsLoader::N2(double beta){
  if(beta == 1) return e3_v2_b1/(e2_b1*e2_b1);
  else return e3_v2_b2/(e2_b2*e2_b2);
}
double ZprimeBitsLoader::N2sd(double beta){
  if(beta == 1) return e3_v2_sdb1/(e2_sdb1*e2_sdb1);
  else return e3_v2_sdb2/(e2_sdb2*e2_sdb2);
}
double ZprimeBitsLoader::M2(double beta){
  if(beta == 1) return e3_v1_b1/e2_b1;
  else return e3_v1_b2/e2_b2;
}
double ZprimeBitsLoader::M2sd(double beta){
  if(beta == 1) return e3_v1_sdb1/e2_sdb1;
  else return e3_v1_sdb2/e2_sdb2;
}
double ZprimeBitsLoader::D2(double beta){
  if(beta == 1) return e3_b1/(e2_b1*e2_b1*e2_b1);
  else return e3_b2/(e2_b2*e2_b2*e2_b2);
}
double ZprimeBitsLoader::D2sd(double beta){
  if(beta == 1) return e3_sdb1/(e2_sdb1*e2_sdb1*e2_sdb1);
  else return e3_sdb2/(e2_sdb2*e2_sdb2*e2_sdb2);
}

