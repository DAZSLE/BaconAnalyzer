#include "../include/BTagWeightLoader.hh"
#include <cmath>
#include <iostream> 
#include <string>
#include <sstream>

using namespace baconhep;

BTagWeightLoader::BTagWeightLoader(TTree *iTree,std::string btagScaleFactorFilename) {
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = btagScaleFactorFilename.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
    btagScaleFactorFilename.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }
  fJetCalib = new BTagCalibration("csvv2",btagScaleFactorFilename);
  freadersL.clear(); freadersM.clear(); freadersT.clear();
  freaders.clear();
  for(auto imtype : measurementTypes) { // freadersL 6 , freadersM 6, freadersT 6
    for(auto ivtype : variationTypes) {
      freadersL.push_back(new BTagCalibrationReader(fJetCalib, BTagEntry::OP_LOOSE,  imtype, ivtype)); // first mujets(HF) then incl(LF) and first central(0,3) then up(1,4) and then down(2,5)
      freadersM.push_back(new BTagCalibrationReader(fJetCalib, BTagEntry::OP_MEDIUM, imtype, ivtype));
      freadersT.push_back(new BTagCalibrationReader(fJetCalib, BTagEntry::OP_TIGHT,  imtype, ivtype));
    }
  }
  freaders.push_back(freadersL); freaders.push_back(freadersM); freaders.push_back(freadersT);
}
BTagWeightLoader::~BTagWeightLoader() {
  delete fJetCalib;
}
void BTagWeightLoader::reset() {
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 1;
}
void BTagWeightLoader::setupTree(TTree *iTree, std::string iJetLabel) {
  reset();
  for(int i0 = 0; i0 < 60; i0++) {float pVar = 1; fVars.push_back(pVar);} 
  int i1 = 0;
  for(auto iwptype : wpTypes) {
    addBTag(iJetLabel.c_str(),iTree,iwptype,fLabels,i1,fVars);
    i1 += 20;
  }
}
void BTagWeightLoader::addBTag(std::string iHeader,TTree *iTree,std::string iLabel,std::vector<std::string> iLabels,int iN,std::vector<float> &iVals) {
  int iBase=iN;
  for(int i0 = 0; i0 < int(iLabels.size()); i0++) {
    std::stringstream pVal0,pVal1,pValminus1,pVal2;
    pVal0       << iHeader << "btagw" << iLabel << "0"      << "_" << iLabels[i0 % iLabels.size()]; //res_PUPPIbtagwL0_CENT -- where iLabel(L,M o T) and iLabels(CENT,MISTAGUP,MISTAGD0,BTAGUP,BTAGD0)
    pVal1       << iHeader << "btagw" << iLabel << "1"      << "_" << iLabels[i0 % iLabels.size()]; //res_PUPPIbtagwL1_CENT
    pValminus1  << iHeader << "btagw" << iLabel << "minus1" << "_" << iLabels[i0 % iLabels.size()]; //res_PUPPIbtagwLminus1_CENT
    pVal2       << iHeader << "btagw" << iLabel << "2"      << "_" << iLabels[i0 % iLabels.size()]; //res_PUPPIbtagwL2_CENT
    iTree->Branch(pVal0      .str().c_str(),&iVals[iBase+0],(pVal0      .str()+"/F").c_str());
    iTree->Branch(pVal1      .str().c_str(),&iVals[iBase+1],(pVal1      .str()+"/F").c_str());
    //iTree->Branch(pValminus1 .str().c_str(),&iVals[iBase+2],(pValminus1 .str()+"/F").c_str());
    //iTree->Branch(pVal2      .str().c_str(),&iVals[iBase+3],(pVal2      .str()+"/F").c_str());
    iBase+=4;
  }
}
void BTagWeightLoader::fillBTag(std::vector<const TJet*> iObjects) {
  // vSFL should contain CENT (), MISTAG(Ms), BTAG(Bs)  - 5 - CENT(vSFL.at(0)),MsUP(vSFL.at(1)),MsDO(vSFL.at(2)),BsUP(vSFL.at(3)),BsDO(vSFL.at(4))
  int iN = 0;
  for(unsigned int j0=0; j0<3; j0++){  // L, M, T
    std::vector<std::vector<float>> vSFL,vSFL_nominal;
    vSFL.clear(); vSFL_nominal.clear();
    vSFL_nominal.push_back(getJetSFs("nominal",iObjects, freaders[j0].at(0), freaders[j0].at(3)));
    for(auto iftype :flavorTypes) {    
      vSFL_nominal.push_back(getJetSFs(iftype,iObjects, freaders[j0].at(0), freaders[j0].at(3))); // 0 and 3 HF and LF respectively - flavor types: Ms,Bs
    }
    vSFL.push_back(vSFL_nominal.at(0));

    for(auto iftype :flavorTypes) {
      for(unsigned int i0=1; i0<3; i0++){
	std::vector<float> vSF0;  vSF0.clear();
	for(unsigned int i1=0; i1<(vSFL.at(0)).size(); i1++) {
	  if(iftype.compare("Ms")==0) vSF0.push_back( (getJetSFs(iftype,iObjects, freaders[j0].at(i0), freaders[j0].at(i0+3))).at(i1) * (vSFL_nominal.at(2).at(i1)));
	  if(iftype.compare("Bs")==0) vSF0.push_back( (getJetSFs(iftype,iObjects, freaders[j0].at(i0), freaders[j0].at(i0+3))).at(i1) * (vSFL_nominal.at(1).at(i1)));
	}
	vSFL.push_back(vSF0);
      }
    }
    
    // Fill btag
    for(unsigned int j1=0; j1<5; j1++){
      int lBase = j1*4+iN;
      fVars[lBase+0] = getBtagEventReweight(0,  iObjects, vSFL.at(j1));
      fVars[lBase+1] = getBtagEventReweight(1,  iObjects, vSFL.at(j1));
      fVars[lBase+2] = getBtagEventReweight(-1, iObjects, vSFL.at(j1));
      fVars[lBase+3] = getBtagEventReweight(2,  iObjects, vSFL.at(j1));
    }
    iN += 20;
  }
}
