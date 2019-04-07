#include "../include/JECLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>

using namespace baconhep;

JECLoader::JECLoader(bool iData,std::string iLabel,std::string iJet) {
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";

  loadJECs(iData,iLabel,iJet);
}
JECLoader::~JECLoader() { 
}
void JECLoader::reset() { 
  isData = false;
  labelEra = "";
  jetType = "";
  correctionParameters.clear();
  JetCorrector.clear();
  jecUnc.clear();
  JetCorrectionsIOV.clear();
}
void JECLoader::loadJECs(bool iData,std::string iLabel,std::string iJet) {
  reset();
  isData = iData;
  labelEra = iLabel;
  jetType = iJet;
  std::cout << "JECLoader: loading " << labelEra << " jet energy correction constants for "<< jetType << std::endl;
  // initialize
  loadCMSSWPath();
  std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
  correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
  JetCorrector = std::vector<FactorizedJetCorrector*>();
  jecUnc = std::vector<JetCorrectionUncertainty*>();
  JetCorrectionsIOV = std::vector<std::pair<int,int> >();

  std::string jerpath = "";
  std::string jerrespath = "";
  std::string jersfpath = "";
  if(labelEra==label2016){
    jerpath = jer2016;
    jerrespath = jerres2016;
    jersfpath = jersf2016;
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      std::stringstream jec2016BCD; jec2016BCD << jecPathname << "/" << jetRecalib2016BCD << "/" << jetRecalib2016BCD;
      correctionParametersBCD.push_back(JetCorrectorParameters(Form("%s_L1FastJet_%s.txt", jec2016BCD.str().c_str(), jetType.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(Form("%s_L2Relative_%s.txt", jec2016BCD.str().c_str(), jetType.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(Form("%s_L3Absolute_%s.txt", jec2016BCD.str().c_str(), jetType.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(Form("%s_L2L3Residual_%s.txt", jec2016BCD.str().c_str(), jetType.c_str())));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(Form("%s_Uncertainty_%s.txt", jec2016BCD.str().c_str(), jetType.c_str()));
      correctionParameters.push_back(correctionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));
      //IOV: 2016EF
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      std::stringstream jec2016EF; jec2016EF << jecPathname << "/" << jetRecalib2016EF << "/" << jetRecalib2016EF;
      correctionParametersEF.push_back(JetCorrectorParameters(Form("%s_L1FastJet_%s.txt", jec2016EF.str().c_str(), jetType.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(Form("%s_L2Relative_%s.txt", jec2016EF.str().c_str(), jetType.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(Form("%s_L3Absolute_%s.txt", jec2016EF.str().c_str(), jetType.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(Form("%s_L2L3Residual_%s.txt", jec2016EF.str().c_str(), jetType.c_str())));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(Form("%s_Uncertainty_%s.txt", jec2016EF.str().c_str(), jetType.c_str()));
      correctionParameters.push_back(correctionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278801 ));
      //IOV: 2016GH
      std::vector<JetCorrectorParameters> correctionParametersGH = std::vector<JetCorrectorParameters> ();
      std::stringstream jec2016GH; jec2016GH << jecPathname << "/" << jetRecalib2016GH << "/" << jetRecalib2016GH;
      correctionParametersGH.push_back(JetCorrectorParameters(Form("%s_L1FastJet_%s.txt", jec2016GH.str().c_str(), jetType.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(Form("%s_L2Relative_%s.txt", jec2016GH.str().c_str(), jetType.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(Form("%s_L3Absolute_%s.txt", jec2016GH.str().c_str(), jetType.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(Form("%s_L2L3Residual_%s.txt", jec2016GH.str().c_str(), jetType.c_str())));
      FactorizedJetCorrector *JetCorrectorGH = new FactorizedJetCorrector(correctionParametersGH);
      JetCorrectionUncertainty *jecUncGH = new JetCorrectionUncertainty(Form("%s_Uncertainty_%s.txt", jec2016GH.str().c_str(), jetType.c_str()));
      correctionParameters.push_back(correctionParametersGH);
      JetCorrector.push_back( JetCorrectorGH );
      jecUnc.push_back(jecUncGH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 99999999 ));
    }
  }
  else if(labelEra==label2017){
    jerpath = jer2017;
    jerrespath = jerres2017;
    jersfpath =jersf2017;
    if (isData) {
      //IOV: 2017B
      std::vector<JetCorrectorParameters> correctionParametersB = std::vector<JetCorrectorParameters> ();
      std::stringstream jec2017B; jec2017B << jecPathname << "/" << jetRecalib2017B << "/" << jetRecalib2017B;
      correctionParametersB.push_back(JetCorrectorParameters(Form("%s_L1FastJet_%s.txt", jec2017B.str().c_str(), jetType.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(Form("%s_L2Relative_%s.txt", jec2017B.str().c_str(), jetType.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(Form("%s_L3Absolute_%s.txt", jec2017B.str().c_str(), jetType.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(Form("%s_L2L3Residual_%s.txt", jec2017B.str().c_str(), jetType.c_str())));
      FactorizedJetCorrector *JetCorrectorB = new FactorizedJetCorrector(correctionParametersB);
      JetCorrectionUncertainty *jecUncB = new JetCorrectionUncertainty(Form("%s_Uncertainty_%s.txt", jec2017B.str().c_str(), jetType.c_str()));
      correctionParameters.push_back(correctionParametersB);
      JetCorrector.push_back( JetCorrectorB );
      jecUnc.push_back(jecUncB);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 299329 ));
      //IOV: 2017C
      std::vector<JetCorrectorParameters> correctionParametersC = std::vector<JetCorrectorParameters> ();
      std::stringstream jec2017C; jec2017C << jecPathname << "/" << jetRecalib2017C << "/" << jetRecalib2017C;
      correctionParametersC.push_back(JetCorrectorParameters(Form("%s_L1FastJet_%s.txt", jec2017C.str().c_str(), jetType.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(Form("%s_L2Relative_%s.txt", jec2017C.str().c_str(), jetType.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(Form("%s_L3Absolute_%s.txt", jec2017C.str().c_str(), jetType.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(Form("%s_L2L3Residual_%s.txt", jec2017C.str().c_str(), jetType.c_str())));
      FactorizedJetCorrector *JetCorrectorC = new FactorizedJetCorrector(correctionParametersC);
      JetCorrectionUncertainty *jecUncC = new JetCorrectionUncertainty(Form("%s_Uncertainty_%s.txt", jec2017C.str().c_str(), jetType.c_str()));
      correctionParameters.push_back(correctionParametersC);
      JetCorrector.push_back( JetCorrectorC );
      jecUnc.push_back(jecUncC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 299368, 302029 ));
      //IOV: 2017DE
      std::vector<JetCorrectorParameters> correctionParametersDE = std::vector<JetCorrectorParameters> ();
      std::stringstream jec2017DE; jec2017DE << jecPathname << "/" << jetRecalib2017DE << "/" << jetRecalib2017DE;
      correctionParametersDE.push_back(JetCorrectorParameters(Form("%s_L1FastJet_%s.txt", jec2017DE.str().c_str(), jetType.c_str())));
      correctionParametersDE.push_back(JetCorrectorParameters(Form("%s_L2Relative_%s.txt", jec2017DE.str().c_str(), jetType.c_str())));
      correctionParametersDE.push_back(JetCorrectorParameters(Form("%s_L3Absolute_%s.txt", jec2017DE.str().c_str(), jetType.c_str())));
      correctionParametersDE.push_back(JetCorrectorParameters(Form("%s_L2L3Residual_%s.txt", jec2017DE.str().c_str(), jetType.c_str())));
      FactorizedJetCorrector *JetCorrectorDE = new FactorizedJetCorrector(correctionParametersDE);
      JetCorrectionUncertainty *jecUncDE = new JetCorrectionUncertainty(Form("%s_Uncertainty_%s.txt", jec2017DE.str().c_str(), jetType.c_str()));
      correctionParameters.push_back(correctionParametersDE);
      JetCorrector.push_back( JetCorrectorDE );
      jecUnc.push_back(jecUncDE);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 302031, 304797));
      //IOV: 2017F
      std::vector<JetCorrectorParameters> correctionParametersF = std::vector<JetCorrectorParameters> ();
      std::stringstream jec2017F; jec2017F << jecPathname << "/" << jetRecalib2017F << "/" << jetRecalib2017F;
      correctionParametersF.push_back(JetCorrectorParameters(Form("%s_L1FastJet_%s.txt", jec2017F.str().c_str(), jetType.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(Form("%s_L2Relative_%s.txt", jec2017F.str().c_str(), jetType.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(Form("%s_L3Absolute_%s.txt", jec2017F.str().c_str(), jetType.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(Form("%s_L2L3Residual_%s.txt", jec2017F.str().c_str(), jetType.c_str())));
      FactorizedJetCorrector *JetCorrectorF = new FactorizedJetCorrector(correctionParametersF);
      JetCorrectionUncertainty *jecUncF = new JetCorrectionUncertainty(Form("%s_Uncertainty_%s.txt", jec2017F.str().c_str(), jetType.c_str()));
      correctionParameters.push_back(correctionParametersF);
      JetCorrector.push_back( JetCorrectorF );
      jecUnc.push_back(jecUncF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 305044, 99999999 ));
    }
    std::stringstream jecMC; jecMC << jecPathname << "/" << jetRecalib2017MC <<"/" << jetRecalib2017MC;
  }
  else if(labelEra==label2018){
    jerpath = jer2018;
    jerrespath = jerres2018;
    jersfpath =jersf2018;
    std::stringstream jecMC; jecMC << jecPathname << "/" << jetRecalib2018MC <<"/" << jetRecalib2017MC;
  }
  else{
    jerpath = jer2017;
    jerrespath = jerres2017;
    jersfpath =jersf2017;
    std::stringstream jecMC; jecMC << jecPathname << "/" << jetRecalib2017MC <<"/" << jetRecalib2017MC;
  }

  std::stringstream respath; respath << jecPathname << "/" << jerpath << "/" << jerrespath << jetType << ".txt";
  std::stringstream ressfpath; ressfpath << jecPathname << "/" << jerpath << "/" <<jerrespath << jetType << ".txt";
  resolution = JME::JetResolution(respath.str().c_str());
  resolution_sf = JME::JetResolutionScaleFactor(ressfpath.str().c_str());
  if (!isData) {
    std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
    std::stringstream jecMC; jecMC << jecPathname << "/" << jetRecalib2017MC << "/" << jetRecalib2017MC;
    correctionParametersMC.push_back(JetCorrectorParameters(Form("%s_L1MCastJet_%s.txt", jecMC.str().c_str(), jetType.c_str())));
    correctionParametersMC.push_back(JetCorrectorParameters(Form("%s_L2Relative_%s.txt", jecMC.str().c_str(), jetType.c_str())));
    correctionParametersMC.push_back(JetCorrectorParameters(Form("%s_L3Absolute_%s.txt", jecMC.str().c_str(), jetType.c_str())));
    correctionParametersMC.push_back(JetCorrectorParameters(Form("%s_L2L3Residual_%s.txt", jecMC.str().c_str(), jetType.c_str())));
    FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
    JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(Form("%s_Uncertainty_%s.txt", jecMC.str().c_str(), jetType.c_str()));
    JetCorrector.push_back( JetCorrectorMC );
    jecUnc.push_back(jecUncMC);
    JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
  }

}
void JECLoader::loadCMSSWPath() {
  char* cmsswPathChar = getenv("CMSSW_BASE");
  if (cmsswPathChar == NULL) {
    std::cout << "Warning in JECLoader::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
    cmsswPath = "";
  }
  cmsswPath = std::string(cmsswPathChar);
}

// Retrieve jet energy uncertainty as a function of pt and eta
double JECLoader::getJecUnc( float pt, float eta , int run) {

  int foundIndex = -1;
  for (uint i=0; i<JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    std::cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
    foundIndex = 0;
  }

  jecUnc[foundIndex]->setJetPt(pt);
  jecUnc[foundIndex]->setJetEta(eta);
  return jecUnc[foundIndex]->getUncertainty(true);
}

//Jet Energy Corrections
double JECLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
					     double rho, double jetArea,
					     int run,
					     bool printDebug,
					     int jetCorrectionLevel){
  int foundIndex = -1;
  for (uint i=0; i<JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    std::cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
    foundIndex = 0;
  }
  
  if (!JetCorrector[foundIndex]) {
    std::cout << "WWARNING: Jet corrector pointer is null. Returning JEC = 0. \n";
    return 0;
  }
  
  JetCorrector[foundIndex]->setJetEta(jetEta);
  JetCorrector[foundIndex]->setJetPt(jetRawPt);
  JetCorrector[foundIndex]->setJetPhi(jetPhi);
  JetCorrector[foundIndex]->setJetE(jetE);
  JetCorrector[foundIndex]->setRho(rho);
  JetCorrector[foundIndex]->setJetA(jetArea);
  
  std::vector<float> corrections;
  corrections = JetCorrector[foundIndex]->getSubCorrections();

  if (printDebug) std::cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jetRawPt << " " << jetEta << " " << jetPhi << "\n";

  double cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {

    //only correct up to the required level. if -1, then do all correction levels
    if (jetCorrectionLevel >= 0 && int(j) > jetCorrectionLevel) continue;

    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) std::cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) std::cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";
  
  return cumulativeCorrection;

}

