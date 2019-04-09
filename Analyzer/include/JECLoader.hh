#ifndef JECLoader_H
#define JECLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Utils.hh"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TRandom3.h"

using namespace baconhep;

class JECLoader { 
public:
  JECLoader(bool iData,std::string iLabel,std::string iJet);
  ~JECLoader();
  void reset();
  void loadJECs(bool isData,std::string iLabel,std::string iJet="AK8PFPuppi");
  std::vector<FactorizedJetCorrector*> getJetCorrector() { return JetCorrector; }
  std::vector<std::pair<int,int> > getJetCorrectionsIOV() { return JetCorrectionsIOV; }
  double getJecUnc( float pt, float eta, int run );
  double JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
				    double rho, double jetArea,
				    int run,
                                    bool printDebug = false,
				    int jetCorrectionLevel = -1);

  std::vector<std::vector<JetCorrectorParameters> > correctionParameters;
  std::vector<FactorizedJetCorrector*> JetCorrector;
  std::vector<JetCorrectionUncertainty*> jecUnc;
  std::vector<std::pair<int,int> > JetCorrectionsIOV;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor resolution_sf;

protected: 
  std::string cmsswPath, labelEra, jetType;
  bool isData;
  void loadCMSSWPath();

  // JECs
  // 2016: 07Aug2017
  std::string jetRecalib2016BCD = "Summer16_07Aug2017BCD_V11_DATA";
  std::string jetRecalib2016EF = "Summer16_07Aug2017EF_V11_DATA";
  std::string jetRecalib2016GH = "Summer16_07Aug2017GH_V11_DATA";
  std::string jetRecalib2016MC = "Summer16_07Aug2017_V11_MC";

  // 2017: 17Nov2017 
  std::string jetRecalib2017B = "Fall17_17Nov2017B_V32_DATA";
  std::string jetRecalib2017C = "Fall17_17Nov2017C_V32_DATA";
  std::string jetRecalib2017DE = "Fall17_17Nov2017DE_V32_DATA";
  std::string jetRecalib2017F = "Fall17_17Nov2017F_V32_DATA";
  std::string jetRecalib2017MC = "Fall17_17Nov2017_V32_MC";

  // 2018: (Prompt?) or whatever is available now
  std::string jetRecalib2018A = "Autumn18_RunA_V8_DATA";
  std::string jetRecalib2018B = "Autumn18_RunB_V8_DATA";
  std::string jetRecalib2018C = "Autumn18_RunC_V8_DATA";
  std::string jetRecalib2018D = "Autumn18_RunD_V8_DATA";
  std::string jetRecalib2018MC = "Autumn18_V8_MC";

  // JERs
  std::string jer2016 = "Summer16_25nsV1_MC";
  std::string jerres2016 = "Summer16_25nsV1_MC_PtResolution_";
  std::string jersf2016 = "Summer16_25nsV1_MC_SF_";

  std::string jer2017 = "Fall17_V3_MC";
  std::string jerres2017 = "Fall17_V3_MC_PtResolution_";
  std::string jersf2017 = "Fall17_V3_MC_SF_";

  std::string jer2018 = "Fall17_V3_MC"; // use this until post-Moriond JERs are out
  std::string jerres2018 = "Fall17_V3_MC_PtResolution_";
  std::string jersf2018 = "Fall17_V3_MC_SF_";
};
#endif
