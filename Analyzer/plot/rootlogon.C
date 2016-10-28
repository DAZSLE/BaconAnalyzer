{
  TString path = gSystem->GetIncludePath();
  path += "-I. -I$ROOTSYS/src -I$ROOFITSYS/include";
  gSystem->SetIncludePath(path.Data());

  if(gSystem->Getenv("CMSSW_VERSION")) {
    gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
    gSystem->AddIncludePath("-I$CMSSW_RELEASE_BASE/src/");
    gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_BASE"))+"/src/");
    gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_RELEASE_BASE"))+"/src/");

  }

  // for plots
  gROOT->Macro("CPlot.cc+");
  gROOT->Macro("KStyle.cc+");
  gROOT->Macro("ZprimeBitsLoader.cc+");
}
