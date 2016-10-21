//============================================================================================
// Photon and Lepton Efficiencies functions:
// - Trigger Eff 
//============================================================================================

double getMuonTriggerDataEff(const double pt, const double eta, const unsigned int runNum=0);
double getMuonTriggerMCEff(const double pt, const double eta);
double getMuonTriggerSF(const double pt, const double eta);

double getEleTriggerDataEff(const double pt, const double eta);
double getEleTriggerMCEff(const double pt, const double eta);
double getEleTriggerSF(const double pt, const double eta);

double getEle27TriggerDataEff(const double pt, const double eta);
double getEle27TriggerMCEff(const double pt, const double eta);
double getEle27TriggerSF(const double pt, const double eta);

double getEle23TriggerDataEff(const double pt, const double eta);
double getEle23TriggerMCEff(const double pt, const double eta);
double getEle23TriggerSF(const double pt, const double eta);

//--------------------------------------------------------------------------------------
double getMuonTriggerDataEff(const double pt, const double eta, const unsigned int runNum)
{
  // Trigger efficiencies are computed for up to Run 257819 and for after Run 257819
  // In Run2015D, there is ~0.393/fb up to 257819

  if(runNum==0) {
    //
    // Take weighted average of data efficiencies of the two run eras weighted by relative 
    // integrated luminosity
    //
    double total_lumi = 2.246;
    double lumi1      = 0.393;
    double frac       = lumi1/total_lumi;

    double data_eff1 = getMuonTriggerDataEff(pt, eta, 1);
    double data_eff2 = getMuonTriggerDataEff(pt, eta, 999999);

    return (frac*data_eff1 + (1.-frac)*data_eff2);

  } else if(runNum <=257819) {
    if(fabs(eta) < 0.9) {
      if     (pt > 22 && pt < 25) { return 0.775; }
      else if(pt > 25 && pt < 30) { return 0.842; }
      else if(pt > 30 && pt < 40) { return 0.905; }
      else if(pt > 40 && pt < 50) { return 0.938; }
      else if(pt > 50 && pt < 60) { return 0.950; }
      else if(pt > 60 && pt <120) { return 0.944; }
      else if(pt >120)            { return 0.951; }

    } else if(fabs(eta) < 1.2) {
      if     (pt > 22 && pt < 25) { return 0.803; }
      else if(pt > 25 && pt < 30) { return 0.869; }
      else if(pt > 30 && pt < 40) { return 0.901; }
      else if(pt > 40 && pt < 50) { return 0.929; }
      else if(pt > 50 && pt < 60) { return 0.934; }
      else if(pt > 60 && pt <120) { return 0.933; }
      else if(pt >120)            { return 0.905; }

    } else if(fabs(eta) < 2.1) {
      if     (pt > 22 && pt < 25) { return 0.855; }
      else if(pt > 25 && pt < 30) { return 0.887; }
      else if(pt > 30 && pt < 40) { return 0.907; }
      else if(pt > 40 && pt < 50) { return 0.926; }
      else if(pt > 50 && pt < 60) { return 0.931; }
      else if(pt > 60 && pt <120) { return 0.932; }
      else if(pt >120)            { return 0.931; }

    } else if(fabs(eta) < 2.4) {
      if     (pt > 22 && pt < 25) { return 0.772; }
      else if(pt > 25 && pt < 30) { return 0.826; }
      else if(pt > 30 && pt < 40) { return 0.860; }
      else if(pt > 40 && pt < 50) { return 0.887; }
      else if(pt > 50 && pt < 60) { return 0.888; }
      else if(pt > 60 && pt <120) { return 0.891; }
      else if(pt >120)            { return 0.836; }
    }
  } else {
    if(fabs(eta) < 0.9) {
      if     (pt > 22 && pt < 25) { return 0.914; }
      else if(pt > 25 && pt < 30) { return 0.935; }
      else if(pt > 30 && pt < 40) { return 0.951; }
      else if(pt > 40 && pt < 50) { return 0.960; }
      else if(pt > 50 && pt < 60) { return 0.961; }
      else if(pt > 60 && pt <120) { return 0.964; }
      else if(pt >120)            { return 0.956; }

    } else if(fabs(eta) < 1.2) {
      if     (pt > 22 && pt < 25) { return 0.913; }
      else if(pt > 25 && pt < 30) { return 0.929; }
      else if(pt > 30 && pt < 40) { return 0.939; }
      else if(pt > 40 && pt < 50) { return 0.944; }
      else if(pt > 50 && pt < 60) { return 0.944; }
      else if(pt > 60 && pt <120) { return 0.941; }
      else if(pt >120)            { return 0.937; }

    } else if(fabs(eta) < 2.1) {
      if     (pt > 22 && pt < 25) { return 0.902; }
      else if(pt > 25 && pt < 30) { return 0.916; }
      else if(pt > 30 && pt < 40) { return 0.930; }
      else if(pt > 40 && pt < 50) { return 0.941; }
      else if(pt > 50 && pt < 60) { return 0.944; }
      else if(pt > 60 && pt <120) { return 0.944; }
      else if(pt >120)            { return 0.936; }

    } else if(fabs(eta) < 2.4) {
      if     (pt > 22 && pt < 25) { return 0.819; }
      else if(pt > 25 && pt < 30) { return 0.855; }
      else if(pt > 30 && pt < 40) { return 0.883; }
      else if(pt > 40 && pt < 50) { return 0.897; }
      else if(pt > 50 && pt < 60) { return 0.899; }
      else if(pt > 60 && pt <120) { return 0.898; }
      else if(pt >120)            { return 0.881; }
    }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getMuonTriggerMCEff(const double pt, const double eta)
{
  if(fabs(eta) < 0.9) {
    if     (pt > 22 && pt < 25) { return 0.891; }
    else if(pt > 25 && pt < 30) { return 0.916; }
    else if(pt > 30 && pt < 40) { return 0.943; }
    else if(pt > 40 && pt < 50) { return 0.960; }
    else if(pt > 50 && pt < 60) { return 0.964; }
    else if(pt > 60 && pt <120) { return 0.968; }
    else if(pt >120)            { return 0.957; }

  } else if(fabs(eta) < 1.2) {
    if     (pt > 22 && pt < 25) { return 0.880; }
    else if(pt > 25 && pt < 30) { return 0.897; }
    else if(pt > 30 && pt < 40) { return 0.929; }
    else if(pt > 40 && pt < 50) { return 0.945; }
    else if(pt > 50 && pt < 60) { return 0.948; }
    else if(pt > 60 && pt <120) { return 0.949; }
    else if(pt >120)            { return 0.939; }

  } else if(fabs(eta) < 2.1) {
    if     (pt > 22 && pt < 25) { return 0.887; }
    else if(pt > 25 && pt < 30) { return 0.911; }
    else if(pt > 30 && pt < 40) { return 0.932; }
    else if(pt > 40 && pt < 50) { return 0.950; }
    else if(pt > 50 && pt < 60) { return 0.952; }
    else if(pt > 60 && pt <120) { return 0.953; }
    else if(pt >120)            { return 0.940; }

  } else if(fabs(eta) < 2.4) {
    if     (pt > 22 && pt < 25) { return 0.811; }
    else if(pt > 25 && pt < 30) { return 0.831; }
    else if(pt > 30 && pt < 40) { return 0.863; }
    else if(pt > 40 && pt < 50) { return 0.889; }
    else if(pt > 50 && pt < 60) { return 0.886; }
    else if(pt > 60 && pt <120) { return 0.888; }
    else if(pt >120)            { return 0.890; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getMuonTriggerSF(const double pt, const double eta)
{
  double data_eff = getMuonTriggerDataEff(pt, eta, 0);
  double mc_eff   = getMuonTriggerMCEff(pt, eta);

  return data_eff / mc_eff;
}

//--------------------------------------------------------------------------------------
double getEleTriggerDataEff(const double pt, const double eta)
{
  if(fabs(eta) < 0.4) {
    if     (pt > 30 && pt < 32) { return 0.853; }
    else if(pt > 32 && pt < 35) { return 0.893; }
    else if(pt > 35 && pt < 40) { return 0.924; }
    else if(pt > 40 && pt < 50) { return 0.945; }
    else if(pt > 50 && pt <200) { return 0.967; }
    else                        { return 0.989; }

  } else if(fabs(eta) < 1.4) {
    if     (pt > 30 && pt < 32) { return 0.755; }
    else if(pt > 32 && pt < 35) { return 0.864; }
    else if(pt > 35 && pt < 40) { return 0.914; }
    else if(pt > 40 && pt < 50) { return 0.944; }
    else if(pt > 50 && pt <200) { return 0.963; }
    else                        { return 0.972; }

  } else if(fabs(eta) < 1.6) {
    if     (pt > 30 && pt < 32) { return 0.477; }
    else if(pt > 32 && pt < 35) { return 0.729; }
    else if(pt > 35 && pt < 40) { return 0.858; }
    else if(pt > 40 && pt < 50) { return 0.895; }
    else if(pt > 50 && pt <200) { return 0.917; }
    else                        { return 0.949; }

  } else if(fabs(eta) < 2.1) {
    if     (pt > 30 && pt < 32) { return 0.774; }
    else if(pt > 32 && pt < 35) { return 0.863; }
    else if(pt > 35 && pt < 40) { return 0.890; }
    else if(pt > 40 && pt < 50) { return 0.907; }
    else if(pt > 50 && pt <200) { return 0.908; }
    else                        { return 0.903; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getEleTriggerMCEff(const double pt, const double eta)
{
  if(fabs(eta) < 0.4) {
    if     (pt > 30 && pt < 32) { return 0.921; }
    else if(pt > 32 && pt < 35) { return 0.943; }
    else if(pt > 35 && pt < 40) { return 0.959; }
    else if(pt > 40 && pt < 50) { return 0.970; }
    else if(pt > 50 && pt <200) { return 0.974; }
    else                        { return 0.996; }

  } else if(fabs(eta) < 1.4) {
    if     (pt > 30 && pt < 32) { return 0.792; }
    else if(pt > 32 && pt < 35) { return 0.879; }
    else if(pt > 35 && pt < 40) { return 0.924; }
    else if(pt > 40 && pt < 50) { return 0.945; }
    else if(pt > 50 && pt <200) { return 0.955; }
    else                        { return 0.940; }

  } else if(fabs(eta) < 1.6) {
    if     (pt > 30 && pt < 32) { return 0.516; }
    else if(pt > 32 && pt < 35) { return 0.782; }
    else if(pt > 35 && pt < 40) { return 0.864; }
    else if(pt > 40 && pt < 50) { return 0.902; }
    else if(pt > 50 && pt <200) { return 0.918; }
    else                        { return 0.950; }

  } else if(fabs(eta) < 2.1) {
    if     (pt > 30 && pt < 32) { return 0.809; }
    else if(pt > 32 && pt < 35) { return 0.901; }
    else if(pt > 35 && pt < 40) { return 0.929; }
    else if(pt > 40 && pt < 50) { return 0.952; }
    else if(pt > 50 && pt <200) { return 0.964; }
    else                        { return 0.961; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getEleTriggerSF(const double pt, const double eta)
{
  return getEleTriggerDataEff(pt,eta) / getEleTriggerMCEff(pt,eta);
}

//--------------------------------------------------------------------------------------
double getEle27TriggerDataEff(const double pt, const double eta)
{
  if(fabs(eta) < 1.4442) {
    if     (pt > 10 && pt < 20) { return 0.000099; }
    else if(pt > 20 && pt < 25) { return 0.00089; }
    else if(pt > 25 && pt < 26) { return 0.0019; }
    else if(pt > 26 && pt < 26.5) { return 0.012; }
    else if(pt > 26.5 && pt < 26.75) { return 0.042; }
    else if(pt > 26.75 && pt < 27) { return 0.099; }
    else if(pt > 27 && pt < 27.25) { return 0.17; }
    else if(pt > 27.25 && pt < 27.5) { return 0.25; }
    else if(pt > 27.5 && pt < 28) { return 0.42; }
    else if(pt > 28 && pt < 28.5) { return 0.52; }
    else if(pt > 28.5 && pt < 29) { return 0.58; }
    else if(pt > 29 && pt < 30) { return 0.66; }
    else if(pt > 30 && pt < 31) { return 0.76; }
    else if(pt > 31 && pt < 32) { return 0.81; }
    else if(pt > 32 && pt < 33) { return 0.85; }
    else if(pt > 33 && pt < 35) { return 0.88; }
    else if(pt > 35 && pt < 40) { return 0.92; }
    else if(pt > 40 && pt < 45) { return 0.94; }
    else if(pt > 45 && pt < 50) { return 0.95; }
    else if(pt > 50 && pt < 60) { return 0.96; }
    else if(pt > 60 && pt < 70) { return 0.97; }
    else if(pt > 70 && pt < 85) { return 0.97; }
    else if(pt > 85 && pt < 100) { return 0.97; }
    else                        { return 0.99; }
  } else if(fabs(eta) > 1.566 && fabs(eta) < 2.5) {
    if     (pt > 10 && pt < 20) { return 0.00024; }
    else if(pt > 20 && pt < 25) { return 0.00086; }
    else if(pt > 25 && pt < 26) { return 0.0030; }
    else if(pt > 26 && pt < 26.5) { return 0.0058; }
    else if(pt > 26.5 && pt < 26.75) { return 0.023; }
    else if(pt > 26.75 && pt < 27) { return 0.021; }
    else if(pt > 27 && pt < 27.25) { return 0.019; }
    else if(pt > 27.25 && pt < 27.5) { return 0.045; }
    else if(pt > 27.5 && pt < 28) { return 0.082; }
    else if(pt > 28 && pt < 28.5) { return 0.20; }
    else if(pt > 28.5 && pt < 29) { return 0.36; }
    else if(pt > 29 && pt < 30) { return 0.56; }
    else if(pt > 30 && pt < 31) { return 0.76; }
    else if(pt > 31 && pt < 32) { return 0.82; }
    else if(pt > 32 && pt < 33) { return 0.85; }
    else if(pt > 33 && pt < 35) { return 0.87; }
    else if(pt > 35 && pt < 40) { return 0.89; }
    else if(pt > 40 && pt < 45) { return 0.90; }
    else if(pt > 45 && pt < 50) { return 0.91; }
    else if(pt > 50 && pt < 60) { return 0.91; }
    else if(pt > 60 && pt < 70) { return 0.91; }
    else if(pt > 70 && pt < 85) { return 0.91; }
    else if(pt > 85 && pt < 100) { return 0.89; }
    else                        { return 0.91; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getEle27TriggerMCEff(const double pt, const double eta)
{
  if(fabs(eta) < 1.4442) {
    if     (pt > 10 && pt < 20) { return 1; }
    else if(pt > 20 && pt < 25) { return 0.00020; }
    else if(pt > 25 && pt < 26) { return 0.0064; }
    else if(pt > 26 && pt < 26.5) { return 0.047; }
    else if(pt > 26.5 && pt < 26.75) { return 0.15; }
    else if(pt > 26.75 && pt < 27) { return 0.22; }
    else if(pt > 27 && pt < 27.25) { return 0.34; }
    else if(pt > 27.25 && pt < 27.5) { return 0.45; }
    else if(pt > 27.5 && pt < 28) { return 0.53; }
    else if(pt > 28 && pt < 28.5) { return 0.62; }
    else if(pt > 28.5 && pt < 29) { return 0.67; }
    else if(pt > 29 && pt < 30) { return 0.74; }
    else if(pt > 30 && pt < 31) { return 0.82; }
    else if(pt > 31 && pt < 32) { return 0.87; }
    else if(pt > 32 && pt < 33) { return 0.89; }
    else if(pt > 33 && pt < 35) { return 0.93; }
    else if(pt > 35 && pt < 40) { return 0.95; }
    else if(pt > 40 && pt < 45) { return 0.97; }
    else if(pt > 45 && pt < 50) { return 0.98; }
    else if(pt > 50 && pt < 60) { return 0.98; }
    else if(pt > 60 && pt < 70) { return 0.98; }
    else if(pt > 70 && pt < 85) { return 0.99; }
    else if(pt > 85 && pt < 100) { return 0.98; }
    else                        { return 0.99; }
  } else if(fabs(eta) > 1.566 && fabs(eta) < 2.5) {
    if     (pt > 10 && pt < 20) { return 0.00032; }
    else if(pt > 20 && pt < 25) { return 0.00093; }
    else if(pt > 25 && pt < 26) { return 0.0056; }
    else if(pt > 26 && pt < 26.5) { return 0.036; }
    else if(pt > 26.5 && pt < 26.75) { return 0.043; }
    else if(pt > 26.75 && pt < 27) { return 0.11; }
    else if(pt > 27 && pt < 27.25) { return 0.10; }
    else if(pt > 27.25 && pt < 27.5) { return 0.19; }
    else if(pt > 27.5 && pt < 28) { return 0.25; }
    else if(pt > 28 && pt < 28.5) { return 0.38; }
    else if(pt > 28.5 && pt < 29) { return 0.46; }
    else if(pt > 29 && pt < 30) { return 0.61; }
    else if(pt > 30 && pt < 31) { return 0.70; }
    else if(pt > 31 && pt < 32) { return 0.78; }
    else if(pt > 32 && pt < 33) { return 0.81; }
    else if(pt > 33 && pt < 35) { return 0.87; }
    else if(pt > 35 && pt < 40) { return 0.92; }
    else if(pt > 40 && pt < 45) { return 0.95; }
    else if(pt > 45 && pt < 50) { return 0.96; }
    else if(pt > 50 && pt < 60) { return 0.97; }
    else if(pt > 60 && pt < 70) { return 0.97; }
    else if(pt > 70 && pt < 85) { return 0.97; }
    else if(pt > 85 && pt < 100) { return 0.98; }
    else                        { return 0.98; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getEle27TriggerSF(const double pt, const double eta)
{
  return getEle27TriggerDataEff(pt,eta) / getEle27TriggerMCEff(pt,eta);
}

//--------------------------------------------------------------------------------------
double getEle23TriggerDataEff(const double pt, const double eta)
{
  if(fabs(eta) < 1.4442) {
    if     (pt > 10 && pt < 15) { return 0.0; }
    else if(pt > 15 && pt < 20) { return 0.00040; }
    else if(pt > 20 && pt < 22) { return 0.0018; }
    else if(pt > 22 && pt < 22.5) { return 0.012; }
    else if(pt > 22.5 && pt < 22.75) { return 0.044; }
    else if(pt > 22.75 && pt < 23) { return 0.088; }
    else if(pt > 23 && pt < 23.25) { return 0.17; }
    else if(pt > 23.25 && pt < 23.5) { return 0.27; }
    else if(pt > 23.5 && pt < 24) { return 0.39; }
    else if(pt > 24 && pt < 24.5) { return 0.49; }
    else if(pt > 24.5 && pt < 25) { return 0.60; }
    else if(pt > 25 && pt < 26) { return 0.66; }
    else if(pt > 26 && pt < 27) { return 0.73; }
    else if(pt > 27 && pt < 28) { return 0.79; }
    else if(pt > 28 && pt < 30) { return 0.84; }
    else if(pt > 30 && pt < 35) { return 0.89; }
    else if(pt > 35 && pt < 40) { return 0.92; }
    else if(pt > 40 && pt < 45) { return 0.94; }
    else if(pt > 45 && pt < 50) { return 0.95; }
    else if(pt > 50 && pt < 60) { return 0.96; }
    else if(pt > 60 && pt < 70) { return 0.97; }
    else if(pt > 70 && pt < 85) { return 0.97; }
    else if(pt > 85 && pt < 100) { return 0.97; }
    else                        { return 0.98; }
  } else if(fabs(eta) > 1.566 && fabs(eta) < 2.5) {
    if     (pt > 10 && pt < 15) { return 0.0; }
    else if(pt > 15 && pt < 20) { return 0.00033; }
    else if(pt > 20 && pt < 22) { return 0.00053; }
    else if(pt > 22 && pt < 22.5) { return 0.0085; }
    else if(pt > 22.5 && pt < 22.75) { return 0.0071; }
    else if(pt > 22.75 && pt < 23) { return 0.029; }
    else if(pt > 23 && pt < 23.25) { return 0.022; }
    else if(pt > 23.25 && pt < 23.5) { return 0.039; }
    else if(pt > 23.5 && pt < 24) { return 0.085; }
    else if(pt > 24 && pt < 24.5) { return 0.20; }
    else if(pt > 24.5 && pt < 25) { return 0.38; }
    else if(pt > 25 && pt < 26) { return 0.57; }
    else if(pt > 26 && pt < 27) { return 0.73; }
    else if(pt > 27 && pt < 28) { return 0.80; }
    else if(pt > 28 && pt < 30) { return 0.84; }
    else if(pt > 30 && pt < 35) { return 0.87; }
    else if(pt > 35 && pt < 40) { return 0.89; }
    else if(pt > 40 && pt < 45) { return 0.91; }
    else if(pt > 45 && pt < 50) { return 0.91; }
    else if(pt > 50 && pt < 60) { return 0.91; }
    else if(pt > 60 && pt < 70) { return 0.91; }
    else if(pt > 70 && pt < 85) { return 0.91; }
    else if(pt > 85 && pt < 100) { return 0.89; }
    else                        { return 0.91; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getEle23TriggerMCEff(const double pt, const double eta)
{
  if(fabs(eta) < 1.4442) {
    if     (pt > 10 && pt < 15) { return 0.00020; }
    else if(pt > 15 && pt < 20) { return 0.00022; }
    else if(pt > 20 && pt < 22) { return 0.0016; }
    else if(pt > 22 && pt < 22.5) { return 0.029; }
    else if(pt > 22.5 && pt < 22.75) { return 0.11; }
    else if(pt > 22.75 && pt < 23) { return 0.21; }
    else if(pt > 23 && pt < 23.25) { return 0.38; }
    else if(pt > 23.25 && pt < 23.5) { return 0.45; }
    else if(pt > 23.5 && pt < 24) { return 0.53; }
    else if(pt > 24 && pt < 24.5) { return 0.62; }
    else if(pt > 24.5 && pt < 25) { return 0.69; }
    else if(pt > 25 && pt < 26) { return 0.76; }
    else if(pt > 26 && pt < 27) { return 0.84; }
    else if(pt > 27 && pt < 28) { return 0.89; }
    else if(pt > 28 && pt < 30) { return 0.93; }
    else if(pt > 30 && pt < 35) { return 0.97; }
    else if(pt > 35 && pt < 40) { return 0.98; }
    else if(pt > 40 && pt < 45) { return 0.99; }
    else if(pt > 45 && pt < 50) { return 0.99; }
    else if(pt > 50 && pt < 60) { return 0.99; }
    else if(pt > 60 && pt < 70) { return 0.99; }
    else if(pt > 70 && pt < 85) { return 0.99; }
    else if(pt > 85 && pt < 100) { return 0.99; }
    else                        { return 0.99; }
  } else if(fabs(eta) > 1.566 && fabs(eta) < 2.5) {
    if     (pt > 10 && pt < 15) { return 0.0; }
    else if(pt > 15 && pt < 20) { return 0.00076; }
    else if(pt > 20 && pt < 22) { return 0.0042; }
    else if(pt > 22 && pt < 22.5) { return 0.035; }
    else if(pt > 22.5 && pt < 22.75) { return 0.044; }
    else if(pt > 22.75 && pt < 23) { return 0.13; }
    else if(pt > 23 && pt < 23.25) { return 0.16; }
    else if(pt > 23.25 && pt < 23.5) { return 0.17; }
    else if(pt > 23.5 && pt < 24) { return 0.33; }
    else if(pt > 24 && pt < 24.5) { return 0.46; }
    else if(pt > 24.5 && pt < 25) { return 0.56; }
    else if(pt > 25 && pt < 26) { return 0.69; }
    else if(pt > 26 && pt < 27) { return 0.80; }
    else if(pt > 27 && pt < 28) { return 0.85; }
    else if(pt > 28 && pt < 30) { return 0.91; }
    else if(pt > 30 && pt < 35) { return 0.95; }
    else if(pt > 35 && pt < 40) { return 0.97; }
    else if(pt > 40 && pt < 45) { return 0.98; }
    else if(pt > 45 && pt < 50) { return 0.98; }
    else if(pt > 50 && pt < 60) { return 0.99; }
    else if(pt > 60 && pt < 70) { return 0.99; }
    else if(pt > 70 && pt < 85) { return 0.99; }
    else if(pt > 85 && pt < 100) { return 0.99; }
    else                        { return 0.99; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getEle23TriggerSF(const double pt, const double eta)
{
  return getEle23TriggerDataEff(pt,eta) / getEle23TriggerMCEff(pt,eta);
}
