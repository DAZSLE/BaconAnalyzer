//============================================================================================
// Photon and Lepton Efficiencies functions:
// - ID Eff for loose and tight muon and ele
// - Trigger Eff 
//============================================================================================

double getLooseMuonSF(const double pt, const double eta);
double getTightMuonSF(const double pt, const double eta);

double getLooseMuonDataIDEff(const double pt, const double eta);
double getLooseMuonMCIDEff(const double pt, const double eta);
double getTightMuonDataIDEff(const double pt, const double eta);
double getTightMuonMCIDEff(const double pt, const double eta);

double getVetoEleSF(const double pt, const double eta);
double getTightEleSF(const double pt, const double eta);

double getVetoEleDataIDEff(const double pt, const double eta);
double getVetoEleMCIDEff(const double pt, const double eta);
double getTightEleDataIDEff(const double pt, const double eta);
double getTightEleMCIDEff(const double pt, const double eta);

double getMediumPhotonDataIDEff(const double pt, const double eta);
double getMediumPhotonMCIDEff(const double pt, const double eta);

double getMediumPhotonSF(const double pt, const double eta);

//--------------------------------------------------------------------------------------
double getLooseMuonSF(const double pt, const double eta)
{
  if(fabs(eta)<0.9) {
    if     (pt >  10 && pt <  15) { return 1.004; }
    else if(pt >  15 && pt <  20) { return 0.977; }
    else if(pt >  20 && pt <  25) { return 0.987; }
    else if(pt>   25 && pt <  30) { return 0.997; }
    else if(pt >  30 && pt <  40) { return 0.999; }
    else if(pt >  40 && pt <  50) { return 1.000; }
    else if(pt >  50 && pt <  60) { return 1.000; }
    else if(pt >  60 && pt <  80) { return 0.999; }
    else if(pt >  80 && pt < 100) { return 0.995; }
    else if(pt > 100)             { return 0.998; }

  } else if(fabs(eta)<1.2) {
    if     (pt >  10 && pt <  15) { return 0.968; }
    else if(pt >  15 && pt <  20) { return 0.976; }
    else if(pt >  20 && pt <  25) { return 0.989; }
    else if(pt>   25 && pt <  30) { return 1.002; }
    else if(pt >  30 && pt <  40) { return 1.000; }
    else if(pt >  40 && pt <  50) { return 1.000; }
    else if(pt >  50 && pt <  60) { return 1.001; }
    else if(pt >  60 && pt <  80) { return 1.004; }
    else if(pt >  80 && pt < 100) { return 0.993; }
    else if(pt > 100)             { return 1.000; }

  } else if(fabs(eta)<2.1) {
    if     (pt >  10 && pt <  15) { return 0.993; }
    else if(pt >  15 && pt <  20) { return 0.991; }
    else if(pt >  20 && pt <  25) { return 0.997; }
    else if(pt>   25 && pt <  30) { return 1.001; }
    else if(pt >  30 && pt <  40) { return 1.001; }
    else if(pt >  40 && pt <  50) { return 1.000; }
    else if(pt >  50 && pt <  60) { return 0.999; }
    else if(pt >  60 && pt <  80) { return 0.997; }
    else if(pt >  80 && pt < 100) { return 1.007; }
    else if(pt > 100)             { return 0.998; }

  } else if(fabs(eta)<2.4) {
    if     (pt >  10 && pt <  15) { return 0.961; }
    else if(pt >  15 && pt <  20) { return 0.986; }
    else if(pt >  20 && pt <  25) { return 0.982; }
    else if(pt>   25 && pt <  30) { return 1.004; }
    else if(pt >  30 && pt <  40) { return 0.998; }
    else if(pt >  40 && pt <  50) { return 1.000; }
    else if(pt >  50 && pt <  60) { return 0.999; }
    else if(pt >  60 && pt <  80) { return 1.002; }
    else if(pt >  80 && pt < 100) { return 1.001; }
    else if(pt > 100)             { return 1.007; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getTightMuonSF(const double pt, const double eta)
{
  if(fabs(eta)<0.9) {
    if     (pt >  10 && pt <  15) { return 1.004; }
    else if(pt >  15 && pt <  20) { return 0.960; }
    else if(pt >  20 && pt <  25) { return 0.977; }
    else if(pt>   25 && pt <  30) { return 0.985; }
    else if(pt >  30 && pt <  40) { return 0.988; }
    else if(pt >  40 && pt <  50) { return 0.991; }
    else if(pt >  50 && pt <  60) { return 0.987; }
    else if(pt >  60 && pt <  80) { return 0.987; }
    else if(pt >  80 && pt < 100) { return 0.982; }
    else if(pt > 100)             { return 0.981; }

  } else if(fabs(eta)<1.2) {
    if     (pt >  10 && pt <  15) { return 0.974; }
    else if(pt >  15 && pt <  20) { return 0.984; }
    else if(pt >  20 && pt <  25) { return 0.982; }
    else if(pt>   25 && pt <  30) { return 0.984; }
    else if(pt >  30 && pt <  40) { return 0.987; }
    else if(pt >  40 && pt <  50) { return 0.986; }
    else if(pt >  50 && pt <  60) { return 0.985; }
    else if(pt >  60 && pt <  80) { return 0.985; }
    else if(pt >  80 && pt < 100) { return 0.977; }
    else if(pt > 100)             { return 0.971; }

  } else if(fabs(eta)<2.1) {
    if     (pt >  10 && pt <  15) { return 0.962; }
    else if(pt >  15 && pt <  20) { return 0.985; }
    else if(pt >  20 && pt <  25) { return 0.982; }
    else if(pt>   25 && pt <  30) { return 0.988; }
    else if(pt >  30 && pt <  40) { return 0.994; }
    else if(pt >  40 && pt <  50) { return 0.993; }
    else if(pt >  50 && pt <  60) { return 0.993; }
    else if(pt >  60 && pt <  80) { return 0.994; }
    else if(pt >  80 && pt < 100) { return 0.988; }
    else if(pt > 100)             { return 0.988; }

  } else if(fabs(eta)<2.4) {
    if     (pt >  10 && pt <  15) { return 0.908; }
    else if(pt >  15 && pt <  20) { return 0.946; }
    else if(pt >  20 && pt <  25) { return 0.957; }
    else if(pt>   25 && pt <  30) { return 0.972; }
    else if(pt >  30 && pt <  40) { return 0.976; }
    else if(pt >  40 && pt <  50) { return 0.981; }
    else if(pt >  50 && pt <  60) { return 0.978; }
    else if(pt >  60 && pt <  80) { return 0.975; }
    else if(pt >  80 && pt < 100) { return 0.982; }
    else if(pt > 100)             { return 0.973; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getLooseMuonDataIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.2) {
    if     (pt >  10 && pt <  20) { return 0.72; }
    else if(pt >  20 && pt <  30) { return 0.85; }
    else if(pt >  30 && pt <  40) { return 0.95; }
    else if(pt>   40 && pt <  50) { return 0.98; }
    else if(pt >  50 && pt <  70) { return 0.98; }
    else if(pt >  70 && pt <  100) { return 0.98; }
    else if(pt > 100)             { return 0.98; }
  } else if(fabs(eta)>1.2 && fabs(eta)<2.5) {
    if     (pt >  10 && pt <  20) { return 0.78; }
    else if(pt >  20 && pt <  30) { return 0.88; }
    else if(pt >  30 && pt <  40) { return 0.96; }
    else if(pt>   40 && pt <  50) { return 0.99; }
    else if(pt >  50 && pt <  70) { return 0.99; }
    else if(pt >  70 && pt <  100) { return 0.99; }
    else if(pt > 100)             { return 0.99; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getLooseMuonMCIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.2) {
    if     (pt >  10 && pt <  20) { return 0.73; }
    else if(pt >  20 && pt <  30) { return 0.86; }
    else if(pt >  30 && pt <  40) { return 0.95; }
    else if(pt>   40 && pt <  50) { return 0.99; }
    else if(pt >  50 && pt <  70) { return 0.99; }
    else if(pt >  70 && pt <  100) { return 0.99; }
    else if(pt > 100)             { return 0.99; }
  } else if(fabs(eta)>1.2 && fabs(eta)<2.5) {
    if     (pt >  10 && pt <  20) { return 0.80; }
    else if(pt >  20 && pt <  30) { return 0.88; }
    else if(pt >  30 && pt <  40) { return 0.96; }
    else if(pt>   40 && pt <  50) { return 0.99; }
    else if(pt >  50 && pt <  70) { return 0.99; }
    else if(pt >  70 && pt <  100) { return 0.99; }
    else if(pt > 100)             { return 0.99; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getTightMuonDataIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.2) {
    if     (pt >  10 && pt <  20) { return 0.58; }
    else if(pt >  20 && pt <  30) { return 0.73; }
    else if(pt >  30 && pt <  40) { return 0.86; }
    else if(pt>   40 && pt <  50) { return 0.92; }
    else if(pt >  50 && pt <  70) { return 0.93; }
    else if(pt >  70 && pt <  100) { return 0.94; }
    else if(pt > 100)             { return 0.93; }
  } else if(fabs(eta)>1.2 && fabs(eta)<2.5) {
    if     (pt >  10 && pt <  20) { return 0.65; }
    else if(pt >  20 && pt <  30) { return 0.78; }
    else if(pt >  30 && pt <  40) { return 0.88; }
    else if(pt>   40 && pt <  50) { return 0.94; }
    else if(pt >  50 && pt <  70) { return 0.95; }
    else if(pt >  70 && pt <  100) { return 0.94; }
    else if(pt > 100)             { return 0.94; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getTightMuonMCIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.2) {
    if     (pt >  10 && pt <  20) { return 0.60; }
    else if(pt >  20 && pt <  30) { return 0.74; }
    else if(pt >  30 && pt <  40) { return 0.87; }
    else if(pt>   40 && pt <  50) { return 0.93; }
    else if(pt >  50 && pt <  70) { return 0.95; }
    else if(pt >  70 && pt <  100) { return 0.95; }
    else if(pt > 100)             { return 0.96; }
  } else if(fabs(eta)>1.2 && fabs(eta)<2.5) {
    if     (pt >  10 && pt <  20) { return 0.67; }
    else if(pt >  20 && pt <  30) { return 0.79; }
    else if(pt >  30 && pt <  40) { return 0.89; }
    else if(pt>   40 && pt <  50) { return 0.95; }
    else if(pt >  50 && pt <  70) { return 0.96; }
    else if(pt >  70 && pt <  100) { return 0.95; }
    else if(pt > 100)             { return 0.96; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getVetoEleSF(const double pt, const double eta)
{
  if(fabs(eta)<1.4) {
    if     (pt > 10 && pt < 20) { return 1.001; }
    else if(pt > 20 && pt < 30) { return 0.985; }
    else if(pt > 30 && pt < 40) { return 0.993; }
    else if(pt > 40 && pt < 50) { return 0.993; }
    else if(pt > 50 && pt < 60) { return 0.996; }
    else if(pt > 60 && pt < 80) { return 0.993; }
    else if(pt > 80)            { return 0.991; }

  } else if(fabs(eta)<1.6) {
    if     (pt > 10 && pt < 20) { return 1.164; }
    else if(pt > 20 && pt < 30) { return 0.989; }
    else if(pt > 30 && pt < 40) { return 0.993; }
    else if(pt > 40 && pt < 50) { return 0.998; }
    else if(pt > 50 && pt < 60) { return 1.001; }
    else if(pt > 60 && pt < 80) { return 1.044; }
    else if(pt > 80)            { return 1.001; }

  } else if(fabs(eta)<2.1) {
    if     (pt > 10 && pt < 20) { return 1.007; }
    else if(pt > 20 && pt < 30) { return 0.987; }
    else if(pt > 30 && pt < 40) { return 1.010; }
    else if(pt > 40 && pt < 50) { return 1.018; }
    else if(pt > 50 && pt < 60) { return 1.014; }
    else if(pt > 60 && pt < 80) { return 1.004; }
    else if(pt > 80)            { return 1.005; }

  } else if(fabs(eta)<2.4) {
    if     (pt > 10 && pt < 20) { return 1.030; }
    else if(pt > 20 && pt < 30) { return 0.995; }
    else if(pt > 30 && pt < 40) { return 1.006; }
    else if(pt > 40 && pt < 50) { return 1.006; }
    else if(pt > 50 && pt < 60) { return 1.027; }
    else if(pt > 60 && pt < 80) { return 1.004; }
    else if(pt > 80)            { return 0.982; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getTightEleSF(const double pt, const double eta)
{
  if(fabs(eta)<1.4) {
    if     (pt > 10 && pt < 20) { return 1.039; }
    else if(pt > 20 && pt < 30) { return 0.965; }
    else if(pt > 30 && pt < 40) { return 0.976; }
    else if(pt > 40 && pt < 50) { return 0.968; }
    else if(pt > 50 && pt < 60) { return 0.979; }
    else if(pt > 60 && pt < 80) { return 0.968; }
    else if(pt > 80)            { return 0.975; }
    
  } else if(fabs(eta)<1.6) {
    if     (pt > 10 && pt < 20) { return 1.143; }
    else if(pt > 20 && pt < 30) { return 0.971; }
    else if(pt > 30 && pt < 40) { return 0.955; }
    else if(pt > 40 && pt < 50) { return 0.976; }
    else if(pt > 50 && pt < 60) { return 0.994; }
    else if(pt > 60 && pt < 80) { return 1.045; }
    else if(pt > 80)            { return 0.987; }

  } else if(fabs(eta)<2.1) {
    if     (pt > 10 && pt < 20) { return 0.899; }
    else if(pt > 20 && pt < 30) { return 0.936; }
    else if(pt > 30 && pt < 40) { return 0.965; }
    else if(pt > 40 && pt < 50) { return 0.984; }
    else if(pt > 50 && pt < 60) { return 0.988; }
    else if(pt > 60 && pt < 80) { return 0.982; }
    else if(pt > 80)            { return 0.982; }

  } else if(fabs(eta)<2.4) {
    if     (pt > 10 && pt < 20) { return 1.010; }
    else if(pt > 20 && pt < 30) { return 0.976; }
    else if(pt > 30 && pt < 40) { return 0.985; }
    else if(pt > 40 && pt < 50) { return 0.990; }
    else if(pt > 50 && pt < 60) { return 0.999; }
    else if(pt > 60 && pt < 80) { return 0.998; }
    else if(pt > 80)            { return 1.059; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getVetoEleDataIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.4442) {
    if     (pt >  10 && pt <  20) { return 0.76; }
    else if(pt >  20 && pt <  30) { return 0.83; }
    else if(pt >  30 && pt <  40) { return 0.91; }
    else if(pt>   40 && pt <  50) { return 0.95; }
    else if(pt >  50 && pt <  70) { return 0.96; }
    else if(pt >  70 && pt <  100) { return 0.95; }
    else if(pt > 100)             { return 0.95; }
  } else if(fabs(eta)>1.566 && fabs(eta)<2.5) {
    if     (pt >  10 && pt <  20) { return 0.73; }
    else if(pt >  20 && pt <  30) { return 0.80; }
    else if(pt >  30 && pt <  40) { return 0.87; }
    else if(pt>   40 && pt <  50) { return 0.91; }
    else if(pt >  50 && pt <  70) { return 0.92; }
    else if(pt >  70 && pt <  100) { return 0.92; }
    else if(pt > 100)             { return 0.93; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getVetoEleMCIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.4442) {
    if     (pt >  10 && pt <  20) { return 0.75; }
    else if(pt >  20 && pt <  30) { return 0.84; }
    else if(pt >  30 && pt <  40) { return 0.92; }
    else if(pt>   40 && pt <  50) { return 0.96; }
    else if(pt >  50 && pt <  70) { return 0.96; }
    else if(pt >  70 && pt <  100) { return 0.96; }
    else if(pt > 100)             { return 0.97; }
  } else if(fabs(eta)>1.566 && fabs(eta)<2.5) {
    if     (pt >  10 && pt <  20) { return 0.72; }
    else if(pt >  20 && pt <  30) { return 0.81; }
    else if(pt >  30 && pt <  40) { return 0.87; }
    else if(pt>   40 && pt <  50) { return 0.91; }
    else if(pt >  50 && pt <  70) { return 0.91; }
    else if(pt >  70 && pt <  100) { return 0.92; }
    else if(pt > 100)             { return 0.92; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getTightEleDataIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.4442) {
    if     (pt >  10 && pt <  20) { return 0.39; }
    else if(pt >  20 && pt <  30) { return 0.52; }
    else if(pt >  30 && pt <  40) { return 0.62; }
    else if(pt>   40 && pt <  50) { return 0.72; }
    else if(pt >  50 && pt <  70) { return 0.76; }
    else if(pt >  70 && pt <  100) { return 0.79; }
    else if(pt > 100)             { return 0.84; }
  } else if(fabs(eta)>1.566 && fabs(eta)<2.5) {
    if     (pt >  10 && pt <  20) { return 0.32; }
    else if(pt >  20 && pt <  30) { return 0.48; }
    else if(pt >  30 && pt <  40) { return 0.60; }
    else if(pt>   40 && pt <  50) { return 0.69; }
    else if(pt >  50 && pt <  70) { return 0.74; }
    else if(pt >  70 && pt <  100) { return 0.79; }
    else if(pt > 100)             { return 0.82; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getTightEleMCIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.4442) {
    if     (pt >  10 && pt <  20) { return 0.38; }
    else if(pt >  20 && pt <  30) { return 0.53; }
    else if(pt >  30 && pt <  40) { return 0.65; }
    else if(pt>   40 && pt <  50) { return 0.74; }
    else if(pt >  50 && pt <  70) { return 0.78; }
    else if(pt >  70 && pt <  100) { return 0.82; }
    else if(pt > 100)             { return 0.87; }
  } else if(fabs(eta)>1.566 && fabs(eta)<2.5) {
    if     (pt >  10 && pt <  20) { return 0.33; }
    else if(pt >  20 && pt <  30) { return 0.49; }
    else if(pt >  30 && pt <  40) { return 0.62; }
    else if(pt>   40 && pt <  50) { return 0.71; }
    else if(pt >  50 && pt <  70) { return 0.75; }
    else if(pt >  70 && pt <  100) { return 0.78; }
    else if(pt > 100)             { return 0.83; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getMediumPhotonDataIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.4442) {
    if     (pt >  10 && pt <  30) { return 0.6741; }
    else if(pt >  30 && pt <  60) { return 0.8001; }
    else if(pt >  60 && pt <  90) { return 0.7723; }
    else if(pt >  90)             { return 0.796; }
  }

  return 1;
}
//--------------------------------------------------------------------------------------
double getMediumPhotonMCIDEff(const double pt, const double eta)
{
  if(fabs(eta)<1.4442) {
    if     (pt >  10 && pt <  30) { return 0.6992; }
    else if(pt >  30 && pt <  60) { return 0.81034; }
    else if(pt >  60 && pt <  90) { return 0.8045; }
    else if(pt >  90)             { return 0.8134; }
  }

  return 1;
}

//--------------------------------------------------------------------------------------
double getMediumPhotonSF(const double pt, const double eta)
{
  return getMediumPhotonDataIDEff(pt,eta) / getMediumPhotonMCIDEff(pt,eta);
}
