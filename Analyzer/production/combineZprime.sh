#!/bin/bash

sample=$1
dir=$2

declare -a vmass=("50" "75" "100" "125" "150" "200" "250" "300" "400" "500")

if [[ ($sample = "All") || ($sample != "All" && $sample = "VectorDiJet1Jet") ]]
 then 
   for k in "${vmass[@]}"
    do 
      rm "$dir"/VectorDiJet1Jet"$k".root; hadd "$dir"/VectorDiJet1Jet_M"$k".root ../zprimebits/VectorDiJet1Jet_M"$k"_*/*.root     
    done
fi

if [[ ($sample = "All") || ($sample != "All" && $sample = "DMSpin0") ]]
 then
   for k in "${vmass[@]}"
    do
      rm "$dir"/DMSpin0_ggPhibb1j_"$k".root; hadd "$dir"/DMSpin0_ggPhibb1j_"$k".root ../zprimebits/DMSpin0_ggPhibb1j_"$k"_*/*.root
    done
fi

declare -a gmass=("50" "75" "100" "125" "150" "200" "300" "400")

if [[ ($sample = "All") || ($sample != "All" && $sample = "VectorDiJet1Gamma") ]]
 then
   for k in "${gmass[@]}"
    do
      rm "$dir"/VectorDiJet1Gamma_M"$k".root; hadd "$dir"/VectorDiJet1Gamma_M"$k".root ../zprimebits/VectorDiJet1Gamma_"$k"*/*.root
    done
fi

if [[ ($sample = "All") || ($sample != "All" && $sample = "QCD") ]]; then 
    #rm "$dir"/QCD.root; 
    #hadd  "$dir"/QCD.root ../zprimebits/QCD_HT*mc/*.root;
    hadd  "$dir"/QCD_HT100to200.root ../zprimebits/QCD_HT100to200_13TeV*/*.root;
    hadd  "$dir"/QCD_HT200to300.root ../zprimebits/QCD_HT200to300_13TeV*/*.root;
    hadd  "$dir"/QCD_HT300to500.root ../zprimebits/QCD_HT300to500_13TeV*/*.root;
    hadd  "$dir"/QCD_HT500to700.root ../zprimebits/QCD_HT500to700_13TeV*/*.root;
    hadd  "$dir"/QCD_HT700to1000.root ../zprimebits/QCD_HT700to1000_13TeV*/*.root;
    hadd  "$dir"/QCD_HT1000to1500.root ../zprimebits/QCD_HT1000to1500_13TeV*/*.root;
    hadd  "$dir"/QCD_HT1500to2000.root ../zprimebits/QCD_HT1500to2000_13TeV*/*.root;
    hadd  "$dir"/QCD_HT2000toInf.root ../zprimebits/QCD_HT2000toInf_13TeV*/*.root;
fi
if [[ ($sample = "All") || ($sample != "All" && $sample = "DY") ]]; then rm "$dir"/DY.root; hadd  "$dir"/DY.root ../zprimebits/*DYJets*/*.root; fi
if [[ ($sample = "All") || ($sample != "All" && $sample = "W") ]]; then rm "$dir"/W.root; hadd  "$dir"/W.root ../zprimebits/*WJetsToQQ*/*.root; fi
if [[ ($sample = "All") || ($sample != "All" && $sample = "TT") ]]; then rm "$dir"/TTbar_madgraphMLM.root; hadd "$dir"/TTbar_madgraphMLM.root ../zprimebits/*TTJets*/*.root; fi
if [[ ($sample = "All") || ($sample != "All" && $sample = "TT_powheg") ]]; then rm  "$dir"/TTbar_powheg.root; hadd "$dir"/TTbar_powheg.root ../zprimebits/TT_13TeV_powheg_pythia8_ext_mc/*.root; fi
if [[ ($sample = "All") || ($sample != "All" && $sample = "SingleTop") ]]; then rm "$dir"/SingleTop.root; hadd "$dir"/SingleTop.root ../zprimebits/ST*/*.root; fi
if [[ ($sample = "All") || ($sample != "All" && $sample = "JetHT") ]];    then 
    rm "$dir"/JetHTRun2016B.root; hadd "$dir"/JetHTRun2016B.root ../zprimebits/JetHTRun2016B*_data/*.root;
    rm "$dir"/JetHTRun2016C.root; hadd "$dir"/JetHTRun2016C.root ../zprimebits/JetHTRun2016C*_data/*.root; 
    rm "$dir"/JetHTRun2016D.root; hadd "$dir"/JetHTRun2016D.root ../zprimebits/JetHTRun2016D*_data/*.root; 
    rm "$dir"/JetHTRun2016E.root; hadd "$dir"/JetHTRun2016E.root ../zprimebits/JetHTRun2016E*_data/*.root; 
    rm "$dir"/JetHTRun2016F.root; hadd "$dir"/JetHTRun2016F.root ../zprimebits/JetHTRun2016F*_data/*.root; 
fi

if [[ ($sample = "All") || ($sample != "All" && $sample = "Diboson") ]]; then
  hadd "$dir"/WWTo4Q_13TeV_amcatnlo.root               ../zprimebits/WWTo4Q_13TeV_amcatnlo_mc/*.root;
  #hadd "$dir"/WWTo4Q_13TeV_powheg                /*.root;
  hadd "$dir"/ZZTo4Q_13TeV_amcatnlo.root               ../zprimebits/ZZTo4Q_13TeV_amcatnlo_mc/*.root;
  hadd "$dir"/WZ_13TeV.root                           ../zprimebits/WZ_13TeV_mc/*.root;
fi

if [[ ($sample = "All") || ($sample != "All" && $sample = "Hbb") ]]; then
  hadd "$dir"/GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8.root      ../zprimebits//GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia_mc/*.root ;
  hadd "$dir"/GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8_ext.root  ../zprimebits//GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8_ext_mc/*.root;
  hadd "$dir"/GluGluHToBB_M125_13TeV_powheg_herwigpp.root           ../zprimebits//GluGluHToBB_M125_13TeV_powheg_herwigpp_mc/*.root;
  hadd "$dir"/GluGluHToBB_M125_13TeV_powheg_pythia8.root            ../zprimebits/GluGluHToBB_M125_13TeV_powheg_pythia8_mc/*.root;
  hadd "$dir"/VBFHToBB_M125_13TeV_amcatnlo_pythia8.root             ../zprimebits/VBFHToBB_M125_13TeV_amcatnlo_pythia8_mc/*.root;
  hadd "$dir"/VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix.root    ../zprimebits/VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_mc/*.root;
  hadd "$dir"/VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_ext.root ../zprimebits/VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_ext_mc/*.root;
  hadd "$dir"/ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8.root    ../zprimebits/ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_mc/*.root;
fi


