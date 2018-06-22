#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
from optparse import OptionParser
from submitZprime import samplesDict, exec_me

EOS = ''#eos root://cmseos.fnal.gov'
cmssw = os.getenv('CMSSW_VERSION', 'CMSSW_9_4_7')
cmssw_base = os.getenv('CMSSW_BASE', 'CMSSW_9_4_7')

import ROOT
normDict = {'DYJetsToQQ_HT180_13TeV': 'DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8',
            'WJetsToQQ_HT180_13TeV': 'WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8',
            'WJetsToQQ_HT400to600_TuneCP5_13TeV_noPF': 'WJetsToQQ_HT400to600_TuneCP5_13TeV',
            'WJetsToQQ_HT600to800_TuneCP5_13TeV_noPF': 'WJetsToQQ_HT600to800_TuneCP5_13TeV',
            'WJetsToQQ_HT_800toInf_TuneCP5_13TeV_noPF': 'WJetsToQQ_HT_800toInf_TuneCP5_13TeV',
            'QCD_HT1000to1500_TuneCP5_13TeV_madgraph_pythia8': 'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT100to200_TuneCP5_13TeV_madgraph_pythia8': 'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1500to2000_TuneCP5_13TeV_madgraph_pythia8': 'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT2000toInf_TuneCP5_13TeV_madgraph_pythia8': 'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT200to300_TuneCP5_13TeV_madgraph_pythia8': 'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT300to500_TuneCP5_13TeV_madgraph_pythia8': 'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT500to700_TuneCP5_13TeV_madgraph_pythia8': 'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT700to1000_TuneCP5_13TeV_madgraph_pythia8': 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', 
            'QCD_HT300to500_TuneCP5_13TeV_madgraph_pythia8_noPF_byLumi': 'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT500to700_TuneCP5_13TeV_madgraph_pythia8_noPF_byLumi': 'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT700to1000_TuneCP5_13TeV_madgraph_pythia8_noPF_byLumi': 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT50to100_13TeV': 'QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT100to200_13TeV': 'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT100to200_13TeV_ext': 'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT100to200_13TeV_all': 'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT200to300_13TeV': 'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT200to300_13TeV_ext': 'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT200to300_13TeV_all': 'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT300to500_13TeV': 'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT300to500_13TeV_ext': 'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT300to500_13TeV_all': 'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT500to700_13TeV': 'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT500to700_13TeV_ext': 'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT500to700_13TeV_all': 'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT700to1000_13TeV': 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT700to1000_13TeV_ext': 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT700to1000_13TeV_all': 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1000to1500_13TeV': 'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1000to1500_13TeV_ext': 'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1000to1500_13TeV_all': 'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1500to2000_13TeV': 'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1500to2000_13TeV_ext': 'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1500to2000_13TeV_all': 'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT2000toInf_13TeV': 'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT2000toInf_13TeV_ext': 'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT2000toInf_13TeV_all': 'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'ST_s_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1': 'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
            'ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin': 'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
            'ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin': 'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
            'ST_tW_antitop_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4': 'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
            'ST_tW_top_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4': 'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
            'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8': 'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
            'ST_s_channel_4f_leptonDecays_TuneCP5_13TeV_amcatnlo_pythia8_noPF': 'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
            'ST_t_channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8': 'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
            'ST_t_channel_top_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8': 'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
            'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8': 'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
            'TTJets_13TeV': 'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'TT_powheg': 'TT_TuneCUETP8M1_13TeV-powheg-pythia8',
            'TT_13TeV_powheg_pythia8_ext': 'TT_TuneCUETP8M1_13TeV-powheg-pythia8',
            'TTToHadronic_TuneCP5_13TeV_powheg_pythia8_byLumi': 'TTToHadronic_TuneCP5_13TeV_powheg_pythia8',
            'TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8_byLumi': 'TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8',
            'WJetsToQQ_HT_600ToInf_13TeV': 'WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',            
            'WJetsToLNu_HT_70To100_13TeV': 'WJetsToLNu_HT_70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_70To100_13TeV_ext': 'WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_100To200_13TeV': 'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_100To200_13TeV_ext': 'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_100To200_13TeV_ext1': 'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_100To200_13TeV_ext2': 'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_200To400_13TeV': 'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_200To400_13TeV_ext': 'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_200To400_13TeV_ext1': 'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_200To400_13TeV_ext2': 'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_400To600_13TeV': 'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_400To600_13TeV_ext': 'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_600To800_13TeV': 'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_600To800_13TeV_ext': 'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_600To800_13TeV_ext1': 'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_800To1200_13TeV': 'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_800To1200_13TeV_ext': 'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_800To1200_13TeV_ext1': 'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_1200To2500_13TeV': 'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_1200To2500_13TeV_ext': 'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_1200To2500_13TeV_ext1': 'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_2500ToInf_13TeV': 'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT_2500ToInf_13TeV_ext': 'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',       
            'WJetsToLNu_HT_2500ToInf_13TeV_ext1': 'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',           
            'ZJetsToQQ_HT600toInf_13TeV_madgraph': 'ZJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'ZJetsToQQ_HT400to600_TuneCP5_13TeV_noPF': 'ZJetsToQQ_HT400to600_TuneCP5_13TeV',
            'ZJetsToQQ_HT600to800_3j_TuneCP5_13TeV_noPF': 'ZJetsToQQ_HT600to800_3j_TuneCP5_13TeV',
            'ZJetsToQQ_HT_800toInf_TuneCP5_13TeV_noPF': 'ZJetsToQQ_HT_800toInf_TuneCP5_13TeV',
            'WWTo4Q_13TeV_amcatnlo': 'WWTo4Q_4f_13TeV_amcatnloFXFX_madspin_pythia8',
            'WWTo4Q_13TeV_powheg': 'WWTo4Q_13TeV_powheg',
            'WW_TuneCP5_13TeV_pythia8': 'WW_TuneCUETP8M1_13TeV-pythia8',
            'WZ_TuneCP5_13TeV_pythia8': 'WZ_TuneCUETP8M1_13TeV-pythia8',
            'ZZ_TuneCP5_13TeV_pythia8': 'ZZ_TuneCUETP8M1_13TeV-pythia8',
            'WZ_13TeV_pythia8': 'WZ_TuneCUETP8M1_13TeV-pythia8',
            'ZZTo4Q_13TeV_amcatnlo': 'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8',
            'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8': 'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8',
            'GluGluHToBB_M125_13TeV_powheg_herwigpp': 'GluGluHToBB_M125_13TeV_powheg_herwigpp',
            'GluGluHToBB_M125_13TeV_powheg_pythia8': 'GluGluHToBB_M125_13TeV_powheg_pythia8',
            'GluGluHToBB_M125_13TeV_powheg_pythia8_ext': 'GluGluHToBB_M125_13TeV_powheg_pythia8',
            'VBFHToBB_M125_13TeV_amcatnlo_pythia8': 'VBFHToBB_M125_13TeV_amcatnlo_pythia8',
            'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix': 'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix',
            'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8',
            'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8': 'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8',
            'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
            'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_ext': 'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
            'Spin0_ggPhi12j_g1_100_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_100',
            'Spin0_ggPhi12j_g1_100_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_100_Pseudo',
            'Spin0_ggPhi12j_g1_125_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_125_Pseudo',
            'Spin0_ggPhi12j_g1_125_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_125',
            'Spin0_ggPhi12j_g1_150_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_150_Pseudo',
            'Spin0_ggPhi12j_g1_150_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_150',
            'Spin0_ggPhi12j_g1_200_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_200_Pseudo',
            'Spin0_ggPhi12j_g1_200_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_200',
            'Spin0_ggPhi12j_g1_250_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_250_Pseudo',
            'Spin0_ggPhi12j_g1_250_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_250',
            'Spin0_ggPhi12j_g1_25_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_25_Pseudo',
            'Spin0_ggPhi12j_g1_25_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_25',
            'Spin0_ggPhi12j_g1_300_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_300_Pseudo',
            'Spin0_ggPhi12j_g1_300_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_300',
            'Spin0_ggPhi12j_g1_350_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_350_Pseudo',
            'Spin0_ggPhi12j_g1_350_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_350',
            'Spin0_ggPhi12j_g1_400_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_400',
            'Spin0_ggPhi12j_g1_400_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_400_Pseudo',
            'Spin0_ggPhi12j_g1_500_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_500_Pseudo',
            'Spin0_ggPhi12j_g1_500_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_500',
            'Spin0_ggPhi12j_g1_600_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_600_Pseudo',
            'Spin0_ggPhi12j_g1_600_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_600',
            'Spin0_ggPhi12j_g1_800_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_800_Pseudo',
            'Spin0_ggPhi12j_g1_800_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_800',
            'Spin0_ggPhi12j_g1_1000_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_1000_Pseudo',
            'Spin0_ggPhi12j_g1_1000_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_1000',
            'Spin0_ggPhi12j_g1_5_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_5',
            'Spin0_ggPhi12j_g1_5_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_5_Pseudo',
            'Spin0_ggPhi12j_g1_50_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_50',
            'Spin0_ggPhi12j_g1_50_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_50_Pseudo',
            'Spin0_ggPhi12j_g1_75_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_75_Pseudo',
            'Spin0_ggPhi12j_g1_75_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_75',
            'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8',
            'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8',
            'WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8': 'WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8',
            'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8': 'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8',
            'ttHTobb_M125_13TeV_powheg_pythia8': 'ttHTobb_M125_13TeV_powheg_pythia8',
            'ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV_powheg_pythia8': 'ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV_powheg_pythia8',
            'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
            'ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8',
            'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8': 'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8',
            'ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8',
            'bbHToBB_M_125_4FS_yb2_13TeV_amcatnlo': 'bbHToBB_M_125_4FS_yb2_13TeV_amcatnlo',
            'bbHToBB_M_125_4FS_ybyt_13TeV_amcatnlo': 'bbHToBB_M_125_4FS_ybyt_13TeV_amcatnlo',
            'VectorDiJet1Jet_100_madgraph_2016':'VectorDiJet1Jet_M100',
            'VectorDiJet1Jet_125_madgraph_2017_noPF':'VectorDiJet1Jet_madgraph_Mphi125Mchi1_13TeV',
            'VectorDiJet1Jet_75_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi75Mchi1_13TeV',
            'VectorDiJet1Jet_50_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi50Mchi1_13TeV',
            'VectorDiJet1Jet_100_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi100Mchi1_13TeV',
            'VectorDiJet1Jet_115_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi115Mchi1_13TeV',
            'VectorDiJet1Jet_150_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi150Mchi1_13TeV',
            'VectorDiJet1Jet_175_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi175Mchi1_13TeV',
            'VectorDiJet1Jet_200_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi200Mchi1_13TeV',
            'VectorDiJet1Jet_250_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi250Mchi1_13TeV',
            'VectorDiJet1Jet_300_madgraph_2017_noPF' : 'VectorDiJet1Jet_madgraph_Mphi300Mchi1_13TeV',
            'VectorDiJet1Jet_100_13TeV_madgraph':'VectorDiJet1Jet_M100',
            'VectorDiJet1Jet_125_13TeV_madgraph':'VectorDiJet1Jet_M125',
            'VectorDiJet1Jet_150_13TeV_madgraph':'VectorDiJet1Jet_M150',
            'VectorDiJet1Jet_200_13TeV_madgraph':'VectorDiJet1Jet_M200',
            'VectorDiJet1Jet_25_13TeV_madgraph':'VectorDiJet1Jet_M25',
            'VectorDiJet1Jet_300_13TeV_madgraph':'VectorDiJet1Jet_M300',
            'VectorDiJet1Jet_400_13TeV_madgraph':'VectorDiJet1Jet_M400',
            'VectorDiJet1Jet_500_13TeV_madgraph':'VectorDiJet1Jet_M500',
            'VectorDiJet1Jet_50_13TeV_madgraph':'VectorDiJet1Jet_M50',
            'VectorDiJet1Jet_600_13TeV_madgraph':'VectorDiJet1Jet_M600',
            'VectorDiJet1Jet_800_13TeV_madgraph':'VectorDiJet1Jet_M800',
            'VectorDiJet1Jet_1000_13TeV_madgraph':'VectorDiJet1Jet_M1000',
            'VectorDiJet1Jet_125_13TeV_madgraph':'VectorDiJet1Jet_M125',
            'VectorDiJet1Jet_75_13TeV_madgraph':'VectorDiJet1Jet_M75',

            'VectorDiJet1Jet_madgraph_Mphi50Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi50Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi75Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi75Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi100Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi100Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi115Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi115Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi125Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi125Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi150Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi150Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi175Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi175Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi200Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi200Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi250Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi250Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi300Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi300Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi350Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi350Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi400Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi400Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi450Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi450Mchi1_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi500Mchi1_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi500Mchi1_13TeV',
            
            'VectorDiJet1Jet_madgraph_Mphi50Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi50Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi75Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi75Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi100Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi100Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi115Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi115Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi125Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi125Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi150Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi150Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi175Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi175Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi200Mchi3000_13TeV_noPF': 'VectorDiJet1Jet_madgraph_Mphi200Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi225Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi225Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi250Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi250Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi300Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi300Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi350Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi350Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi400Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi400Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi450Mchi3000_13TeV_noPF_8': 'VectorDiJet1Jet_madgraph_Mphi450Mchi3000_13TeV',
            'VectorDiJet1Jet_madgraph_Mphi500Mchi3000_13TeV_noPF': 'VectorDiJet1Jet_madgraph_Mphi500Mchi3000_13TeV',
            
            'DYJetsToLL_M_50_13TeV_ext': 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8': 'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8',
            'WW_13TeV_pythia8': 'WW_TuneCUETP8M1_13TeV-pythia8',
            'ZZ_13TeV_pythia8': 'ZZ_TuneCUETP8M1_13TeV-pythia8',

            'BulkGravTohhTohVVhbb_narrow_M_1000_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_1000_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_1200_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_1200_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_1400_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_1400_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_1600_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_1600_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_1800_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_1800_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_2000_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_2000_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_2500_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_2500_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_3000_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_3000_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_3500_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_3500_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_4000_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_4000_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_4500_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_4500_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_600_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_600_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_650_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_650_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_700_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_700_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_800_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_800_13TeV_madgraph',
            'BulkGravTohhTohVVhbb_narrow_M_900_13TeV_madgraph_8X': 'BulkGravTohhTohVVhbb_narrow_M_900_13TeV_madgraph',
            
            }

filesToTransfer = "%s.tgz, %s/bin/slc6_amd64_gcc630/NormalizeNtuple, %s/src/BaconAnalyzer/Analyzer/data.tgz, %s/src/BaconAnalyzer/Analyzer/production/skimmer.py, %s/src/BaconAnalyzer/Analyzer/production/skimmerDDT.py, %s/src/BaconAnalyzer/Analyzer/production/skimmerN2.py, %s/src/BaconAnalyzer/Analyzer/production/skimmerWtag.py, %s/src/BaconAnalyzer/Analyzer/production/submitZprime.py"%(cmssw_base,cmssw_base,cmssw_base,cmssw_base,cmssw_base,cmssw_base,cmssw_base,cmssw_base)

def justHadd(options,args):    
    DataDir = options.idir
    OutDir = options.idir
    
    samples = samplesDict[options.sample]
    
    postfix = ''
    exec_me('mkdir -p $PWD/hadd_jobs/',options.dryRun)
    exec_me('mkdir -p /%s/hadd'%(OutDir),options.dryRun) 
    for label, isMc in samples.iteritems():
        basename = label + '.root'
        
        filesToConvert = []
        badFiles = []
        filesToConvert, badFiles = getFilesRecursively(DataDir,label+'/',None,None)
        print "files To Convert = ",filesToConvert
        print "bad files = ", badFiles
        cwd = os.getcwd()

        haddAll = False
        haddOutExistsList = []
        haddOutList = []
        for i in range(0,len(filesToConvert)/50+1):
            if not os.path.isfile(OutDir+'/hadd/'+basename.replace('.root','_%i.root'%i)):
                haddOutExistsList.append(False)
                haddOutList.append(OutDir+'/hadd/'+basename.replace('.root','_%i.root'%i))
                haddCommand = '#!/bin/bash\n'
                haddCommand += 'source /cvmfs/cms.cern.ch/cmsset_default.sh\n'
                haddCommand += 'pwd\n'
                haddCommand += 'tar -xf %s.tgz\n'% (cmssw)
                haddCommand += 'rm %s.tgz\n'% (cmssw)
                haddCommand += 'export SCRAM_ARCH=slc6_amd64_gcc630\n'
                haddCommand += 'mkdir -p %s/src\n'% (cmssw)
                haddCommand += 'scramv1 project CMSSW CMSSW_9_2_12\n'
                haddCommand += 'cd CMSSW_9_2_12/src\n'
                haddCommand += 'scram b ProjectRename\n'
                haddCommand += 'eval `scramv1 runtime -sh`\n'
                haddCommand += 'pwd\n'
                haddCommand += 'mkdir -p $PWD/hadd\n'       
                haddCommand += 'hadd -f hadd/%s %s\n'%(basename.replace('.root','_%i.root'%i),(' '.join(filesToConvert[i*50:(i+1)*50])))
                haddCommand += 'xrdcp -s $PWD/hadd/%s root://cmseos.fnal.gov//%s/hadd/%s\n'%(basename.replace('.root','_%i.root'%i),OutDir,basename.replace('.root','_%i.root'%i))
                haddCommand += 'rm -r $PWD/hadd\n'
                with open('hadd_jobs/hadd_command_%s.sh'%(basename.replace('.root','_%i.root'%i)),'w') as f:
                    f.write(haddCommand)
                os.system('rm -f hadd_jobs/hadd_command_%s.stdout' % basename.replace('.root','_%i.root'%i))
                os.system('rm -f hadd_jobs/hadd_command_%s.stderr' % basename.replace('.root','_%i.root'%i))
                os.system('rm -f hadd_jobs/hadd_command_%s.log' % basename.replace('.root','_%i.root'%i))
                os.system('rm -f hadd_jobs/hadd_command_%s.jdl'% basename.replace('.root','_%i.root'%i))
                condor_file = open('hadd_jobs/hadd_command_%s.jdl' % basename.replace('.root','_%i.root'%i), 'w')
                condor_file.write('universe = vanilla\n')
                condor_file.write('Executable = hadd_jobs/hadd_command_%s.sh\n'% basename.replace('.root','_%i.root'%i))
                condor_file.write('Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )\n')
                #condor_file.write('request_disk = 3000000\n')
                condor_file.write('request_disk = 5000000\n')
                #condor_file.write('request_memory = 5000\n')
                condor_file.write('request_memory = 6000\n')
                condor_file.write('Should_Transfer_Files = YES\n')
                condor_file.write('WhenToTransferOutput = ON_EXIT\n')
                condor_file.write('Transfer_Input_Files = %s\n'%filesToTransfer)
                condor_file.write('use_x509userproxy = true\n')
                condor_file.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
                condor_file.write('Output = %s.stdout\n' % os.path.abspath(condor_file.name))
                condor_file.write('Error = %s.stderr\n' % os.path.abspath(condor_file.name))
                condor_file.write('Log = %s.log\n' % os.path.abspath(condor_file.name))
                condor_file.write('Queue 1\n')
                condor_file.close()
                os.system('chmod +x %s'% os.path.abspath(condor_file.name))
                exec_me('condor_submit %s'%(os.path.abspath(condor_file.name)),options.dryRun)
            else:
                haddOutExistsList.append(True)
                haddOutList.append(OutDir+'/hadd/'+basename.replace('.root','_%i.root'%i))
                
        haddAll = all(haddOutExistsList)
        if haddAll:
            haddCommand = '#!/bin/bash\n'
            haddCommand += 'source /cvmfs/cms.cern.ch/cmsset_default.sh\n'
            haddCommand += 'pwd\n'
            haddCommand += 'tar -xf %s.tgz\n'% (cmssw)
            haddCommand += 'rm %s.tgz\n'% (cmssw)
            haddCommand += 'export SCRAM_ARCH=slc6_amd64_gcc630\n'
            haddCommand += 'mkdir -p %s/src\n'% (cmssw)
            haddCommand += 'scramv1 project CMSSW %s\n'% (cmssw)
            haddCommand += 'cd %s/src\n'% (cmssw)
            haddCommand += 'scram b ProjectRename\n'
            haddCommand += 'eval `scramv1 runtime -sh`\n'
            haddCommand += 'pwd\n'
            haddCommand += 'mkdir -p $PWD/hadd\n'        
            haddCommand += 'hadd -f hadd/%s %s\n'%(basename,(' '.join(haddOutList)))
            haddCommand += 'xrdcp -s $PWD/hadd/%s root://cmseos.fnal.gov/%s/hadd/%s\n'%(basename.replace('.root','_%i.root'%i),OutDir,basename.replace('.root','_%i.root'%i))
            haddCommand += 'rm -r $PWD/hadd\n'
            with open('hadd_jobs/hadd_command_%s.sh'%(basename),'w') as f:
                f.write(haddCommand)
            os.system('rm -f hadd_jobs/hadd_command_%s.stdout' % basename)
            os.system('rm -f hadd_jobs/hadd_command_%s.stderr' % basename)
            os.system('rm -f hadd_jobs/hadd_command_%s.log' % basename)
            os.system('rm -f hadd_jobs/hadd_command_%s.jdl'% basename)
            condor_file = open('hadd_jobs/hadd_command_%s.jdl' % basename, 'w')
            condor_file.write('universe = vanilla\n')
            condor_file.write('Executable = hadd_jobs/hadd_command_%s.sh\n'% basename)
            condor_file.write('Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )\n')
            #condor_file.write('request_disk = 3000000\n')
            condor_file.write('request_disk = 5000000\n')
            #condor_file.write('request_memory = 5000\n')
            condor_file.write('request_memory = 8000\n')
            condor_file.write('Should_Transfer_Files = YES\n')
            condor_file.write('WhenToTransferOutput = ON_EXIT\n')
            condor_file.write('Transfer_Input_Files = %s\n'%filesToTransfer)
            condor_file.write('use_x509userproxy = true\n')
            condor_file.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
            condor_file.write('Output = %s.stdout\n' % os.path.abspath(condor_file.name))
            condor_file.write('Error = %s.stderr\n' % os.path.abspath(condor_file.name))
            condor_file.write('Log = %s.log\n' % os.path.abspath(condor_file.name))
            condor_file.write('Queue 1\n')
            condor_file.close()
            os.system('chmod +x %s'% os.path.abspath(condor_file.name))
            exec_me('condor_submit %s'%(os.path.abspath(condor_file.name)),options.dryRun)

def addSklimCommand(iBasename,isMc,options):
    haddCommand = ''
    basename = iBasename

    DDTOn = False
    N2On = False
    WtagOn = False
    if options.justDDT:
        DDTOn = True
    elif options.justN2:
        N2On = True
    elif options.justWtag:
        WtagOn = True
    else:
        DDTOn = False
        N2On = False
        WtagOn = False
    if DDTOn:
        haddCommand += 'python skimmerDDT.py -i $PWD/hadd/ -o $PWD/sklim/ -s %s\n'%(basename.replace('.root',''))
    elif N2On:
        haddCommand += 'python skimmerN2.py -i $PWD/hadd/ -o $PWD/sklim/ -s %s \n'%(basename.replace('.root',''))
    elif WtagOn:
        if isMc=='mc':
            #haddCommand += 'python skimmerWtag.py -i $PWD/hadd/ -o $PWD/sklim/ -s %s --isMc --isPu --jet %s --ddt %s --iddt %s \n'%(basename.replace('.root',''),options.jet,options.ddt,options.iddt)
            haddCommand += 'python skimmerWtag.py -i $PWD/hadd/ -o $PWD/sklim/ -s %s --isMc --jet %s --ddt %s --iddt %s \n'%(basename.replace('.root',''),options.jet,options.ddt,options.iddt)
        else:
            haddCommand += 'python skimmerWtag.py -i $PWD/hadd/ -o $PWD/sklim/ -s %s --isPu --jet %s --ddt %s --iddt %s \n'%(basename.replace('.root',''),options.jet,options.ddt,options.iddt)
    else:
        haddCommand += 'python skimmer.py -i $PWD/hadd/ -o $PWD/sklim/ -s %s\n'%(basename.replace('.root',''))
    return haddCommand

def main(options,args):
    global DataDir,OutDir;
    DataDir = options.idir
    OutDir = options.idir
    
    samples = samplesDict[options.sample]
    
    loopOverFiles = options.loop
    if options.haddlarge: loopOverFiles = True

    EOS = ''#eos root://cmseos.fnal.gov'
    postfix = ''

    exec_me('mkdir -p $PWD/hadd_jobs/',options.dryRun)
    exec_me('mkdir -p /%s/%s'%(OutDir,haddDir),options.dryRun)
    exec_me('mkdir -p /%s/%s'%(OutDir,sklimDir),options.dryRun)
    exec_me('mkdir -p /%s/%s'%(OutDir,normDir),options.dryRun)
    if options.justDDT:
        exec_me('mkdir -p /%s/sklimDDT'%(OutDir),options.dryRun)
        exec_me('mkdir -p /%s/normDDT'%(OutDir),options.dryRun)
    if options.justN2:
        exec_me('mkdir -p /%s/sklimN2'%(OutDir),options.dryRun)
        exec_me('mkdir -p /%s/normN2'%(OutDir),options.dryRun)
    # if options.justWtag:
    #     exec_me('mkdir -p /%s/sklimWtag%s'%(OutDir,options.jet),options.dryRun)

    for label, isMc in samples.iteritems():
        basename = label + '.root'
        if options.job>-1:
            basename = label + '_%i.root'%options.job

        haddOn = True
        sklimOn = True
        normOn = True
        DDTOn = False
        N2On = False
        WtagOn = False            
        if os.path.isfile(OutDir+'/'+haddDir+'/'+basename):
            haddOn = False
        if os.path.isfile(OutDir+'/'+sklimDir+'/'+basename):
            sklimOn = False
        if os.path.isfile(OutDir+'/'+normDir+'/'+basename.replace('.root','_1000pb_weighted.root')):
            normOn = False

        if options.justNorm:
            haddOn = False
            sklimOn = False
            normOn = True
        elif options.justSklim:
            haddOn = False
            sklimOn = True
            normOn = False
        elif 'runPu' in options.executable:
            haddOn = True
            sklimOn = False
            normOn = False

        if options.justDDT:
            DDTOn = True
        elif options.justN2:
            haddOn = False
            normOn = True
            N2On = True
        elif options.justWtag:
            haddOn = False
            #sklimOn = True
            WtagOn = True
            normOn = False
        else:
            DDTOn = False
            N2On = False
            WtagOn = False

        filesToConvert = []
        badFiles = []
        if haddOn or loopOverFiles:
            filesToConvert, badFiles = getFilesRecursively(DataDir,label+'/',None,None)
            #filesToConvert, badFiles = getFilesRecursively(DataDir+'/'+haddDir+'/',label+'/',None,None)
        print "files To Convert = ",filesToConvert
        print "bad files = ", badFiles
        cwd = os.getcwd()
        haddCommand = '#!/bin/bash\n'
        haddCommand += 'source /cvmfs/cms.cern.ch/cmsset_default.sh\n'
        haddCommand += 'pwd\n'
        haddCommand += 'tar -xf %s.tgz\n'% (cmssw)
        haddCommand += 'rm %s.tgz\n'% (cmssw)
        haddCommand += 'export SCRAM_ARCH=slc6_amd64_gcc630\n'
        haddCommand += 'mkdir -p %s/src\n'% (cmssw)
        haddCommand += 'scramv1 project CMSSW CMSSW_9_2_12\n'
        haddCommand += 'cd CMSSW_9_2_12/src\n'
        haddCommand += 'scram b ProjectRename\n'
        haddCommand += 'eval `scramv1 runtime -sh`\n'
        haddCommand += 'pwd\n'
        haddCommand += 'cp ../../data.tgz .\n'
        haddCommand += 'mkdir -p ${PWD}/BaconAnalyzer/Analyzer/\n'
        haddCommand += 'tar -xvzf data.tgz -C ${PWD}/BaconAnalyzer/Analyzer/\n'
        haddCommand += 'mkdir -p $PWD/sklim\n'
        haddCommand += 'cp ../../submitZprime.py .\n'
        haddCommand += 'cp ../../skimmer.py .\n'
        if sklimOn:
            if DDTOn:
                haddCommand += 'cp ../../skimmerDDT.py .\n'
            if N2On:
                haddCommand += 'cp ../../skimmerN2.py .\n'
            if WtagOn:
                haddCommand += 'cp ../../skimmerWtag.py .\n'
        haddCommand += 'eval `scramv1 runtime -sh`\n'
        haddCommand += 'mkdir -p $PWD/hadd\n'
        haddCommand += 'mkdir -p $PWD/sklim\n'
        if haddOn:
            print "files len = ",len(filesToConvert)/50+1
            for i in range(0,len(filesToConvert)/50+1):
                haddCommand += 'mkdir -p $PWD/hadd\n'
                haddCommand += 'mkdir -p $PWD/sklim\n'
                print 'hadd -f hadd/%s %s\n'%(basename.replace('.root','_%i.root'%i),(' '.join(filesToConvert[i*50:(i+1)*50])))
                haddCommand += 'hadd -f hadd/%s %s\n'%(basename.replace('.root','_%i.root'%i),(' '.join(filesToConvert[i*50:(i+1)*50])))
                if options.haddlarge:
                    haddCommand += 'xrdcp $PWD/hadd/%s root://cmseos.fnal.gov/%s/%s/%s\n'%(basename.replace('.root','_%i.root'%i),OutDir,haddDir,basename.replace('.root','_%i.root'%i))
                    haddCommand += 'rm -r -f hadd/*'
            else:
                haddCommand += 'hadd -f $PWD/hadd/%s $PWD/hadd/%s\n'%(basename,basename.replace('.root','_*.root'))
                haddCommand += 'rm $PWD/hadd/%s\n'%(basename.replace('.root','_*.root'))
                haddCommand += 'xrdcp $PWD/hadd/%s root://cmseos.fnal.gov/%s/%s/%s\n'%(basename,OutDir,haddDir,basename)
        if sklimOn:
            if loopOverFiles:
                if not haddOn:
                    print "files len = ",len(filesToConvert)/50+1
                    for i in range(0,len(filesToConvert)/50+1):
                        haddCommand += 'mkdir -p $PWD/sklim\n'
                        haddCommand += 'mkdir -p $PWD/hadd\n'
                        haddCommand += 'xrdcp root://cmseos.fnal.gov//%s/%s/%s $PWD/hadd/%s\n'%(OutDir,haddDir,basename.replace('.root','_%i.root'%i),basename.replace('.root','_%i.root'%i))
                    haddCommand += addSklimCommand(basename,isMc,options)
                    for i in range(0,len(filesToConvert)/50+1):
                        haddCommand += 'xrdcp -s $PWD/sklim/%s root://cmseos.fnal.gov//%s/%s/%s\n'%(basename.replace('.root','_%i.root'%i),OutDir,sklimDir,basename.replace('.root','_%i.root'%i))
                else:
                    haddCommand += addSklimCommand(basename,isMc,options)
                    for i in range(0,len(filesToConvert)/50+1):
                        haddCommand += 'xrdcp -s $PWD/sklim/%s root://cmseos.fnal.gov//%s/%s/%s\n'%(basename.replace('.root','_%i.root'%i),OutDir,sklimDir,basename.replace('.root','_%i.root'%i))
            else:
                print 'sklimon'
                if not haddOn:
                    if not WtagOn:
                        haddCommand += 'xrdcp root://cmseos.fnal.gov//%s/%s/%s $PWD/hadd/%s\n'%(OutDir,haddDir,basename,basename)
                    else:
                        print  'xrdcp root://cmseos.fnal.gov//%s/normN2/%s $PWD/hadd/%s\n'%(OutDir,basename.replace('.root','_1000pb_weighted.root'),basename)
                        if isMc=='mc':
                            haddCommand += 'xrdcp root://cmseos.fnal.gov//%s/normN2/%s $PWD/hadd/%s\n'%(OutDir,basename.replace('.root','_1000pb_weighted.root'),basename)
                        else:
                            haddCommand += 'xrdcp root://cmseos.fnal.gov//%s/sklimN2/%s $PWD/hadd/%s\n'%(OutDir,basename,basename)
                haddCommand += addSklimCommand(basename,isMc,options)
                haddCommand += 'xrdcp -s $PWD/sklim/%s root://cmseos.fnal.gov//%s/%s/%s\n'%(basename,OutDir,sklimDir,basename)
        if isMc=='mc' and normOn:            
            haddCommand += 'cp ../../NormalizeNtuple .\n'
            #if not sklimOn:
                #haddCommand += 'mkdir -p $PWD/sklim\n'
                #haddCommand += 'xrdcp root://cmseos.fnal.gov//%s/sklim/%s $PWD/sklim/%s\n'%(OutDir,basename,basename)                
            #haddCommand += 'mkdir -p /uscmst1b_scratch/lpc1/3DayLifetime/cmantill/sklim\n'
            #haddCommand += 'xrdcp root://cmseos.fnal.gov//%s/%s/%s /uscmst1b_scratch/lpc1/3DayLifetime/cmantill/sklim/%s\n'%(OutDir,sklimDir,basename,basename)
            #haddCommand += 'echo "%s\t/uscmst1b_scratch/lpc1/3DayLifetime/cmantill/sklim/%s" > normlist.txt\n'%(normDict[basename.replace('.root','')],basename)
            haddCommand += 'echo "%s\troot://cmseos.fnal.gov//%s/%s/%s" > normlist.txt\n'%(normDict[basename.replace('.root','')],OutDir,sklimDir,basename)
            haddCommand += 'NormalizeNtuple normlist.txt\n'
            haddCommand += 'xrdcp root://cmseos.fnal.gov//%s/%s/%s root://cmseos.fnal.gov//%s/%s/%s\n'%(OutDir,sklimDir,basename.replace('.root','_1000pb_weighted.root'),OutDir,normDir,basename.replace('.root','_1000pb_weighted.root'))
            haddCommand += 'rm -f root://cmseos.fnal.gov//%s/%s/%s\n'%(OutDir,sklimDir,basename.replace('.root','_1000pb_weighted.root'))
            #haddCommand += 'xrdcp /uscmst1b_scratch/lpc1/3DayLifetime/cmantill/sklim/%s root://cmseos.fnal.gov//%s/%s/%s\n'%(basename.replace('.root','_1000pb_weighted.root'),OutDir,normDir,basename.replace('.root','_1000pb_weighted.root'))
        haddCommand += 'rm -r $PWD/hadd\n'
        haddCommand += 'rm -r $PWD/sklim\n'
            
            
        with open('hadd_jobs/hadd_command_%s.sh'%(basename),'w') as f:
            f.write(haddCommand)

        if haddOn or sklimOn or (isMc=='mc' and normOn):
            os.system('rm -f hadd_jobs/hadd_command_%s.stdout' % basename)
            os.system('rm -f hadd_jobs/hadd_command_%s.stderr' % basename)
            os.system('rm -f hadd_jobs/hadd_command_%s.log' % basename)
            os.system('rm -f hadd_jobs/hadd_command_%s.jdl'% basename)
            condor_file = open('hadd_jobs/hadd_command_%s.jdl' % basename, 'w')
            condor_file.write('universe = vanilla\n')
            condor_file.write('Executable = hadd_jobs/hadd_command_%s.sh\n'% basename)
            condor_file.write('Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )\n')
            condor_file.write('request_disk = 3000000\n')
            #condor_file.write('request_disk = 4000000\n')
            condor_file.write('request_memory = 10000\n')
            #condor_file.write('request_memory = 10000\n')
            condor_file.write('Should_Transfer_Files = YES\n')
            condor_file.write('WhenToTransferOutput = ON_EXIT\n')
            condor_file.write('Transfer_Input_Files = %s\n'%filesToTransfer)
            condor_file.write('use_x509userproxy = true\n')
            condor_file.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
            condor_file.write('Output = %s.stdout\n' % os.path.abspath(condor_file.name))
            condor_file.write('Error = %s.stderr\n' % os.path.abspath(condor_file.name))
            condor_file.write('Log = %s.log\n' % os.path.abspath(condor_file.name))
            condor_file.write('Queue 1\n')
            condor_file.close()
            os.system('chmod +x %s'% os.path.abspath(condor_file.name))
            exec_me('condor_submit %s'%(os.path.abspath(condor_file.name)),options.dryRun)
                
        with open('hadd_jobs/bad_files_%s.txt'%basename,'w') as f:
            for badFile in badFiles:
                f.write(badFile+'\n')


def getFilesRecursively(dir,searchstring,additionalstring = None, skipString = None):
    
    thesearchstring = searchstring
    theadditionalstring = None
    if not additionalstring == None: 
        theadditionalstring = additionalstring

    cfiles = []
    badfiles = []
    files = []
    print 'ls %s/%s'%(dir,thesearchstring)
    os.system('ls %s/%s > tmp.txt'%(dir,thesearchstring))
    with open("tmp.txt", 'r') as mylist:
        files = [(myfile.replace('\n', ''), True) for myfile in mylist.readlines()]

    nfiles = len(files)
    for ifile, fi in enumerate(files):
        if ifile%100==0:
            print '%i/%i files checked in %s'%(ifile,nfiles,dir+'/'+thesearchstring)
        if 'runPu' not in options.executable:
            try:
                filename = 'root://cmseos.fnal.gov//%s/%s/%s'%(dir,thesearchstring,fi[0])
                f = ROOT.TFile.Open(filename)
                if f.IsZombie():
                    print 'file is zombie'
                    f.Close()
                    badfiles.append(filename)
                    continue
                elif not f.Get('Events'):
                    print 'tree is false'
                    f.Close()
                    badfiles.append(filename)
                    continue
                elif not f.Get('Events').InheritsFrom('TTree'):
                    print 'tree is not a tree'
                    f.Close()
                    badfiles.append(filename)
                    continue
                else:
                    f.Close()
                    cfiles.append(filename)
            except:
                print 'could not open file or tree'
                badfiles.append(filename)
                continue
        else:
            filename = 'root://cmseos.fnal.gov//%s/%s/%s'%(dir,thesearchstring,fi[0])
            cfiles.append(filename)
                
    return cfiles, badfiles



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('--train', action='store_true', dest='train', default=False, help='train')
    parser.add_option("--lumi", dest="lumi", default = 30,type='float',help="luminosity", metavar="lumi")
    parser.add_option('-i','--idir', dest='idir', default = 'data/',help='directory with bacon bits', metavar='idir')
    parser.add_option('-j','--job', dest='job', default =-1,type='int',help='just do part of it', metavar='idir')
    #parser.add_option('-o','--odir', dest='odir', default = 'skim/',help='directory to write hadded bits', metavar='odir')
    parser.add_option('-s','--sample',dest="sample", default="All",type='string',
                      #choices=['All','Hbb','QCD','JetHT','SingleMuon','DMSpin0','TT','DY','W','Diboson','Triboson','SingleTop','VectorDiJet1Jet','VectorDiJet1Gamma','MC','Data'],
                      help="samples to produces")
    parser.add_option('-e','--executable',dest="executable", default="runZprime", help = "executable name")
    parser.add_option("--tagddt", type=str, default='', help="tagddt")
    parser.add_option("--ddt", type=str, default='GridOutput_v13.root', help="ddt")
    parser.add_option("--iddt", type=str, default='Rho2D', help="iddt")
    parser.add_option('--jet', dest='jet', default='AK8', help='jet type')
    parser.add_option('--loop',dest="loop",default=False,action='store_true', help='loop over dir to get list of files')
    parser.add_option('--haddlarge',dest="haddlarge",default=False,action='store_true', help='hadd large files')
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                      help="Just print out commands to run")
    parser.add_option('--just-hadd',dest="justHadd",default=False,action='store_true',
                      help="Just run hadd (two-step)")
    parser.add_option('--just-sklim',dest="justSklim",default=False,action='store_true',
                      help="Just run sklim step")
    parser.add_option('--just-norm',dest="justNorm",default=False,action='store_true',
                      help="Just run norm step")
    parser.add_option('--just-DDT',dest="justDDT",default=False,action='store_true',
                      help="Run DDT skimmer")
    parser.add_option('--just-N2',dest="justN2",default=False,action='store_true',
                      help="Run N2 skimmer")
    parser.add_option('--just-Wtag',dest="justWtag",default=False,action='store_true',
                      help="Run Wtag skimmer")
    (options, args) = parser.parse_args()

    
    if options.justHadd:
        justHadd(options,args)
    else:
        global haddDir,sklimDir,normDir;
        haddDir = "hadd"
        sklimDir = "sklim"
        normDir = "norm"
        if options.justDDT:
            haddDir = "hadd"
            sklimDir +='DDT'
            normDir +='DDT'
        if options.justN2:
            haddDir = "hadd"
            sklimDir +='N2'
            normDir +='N2'
        if options.justWtag:
            haddDir = "hadd"
            sklimDir +='Wtag%s%s'%(options.jet,options.tagddt)

        if '2016' in options.sample:
            haddDir+='2016'
            sklimDir+='2016'
            normDir+='2016'
        main(options,args)
    
