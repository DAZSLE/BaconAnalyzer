#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
from optparse import OptionParser
from submitZprime import samplesDict, exec_me

import ROOT
normDict = {'DYJetsToQQ_HT180_13TeV': 'DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8',
            'WJetsToQQ_HT180_13TeV': 'WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8',
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
            'TTJets_13TeV': 'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'TT_powheg': 'TT_TuneCUETP8M1_13TeV-powheg-pythia8',
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
            'WWTo4Q_13TeV_amcatnlo': 'WWTo4Q_4f_13TeV_amcatnloFXFX_madspin_pythia8',
            'WWTo4Q_13TeV_powheg': 'WWTo4Q_13TeV_powheg',
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
            'DYJetsToLL_M_50_13TeV_ext': 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8': 'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8',
            'WW_13TeV_pythia8': 'WW_TuneCUETP8M1_13TeV-pythia8',
            'ZZ_13TeV_pythia8': 'ZZ_TuneCUETP8M1_13TeV-pythia8'
            }

def justHadd(options,args):    
    DataDir = options.idir
    OutDir = options.idir
    
    samples = samplesDict[options.sample]
    
    EOS = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'
    postfix = ''
    exec_me('mkdir -p $PWD/hadd_jobs/',options.dryRun)
    exec_me('%s mkdir -p /%s/hadd'%(EOS,OutDir),options.dryRun)
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
        for i in range(0,len(filesToConvert)/500+1):
            if not os.path.isfile(OutDir+'/hadd/'+basename.replace('.root','_%i.root'%i)):
                haddOutExistsList.append(False)
                haddOutList.append(OutDir+'/hadd/'+basename.replace('.root','_%i.root'%i))
                haddCommand = '#!/bin/bash\n'
                haddCommand += 'pwd\n'
                haddCommand += 'cd %s\n'%cwd
                haddCommand += 'pwd\n'
                haddCommand += 'eval `scramv1 runtime -sh`\n'
                haddCommand += 'cd -\n'
                haddCommand += 'pwd\n'
                haddCommand += 'mkdir -p $PWD/hadd\n'        
                haddCommand += 'hadd -f hadd/%s %s\n'%(basename.replace('.root','_%i.root'%i),(' '.join(filesToConvert[i*500:(i+1)*500])).replace('eos','root://eoscms.cern.ch//eos'))
                haddCommand += '%s cp $PWD/hadd/%s /%s/hadd/%s\n'%(EOS,basename.replace('.root','_%i.root'%i),OutDir,basename.replace('.root','_%i.root'%i))
                haddCommand += 'rm -r $PWD/hadd\n'
                with open('hadd_jobs/hadd_command_%s.sh'%(basename.replace('.root','_%i.root'%i)),'w') as f:
                    f.write(haddCommand)
                exec_me('bsub -q 8nh -o $PWD/hadd_jobs/hadd_command_%s.log source $PWD/hadd_jobs/hadd_command_%s.sh'%(basename.replace('.root','_%i.root'%i),basename.replace('.root','_%i.root'%i)),options.dryRun)
            else:
                haddOutExistsList.append(True)
                haddOutList.append(OutDir+'/hadd/'+basename.replace('.root','_%i.root'%i))
                
        haddAll = all(haddOutExistsList)

        print 'haddAll', haddAll
        print 'haddOutExistsList', haddOutList
        print 'haddOutList', haddOutList
        
        if haddAll:
            haddCommand = '#!/bin/bash\n'
            haddCommand += 'pwd\n'
            haddCommand += 'cd %s\n'%cwd
            haddCommand += 'pwd\n'
            haddCommand += 'eval `scramv1 runtime -sh`\n'
            haddCommand += 'cd -\n'
            haddCommand += 'pwd\n'
            haddCommand += 'mkdir -p $PWD/hadd\n'        
            haddCommand += 'hadd -f hadd/%s %s\n'%(basename,(' '.join(haddOutList)).replace('eos','root://eoscms.cern.ch//eos'))
            haddCommand += '%s cp $PWD/hadd/%s /%s/hadd/%s\n'%(EOS,basename,OutDir,basename)
            haddCommand += 'rm -r $PWD/hadd\n'
            with open('hadd_jobs/hadd_command_%s.sh'%(basename),'w') as f:
                f.write(haddCommand)
            exec_me('bsub -q 8nh -o $PWD/hadd_jobs/hadd_command_%s.log source $PWD/hadd_jobs/hadd_command_%s.sh'%(basename,basename),options.dryRun)
            

    
def main(options,args):

    DataDir = options.idir
    OutDir = options.idir
    
    samples = samplesDict[options.sample]
    
    EOS = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'
    postfix = ''
    exec_me('mkdir -p $PWD/hadd_jobs/',options.dryRun)
    exec_me('%s mkdir -p /%s/hadd'%(EOS,OutDir),options.dryRun)
    exec_me('%s mkdir -p /%s/sklim'%(EOS,OutDir),options.dryRun)
    exec_me('%s mkdir -p /%s/norm'%(EOS,OutDir),options.dryRun)
    for label, isMc in samples.iteritems():
        basename = label + '.root'
        if options.job>-1:
            basename = label + '_%i.root'%options.job
            
        haddOn = True
        sklimOn = True
        normOn = True
        if os.path.isfile(OutDir+'/hadd/'+basename):
            haddOn = False
        if os.path.isfile(OutDir+'/sklim/'+basename):
            sklimOn = False
        if os.path.isfile(OutDir+'/norm/'+basename.replace('.root','_1000pb_weighted.root')):
            normOn = False

        if options.justNorm:            
            haddOn = False
            sklimOn = False
            normOn = True
        elif options.justSklim:  
            haddOn = False
            sklimOn = True
            normOn = False
            
            
        
        filesToConvert = []
        badFiles = []
        if haddOn:
            filesToConvert, badFiles = getFilesRecursively(DataDir,label+'/',None,None)
        print "files To Convert = ",filesToConvert
        print "bad files = ", badFiles
        cwd = os.getcwd()
        haddCommand = '#!/bin/bash\n'
        haddCommand += 'pwd\n'
        haddCommand += 'cd %s\n'%cwd
        haddCommand += 'pwd\n'
        haddCommand += 'eval `scramv1 runtime -sh`\n'
        haddCommand += 'cd -\n'
        haddCommand += 'pwd\n'
        haddCommand += 'mkdir -p $PWD/hadd\n'        
        for i in range(0,len(filesToConvert)/500+1):         
            if haddOn:
                haddCommand += 'hadd -f hadd/%s %s\n'%(basename.replace('.root','_%i.root'%i),(' '.join(filesToConvert[i*500:(i+1)*500])).replace('eos','root://eoscms.cern.ch//eos'))
        if haddOn:
            haddCommand += 'hadd -f $PWD/hadd/%s $PWD/hadd/%s\n'%(basename,basename.replace('.root','_*.root'))
            haddCommand += 'rm $PWD/hadd/%s\n'%(basename.replace('.root','_*.root'))     
            haddCommand += '%s cp $PWD/hadd/%s /%s/hadd/%s\n'%(EOS,basename,OutDir,basename)
        if sklimOn:
            haddCommand += 'mkdir -p $PWD/sklim\n'
            haddCommand += 'export PYTHONPATH=${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/:${PYTHONPATH}\n'
            if not haddOn:
                haddCommand += '%s cp /%s/hadd/%s $PWD/hadd/%s\n'%(EOS,OutDir,basename,basename)                
            haddCommand += 'python ${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/skimmer.py -i $PWD/hadd/ -o $PWD/sklim/ -s %s\n'%(basename.replace('.root',''))
            haddCommand += '%s cp $PWD/sklim/%s /%s/sklim/%s\n'%(EOS,basename,OutDir,basename)        
        if isMc=='mc' and normOn:            
            if not sklimOn:
                haddCommand += '%s cp /%s/sklim/%s $PWD/sklim/%s\n'%(EOS,OutDir,basename,basename)                
            haddCommand += 'echo "%s\t${PWD}/sklim/%s" > normlist.txt\n'%(normDict[basename.replace('.root','')],basename)
            haddCommand += 'NormalizeNtuple normlist.txt\n'
            haddCommand += '%s cp $PWD/sklim/%s /%s/norm/%s\n'%(EOS,basename.replace('.root','_1000pb_weighted.root'),OutDir,basename.replace('.root','_1000pb_weighted.root'))
        haddCommand += 'rm -r $PWD/hadd\n'
        haddCommand += 'rm -r $PWD/sklim\n'
            
            
        with open('hadd_jobs/hadd_command_%s.sh'%(basename),'w') as f:
            f.write(haddCommand)

            if haddOn or sklimOn or (isMc=='mc' and normOn):
                exec_me('bsub -q 8nh -o $PWD/hadd_jobs/hadd_command_%s.log source $PWD/hadd_jobs/hadd_command_%s.sh'%(basename,basename),options.dryRun)

        with open('hadd_jobs/bad_files_%s.txt'%basename,'w') as f:
            for badFile in badFiles:
                f.write(badFile+'\n')


def getFilesRecursively(dir,searchstring,additionalstring = None, skipString = None):
    
    # thesearchstring = "_"+searchstring+"_"
    thesearchstring = searchstring

    theadditionalstring = None
    if not additionalstring == None: 
        theadditionalstring = additionalstring

    cfiles = []
    badfiles = []
    for root, dirs, files in os.walk(dir+'/'+thesearchstring):
        nfiles = len(files)
        for ifile, file in enumerate(files):
            
            if ifile%100==0:
                print '%i/%i files checked in %s'%(ifile,nfiles,dir+'/'+thesearchstring)
            try:
                #f = ROOT.TFile.Open((os.path.join(root, file)).replace('eos','root://eoscms.cern.ch//eos'))
                f = ROOT.TFile.Open((os.path.join(root, file)))
                if f.IsZombie():
                    print 'file is zombie'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                elif not f.Get('Events'):
                    print 'tree is false'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                elif not f.Get('Events').InheritsFrom('TTree'):
                    print 'tree is not a tree'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                else:
                    f.Close()
                    cfiles.append(os.path.join(root, file))                    
            except:
                print 'could not open file or tree'
                badfiles.append(os.path.join(root, file))                    
                continue
                
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
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")
    parser.add_option('--just-hadd',dest="justHadd",default=False,action='store_true',
                  help="Just run hadd (two-step)")
    parser.add_option('--just-sklim',dest="justSklim",default=False,action='store_true',
                  help="Just run sklim step")
    parser.add_option('--just-norm',dest="justNorm",default=False,action='store_true',
                  help="Just run norm step")
    (options, args) = parser.parse_args()

    
    if options.justHadd:
        justHadd(options,args)
    else:
        main(options,args)
    
