#!/usr/bin/env python

import sys, commands, os, fnmatch
from optparse import OptionParser

EOS = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'

def exec_me(command, dryRun=False):
    print command
    if not dryRun:
        os.system(command)
        
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")    
    parser.add_option("--monitor",default='',help="Monitor mode (sub/resub/check directory of jobs)")
    parser.add_option('-s','--sample',dest="sample", default="All",
                      choices=['All','Hbb','QCD','JetHT','DMSpin0','TT','DY','W','Diboson','Triboson','SingleTop','VectorDiJet1Jet','VectorDiJet1Gamma','MC','Data'],
                      help="samples to produces")
    
    (options,args) = parser.parse_args()

    monitorOption = ''
    if options.monitor is not '':
        monitorOption = '--monitor %s'%options.monitor
    
    json = "$PWD/../data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"
    xsec = 1

    eosOutDir = '/eos/cms/store/group/phys_exotica/dijet/dazsle'

    
    optionsDataMc = {
        'mc': "Output.root --passSumEntries 5:Events  -a 2:mc -a 3:none  -n 7000 -q 2nw4cores",
        'data': "Output.root -a 5:1  -a 2:data -a 3:%s -n 7000 -q 1nh"%(json)
    }
        
    analysisDir = "zprimebits-v11.051"
    executable = "runZprime"
    
    samplesDict = {}
    
    
    samplesDict['JetHT'] = {
        'JetHTRun2016B_PromptReco_v2': 'data',
        'JetHTRun2016C_PromptReco_v2': 'data',
        'JetHTRun2016D_PromptReco_v2': 'data',
        'JetHTRun2016E_PromptReco_v2': 'data',
        'JetHTRun2016F_PromptReco_v1': 'data',
        'JetHTRun2016G_PromptReco_v1': 'data'
        }
    samplesDict['Hbb'] = {
        'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8': 'mc',  
        'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8_ext': 'mc',        
        'GluGluHToBB_M125_13TeV_powheg_herwigpp': 'mc',
        'GluGluHToBB_M125_13TeV_powheg_pythia8': 'mc',
        'VBFHToBB_M125_13TeV_amcatnlo_pythia8': 'mc',
        'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix': 'mc',
        'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_ext': 'mc',
        'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'mc',
        'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'mc',
        'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'mc',
        'ttHTobb_M125_13TeV_powheg_pythia8': 'mc',
        'ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV_powheg_pythia8': 'mc'
        }
    samplesDict['QCD'] = {
        'QCD_HT100to200_13TeV':'mc',
        'QCD_HT200to300_13TeV_ext': 'mc',
        'QCD_HT300to500_13TeV_ext': 'mc',
        'QCD_HT500to700_13TeV_ext': 'mc',
        'QCD_HT700to1000_13TeV_ext':'mc',
        'QCD_HT1000to1500_13TeV_ext':'mc',
        'QCD_HT1500to2000_13TeV_ext':'mc',
        'QCD_HT2000toInf_13TeV_ext':'mc'
        }
    samplesDict['SingleTop'] = {
        'ST_t-channel_antitop_4f_inclusiveDecays_13TeV_powheg': 'mc',
        'ST_t-channel_top_4f_inclusiveDecays_13TeV_powheg': 'mc',
        'ST_tW_antitop_5f_inclusiveDecays_13TeV': 'mc',
        'ST_tW_top_5f_inclusiveDecays_13TeV': 'mc'
        }
    samplesDict['W'] = {
        'WJetsToQQ_HT180_13TeV': 'mc',
        'WJetsToQQ_HT_600ToInf_13TeV': 'mc'
        }
    samplesDict['DY'] = {
        'DYJetsToQQ_HT180_13TeV': 'mc',
        'ZJetsToQQ_HT600toInf_13TeV_madgraph': 'mc',
        }
    samplesDict['TT'] = {
        'TTJets_13TeV':'mc',
        }
    samplesDict['Diboson'] = {
        'WWTo4Q_13TeV_amcatnlo': 'mc',
        'WWTo4Q_13TeV_powheg': 'mc',
        'ZZTo4Q_13TeV_amcatnlo':'mc',
        'WZ_13TeV': 'mc'
        }
    samplesDict['Triboson'] = {
        'TTWJetsToQQ_13TeV': 'mc',
        'TTGJets_13TeV': 'mc',
        'TTZToQQ_13TeV': 'mc',
        }
    samplesDict['VectorDiJet1Jet'] = {
        'VectorDiJet1Jet_M50':'mc',
        'VectorDiJet1Jet_M75':'mc',
        'VectorDiJet1Jet_M100':'mc',
        'VectorDiJet1Jet_M125':'mc',
        'VectorDiJet1Jet_M150':'mc',
        'VectorDiJet1Jet_M200':'mc',
        'VectorDiJet1Jet_M250':'mc',
        'VectorDiJet1Jet_M300':'mc',
        'VectorDiJet1Jet_M400':'mc',
        'VectorDiJet1Jet_M500':'mc'
        }
    samplesDict['VectorDiJet1Gamma'] = {
        'VectorDiJet1Gamma_20_1_800_v2':'mc',
        'VectorDiJet1Gamma_50_1_800_v2':'mc',
        'VectorDiJet1Gamma_75_1_800_v2':'mc',
        'VectorDiJet1Gamma_100_1_800_v2':'mc',
        'VectorDiJet1Gamma_125_1_800_v2':'mc',
        'VectorDiJet1Gamma_150_1_800_v2':'mc',
        'VectorDiJet1Gamma_200_1_800_v2':'mc',
        'VectorDiJet1Gamma_300_1_800_v2':'mc',
        'VectorDiJet1Gamma_400_1_800_v2':'mc'
    }
    samplesDict['DMSpin0'] = {
        'DMSpin0_ggPhibb1j_50':'mc',
        'DMSpin0_ggPhibb1j_75':'mc',
        'DMSpin0_ggPhibb1j_100':'mc',
        'DMSpin0_ggPhibb1j_125':'mc',
        'DMSpin0_ggPhibb1j_150':'mc',
        'DMSpin0_ggPhibb1j_200':'mc',
        'DMSpin0_ggPhibb1j_250':'mc',
        'DMSpin0_ggPhibb1j_300':'mc',
        'DMSpin0_ggPhibb1j_400':'mc',
        'DMSpin0_ggPhibb1j_500':'mc'
        }
    
    samplesDict['MC'] = dict(samplesDict['Hbb'].items() +
                             samplesDict['QCD'].items() +
                             samplesDict['SingleTop'].items() +
                             samplesDict['W'].items() +
                             samplesDict['DY'].items() +
                             samplesDict['TT'].items() +
                             samplesDict['Diboson'].items() +
                             samplesDict['Triboson'].items() +
                             samplesDict['VectorDiJet1Jet'].items() +
                             samplesDict['VectorDiJet1Gamma'].items() +
                             samplesDict['DMSpin0'].items())
    
    samplesDict['Data'] = dict(samplesDict['JetHT'].items())
                              
    samplesDict['All'] = dict(samplesDict['MC'].items() + samplesDict['Data'].items())

    samples = samplesDict[options.sample]

    exec_me('%s mkdir -p %s/%s'%(EOS,eosOutDir,analysisDir))  
    for label, isMc in samples.iteritems():
        exec_me('%s mkdir -p %s/%s/%s'%(EOS,eosOutDir,analysisDir,label))  
        exec_me("python baconBatch.py %s %s -a 4:%f --list 1:../lists/production11/%s.txt --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s  %s"%(executable,optionsDataMc[isMc],xsec,label,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
