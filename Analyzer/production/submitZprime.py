#!/usr/bin/env python

import sys, commands, os, fnmatch
from optparse import OptionParser

samplesDict = {}
samplesDict['JetHTrereco'] = {
    'JetHTRun2017B_17Nov2017_v1_noPF_old': 'prompt17',
    'JetHTRun2017C_17Nov2017_v1_noPF_old': 'prompt17',
    'JetHTRun2017D_17Nov2017_v1_noPF_old': 'prompt17',
    'JetHTRun2017E_17Nov2017_v1_noPF_old': 'prompt17',
    'JetHTRun2017F_17Nov2017_v1_noPF_old': 'prompt17',
}
samplesDict['JetHT'] = {
    # 'JetHTRun2017A_12Sep2017_v1': 'prompt17',
    # 'JetHTRun2017A_PromptReco_v1': 'prompt17',
    # 'JetHTRun2017A_PromptReco_v2': 'prompt17',
    # 'JetHTRun2017A_PromptReco_v3': 'prompt17',
    #'JetHTRun2017B_12Sep2017_v1': 'prompt17',
    #'JetHTRun2017B_22Jun2017_v1': 'prompt17',
    #'JetHTRun2017B_23Jun2017_v1': 'prompt17',
    #'JetHTRun2017B_12Sep2017_v1_noPF' : 'prompt17',
    #'JetHTRun2017B_22Jun2017_v1_noPF' : 'prompt17',
    #'JetHTRun2017B_23Jun2017_v1_noPF' : 'prompt17',
    'JetHTRun2017B_PromptReco_v1_noPF' : 'prompt17',
    'JetHTRun2017B_PromptReco_v2_noPF' : 'prompt17',
    #'JetHTRun2017C_12Sep2017_v1_noPF' : 'prompt17',
    'JetHTRun2017C_PromptReco_v1_noPF' : 'prompt17',
    'JetHTRun2017C_PromptReco_v2_noPF' : 'prompt17',
    'JetHTRun2017C_PromptReco_v3_noPF' : 'prompt17',
    'JetHTRun2017D_PromptReco_v1_noPF' : 'prompt17',
    'JetHTRun2017E_PromptReco_v1_noPF' : 'prompt17',
    'JetHTRun2017F_PromptReco_v1_noPF' : 'prompt17',
    #'JetHTRun2017B_PromptReco_v1': 'prompt17',
    #'JetHTRun2017B_PromptReco_v2_noPF': 'prompt17',
    #'JetHTRun2017C_12Sep2017_v1': 'prompt17',
    #'JetHTRun2017C_PromptReco_v1': 'prompt17',
    #'JetHTRun2017C_PromptReco_v2': 'prompt17',
    #'JetHTRun2017C_PromptReco_v3': 'prompt17',
    #'JetHTRun2017D_PromptReco_v1': 'prompt17',
    #'JetHTRun2017E_PromptReco_v1_noPF': 'prompt17',
    #'JetHTRun2017F_PromptReco_v1_noPF': 'prompt17',
    # 'JetHTRun2016B_03Feb2017_ver1_v1_v3': 'rereco',
    # 'JetHTRun2016B_03Feb2017_ver2_v2_v3': 'rereco',
    # 'JetHTRun2016C_03Feb2017_v1_v3': 'rereco',
    # 'JetHTRun2016D_03Feb2017_v1_v3': 'rereco',
    # 'JetHTRun2016E_03Feb2017_v1_v3': 'rereco',
    # 'JetHTRun2016F_03Feb2017_v1_v3': 'rereco',
    # 'JetHTRun2016G_03Feb2017_v1_v3': 'rereco',
    # 'JetHTRun2016H_03Feb2017_ver2_v1_v3': 'data',
    # 'JetHTRun2016H_03Feb2017_ver3_v1_v3': 'data'
}
samplesDict['SingleMuonrereco'] = {
    'SingleMuonRun2017B_17Nov2017_v1_noPF': 'prompt17',
    #'SingleMuonRun2017C_17Nov2017_v1_noPF': 'prompt17',
    'SingleMuonRun2017D_17Nov2017_v1_noPF': 'prompt17',
    'SingleMuonRun2017E_17Nov2017_v1_noPF_304119': 'prompt17',
    'SingleMuonRun2017E_17Nov2017_v1_noPF_304507': 'prompt17',
    'SingleMuonRun2017E_17Nov2017_v1_noPF_305185': 'prompt17',
    'SingleMuonRun2017E_17Nov2017_v1_noPF_306029': 'prompt17',
    'SingleMuonRun2017E_17Nov2017_v1_noPF_306126': 'prompt17',
    'SingleMuonRun2017F_17Nov2017_v1_noPF_305185': 'prompt17',
    'SingleMuonRun2017F_17Nov2017_v1_noPF_305248': 'prompt17',
    'SingleMuonRun2017F_17Nov2017_v1_noPF_305364': 'prompt17',
    'SingleMuonRun2017F_17Nov2017_v1_noPF_305636': 'prompt17',
    'SingleMuonRun2017F_17Nov2017_v1_noPF_306029': 'prompt17',
    'SingleMuonRun2017F_17Nov2017_v1_noPF_306126': 'prompt17',
}
samplesDict['SingleMuon'] = {
    'SingleMuonRun2017B_PromptReco_v2_noPF' : 'prompt17',
    'SingleMuonRun2017C_PromptReco_v1_noPF' : 'prompt17',
    'SingleMuonRun2017C_PromptReco_v2_noPF' : 'prompt17',
    'SingleMuonRun2017C_PromptReco_v3_noPF' : 'prompt17',
    'SingleMuonRun2017D_PromptReco_v1_noPF' : 'prompt17',
    #'SingleMuonRun2017E_PromptReco_v1_noPF' : 'prompt17',
    #'SingleMuonRun2017F_PromptReco_v1_noPF' : 'prompt17',
    'SingleMuonRun2017G_PromptReco_v1_noPF' : 'prompt17',
    #'SingleMuonRun2017H_PromptReco_v1_noPF' : 'prompt17',
    # 'SingleMuonRun2016B_03Feb2017_ver1_v1': 'rereco',
    # 'SingleMuonRun2016B_03Feb2017_ver2_v2': 'rereco',
    # 'SingleMuonRun2016C_03Feb2017_v1': 'rereco',
    # 'SingleMuonRun2016D_03Feb2017_v1': 'rereco',
    # 'SingleMuonRun2016E_03Feb2017_v1': 'rereco',
    # 'SingleMuonRun2016F_03Feb2017_v1': 'rereco',
    # 'SingleMuonRun2016G_03Feb2017_v1': 'rereco',
    # 'SingleMuonRun2016H_03Feb2017_ver2_v1': 'data',
    # 'SingleMuonRun2016H_03Feb2017_ver3_v1': 'data'
}
samplesDict['BulkGrav'] = {
   'BulkGravToWW_narrow_M_4000_13TeV_madgraph': 'mc',
   'BulkGravToWW_narrow_M_4000_13TeV_madgraph_herwigpp': 'mc', 
   'BulkGravToWW_narrow_M_600_13TeV_madgraph': 'mc',
   'BulkGravToWW_narrow_M_600_13TeV_madgraph_herwigpp': 'mc',	
   'BulkGravToWW_narrow_M_1000_13TeV_madgraph': 'mc', 
   'BulkGravToWW_narrow_M_1000_13TeV_madgraph_herwigpp': 'mc', 
   'BulkGravToWW_narrow_M_2000_13TeV_madgraph': 'mc', 
   'BulkGravToWW_narrow_M_2000_13TeV_madgraph_herwigpp': 'mc', 
   'BulkGravToWW_narrow_M_3000_13TeV_madgraph': 'mc', 
   'BulkGravToWW_narrow_M_3000_13TeV_madgraph_herwigpp': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_1000_13TeV_madgraph': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_1000_13TeV_madgraph_herwig_ext': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_2000_13TeV_madgraph': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_2000_13TeV_madgraph_herwig_ext': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_3000_13TeV_madgraph': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_3000_13TeV_madgraph_herwig': 'mc'

}
samplesDict['Hbb'] = {
    #'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8': 'mc', 
    'GluGluHToBB_M125_13TeV_powheg_pythia8': 'mc', 
    'GluGluHToBB_M125_13TeV_powheg_pythia8_ext': 'mc',
    #'VBFHToBB_M125_13TeV_amcatnlo_pythia8': 'mc', 
    #'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix': 'mc', 
    #'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_ext': 'mc',
    #'WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8': 'mc', 
    #'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'mc', 
    #'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8': 'mc', 
    #'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'mc', 
    #'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8': 'mc', 
    #'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_ext': 'mc', 
    #'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'mc', 
    #'bbHToBB_M_125_4FS_yb2_13TeV_amcatnlo': 'mc', 
    #'bbHToBB_M_125_4FS_ybyt_13TeV_amcatnlo': 'mc', 
    #'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8': 'mc', 
    #'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'mc', 
    #'ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8':'mc',
    #'ttHTobb_M125_13TeV_powheg_pythia8': 'mc', 
   # 'ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV_powheg_pythia8': 'mc'
    }
samplesDict['QCD'] = {       
    #'QCD_HT1000to1500_13TeV': 'mc',
    'QCD_HT1000to1500_13TeV_ext': 'mc',
   # 'QCD_HT1000to1500_13TeV_all': 'mc',
    'QCD_HT100to200_13TeV': 'mc',
  #  'QCD_HT1500to2000_13TeV': 'mc',
    'QCD_HT1500to2000_13TeV_ext': 'mc',
   # 'QCD_HT1500to2000_13TeV_all': 'mc',
  #  'QCD_HT2000toInf_13TeV': 'mc',
    'QCD_HT2000toInf_13TeV_ext': 'mc',
   # 'QCD_HT2000toInf_13TeV_all': 'mc',
    'QCD_HT200to300_13TeV': 'mc',
  #  'QCD_HT200to300_13TeV_ext': 'mc',
  #  'QCD_HT200to300_13TeV_all': 'mc',
  #  'QCD_HT300to500_13TeV': 'mc',
    'QCD_HT300to500_13TeV_ext': 'mc',
 #   'QCD_HT300to500_13TeV_all': 'mc',
  #  'QCD_HT500to700_13TeV': 'mc',
    'QCD_HT500to700_13TeV_ext': 'mc',
 #   'QCD_HT500to700_13TeV_all': 'mc',
  #  'QCD_HT50to100_13TeV': 'mc',
  #  'QCD_HT700to1000_13TeV': 'mc',
    'QCD_HT700to1000_13TeV_ext': 'mc',
 #   'QCD_HT700to1000_13TeV_all': 'mc',
    }
samplesDict['SingleTop'] = {
    'ST_s_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1': 'mc',
    #'ST_t_channel_antitop_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1': 'mc',
    'ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin': 'mc',
    'ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin': 'mc',
    'ST_tW_antitop_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4': 'mc',
    'ST_tW_top_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4': 'mc'
    }
samplesDict['Diboson'] = {
    'WWTo4Q_13TeV_amcatnlo': 'mc',
    'WWTo4Q_13TeV_powheg': 'mc',
    'WZ_13TeV_pythia8':'mc',
    'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8':'mc',
    'WW_13TeV_pythia8':'mc',
    'ZZ_13TeV_pythia8':'mc' 
    }
samplesDict['Triboson'] = {
    #'TTWJetsToQQ_13TeV': 'mc',
    #'TTGJets_13TeV': 'mc',
    #'TTGJets_13TeV_ext1': 'mc',
    #'TTZToQQ_13TeV': 'mc',
    #'TTToSemilepton_powheg': 'mc',
    #'WWW_4F_13TeV_amcatnlo_pythia8': 'mc',
    #'WWZ_13TeV_amcatnlo_pythia8': 'mc',
    #'ZZZ_13TeV_amcatnlo_pythia8': 'mc'
    }
samplesDict['VectorDiJet1Jet'] = {
    #'VectorDiJet1Jet_100_madgraph_2016': 'mc',
    #'VectorDiJet1Jet_125_madgraph_2017_noPF': 'mc',
    #'VectorDiJet1Jet_200_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_75_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_50_madgraph_2017_noPF' : 'mc',
    'VectorDiJet1Jet_100_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_115_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_150_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_175_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_250_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_300_madgraph_2017_noPF' : 'mc',
    # 'VectorDiJet1Jet_100_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_150_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_200_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_25_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_300_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_400_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_500_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_50_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_600_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_800_13TeV_madgraph': 'mc', 
    # 'VectorDiJet1Jet_1000_13TeV_madgraph': 'mc',
    # 'VectorDiJet1Jet_125_13TeV_madgraph': 'mc',
    # 'VectorDiJet1Jet_75_13TeV_madgraph': 'mc'
    }
samplesDict['VectorDiJet1Gamma'] = {
    #'VectorDiJet1Gamma_1000_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_100_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_10_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_125_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_200_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_25_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_300_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_400_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_500_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_50_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_600_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_75_13TeV_madgraph': 'mc', 
    #'VectorDiJet1Gamma_800_13TeV_madgraph': 'mc', 
    }
samplesDict['W'] = {
    'WJetsToQQ_HT180_13TeV': 'mc',
    #'WJetsToQQ_HT_600ToInf_13TeV': 'mc',
    #'WJetsToLNu_HT_70To100_13TeV': 'mc',
    #'WJetsToLNu_HT_100To200_13TeV': 'mc',
    #'WJetsToLNu_HT_100To200_13TeV_ext1': 'mc',
    #'WJetsToLNu_HT_100To200_13TeV_ext2': 'mc',
    #'WJetsToLNu_HT_200To400_13TeV': 'mc',
    #'WJetsToLNu_HT_200To400_13TeV_ext1': 'mc',
    #'WJetsToLNu_HT_200To400_13TeV_ext2': 'mc',
    #'WJetsToLNu_HT_400To600_13TeV': 'mc',
    #'WJetsToLNu_HT_400To600_13TeV_ext1': 'mc',
    #'WJetsToLNu_HT_600To800_13TeV': 'mc',
    #'WJetsToLNu_HT_600To800_13TeV_ext1': 'mc',
    #'WJetsToLNu_HT_800To1200_13TeV': 'mc',
    #'WJetsToLNu_HT_800To1200_13TeV_ext1': 'mc',
    #'WJetsToLNu_HT_1200To2500_13TeV': 'mc',
    #'WJetsToLNu_HT_1200To2500_13TeV_ext1': 'mc',
    #'WJetsToLNu_HT_2500ToInf_13TeV': 'mc',
    #'WJetsToLNu_HT_2500ToInf_13TeV_ext1': 'mc',
    }
samplesDict['DY'] = {
    #'DYJetsToLL_M_50_13TeV_ext': 'mc',
    'DYJetsToQQ_HT180_13TeV': 'mc',
    #'ZJetsToQQ_HT600toInf_13TeV_madgraph': 'mc',
    }
samplesDict['TT'] = {
    'TT_13TeV_powheg_pythia8_ext':'mc',
    #'TT_powheg':'mc',
    #'TT_TuneEE5C_13TeV_powheg_herwigpp':'mc'
    }
samplesDict['DMSpin0'] = {
    'Spin0_ggPhibb1j_g1_1000_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_1000_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_100_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_100_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_10_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_10_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_125_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_200_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_20_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_300_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_300_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_350_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_350_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_400_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_400_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_500_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_500_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_50_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_50_Scalar': 'mc',
    'Spin0_ggPhibb1j_g1_750_PseudoScalar': 'mc',
    'Spin0_ggPhibb1j_g1_750_Scalar': 'mc', 
    }

samplesDict['MC'] = dict(samplesDict['Hbb'].items() +
                         samplesDict['QCD'].items() +
                         samplesDict['SingleTop'].items() +
                         samplesDict['W'].items() +
                         samplesDict['DY'].items() +
                         samplesDict['TT'].items() +
                         samplesDict['Diboson'].items() +
                        #samplesDict['Triboson'].items() +
                         samplesDict['VectorDiJet1Jet'].items() +
                       # samplesDict['VectorDiJet1Gamma'].items() +
                         samplesDict['DMSpin0'].items())

samplesDict['Data'] = dict(samplesDict['JetHT'].items() +
                           samplesDict['SingleMuon'].items())

samplesDict['All'] = dict(samplesDict['MC'].items() + samplesDict['Data'].items())

for label, isMc in samplesDict['All'].iteritems():
    samplesDict[label] = {label: isMc}
        
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
                      #choices=['All','Hbb','QCD','JetHT','SingleMuon','DMSpin0','TT','DY','W','Diboson','Triboson','SingleTop','VectorDiJet1Jet','VectorDiJet1Gamma','MC','Data'],
                      help="samples to produce")
    parser.add_option('-t','--tag',dest="tag", default = "zprimebits-v12.04", help = "tag, which is the same as folder") 
    parser.add_option('-b','--batch',dest="sub", default = False, help = "use condor or batch system")
    parser.add_option("--njobs-per-file",dest="njobs_per_file",type='int',default=1,help="Split into n jobs per file, will automatically produce submission scripts")
    
    (options,args) = parser.parse_args()

    monitorOption = ''
    if options.monitor is not '':
        monitorOption = '--monitor %s'%options.monitor
    
    jsonPrompt = "$PWD/../data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
    jsonRereco = "$PWD/../data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    jsonPrompt17 = "Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"

    xsec = 1
    analysisDir = options.tag
    executable = "runZprime"
    samples = samplesDict[options.sample]

    if options.sub:
        eosOutDir = '/eos/cms/store/group/phys_exotica/dijet/dazsle'
        execPython = 'baconBatch.py'
        EOS = ''
        optionsDataMc = {
            'mc': "Output.root --passSumEntries 5:Events -a 6:subjob_i -a 7:%i -a 2:mc -a 3:none  -n 8000 -q 2nw4cores --njobs-per-file %d"%(options.njobs_per_file,options.njobs_per_file),
            'data': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonPrompt,options.njobs_per_file),        
            'rereco': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonRereco,options.njobs_per_file),
            }
    else:
        eosOutDir = '/store/user/lpcbacon/dazsle'
        execPython = 'baconCondor.py'
        EOS = 'eos root://cmseos.fnal.gov'
        optionsDataMc = {
            'prompt17': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonPrompt17,options.njobs_per_file),
            'mc': "Output.root --passSumEntries 5:Events -a 6:subjob_i -a 7:%i -a 2:mc -a 3:none  -n 8000 -q 2nw4cores --njobs-per-file %d"%(options.njobs_per_file,options.njobs_per_file),
            }

    exec_me('%s mkdir -p %s/%s'%(EOS,eosOutDir,analysisDir))  
    for label, isMc in samples.iteritems():
        exec_me('%s mkdir -p %s/%s/%s'%(EOS,eosOutDir,analysisDir,label))
        if isMc in ['prompt17']:
            #if 'JetHTRun2017B_12Sep2017_v1' in label or 'JetHTRun2017B_PromptReco_v1' in label or 'SingleMuon' in label:
            #if 'SingleMuon' in label:
            #    exec_me("python %s %s %s -a 4:%f --list 1:../lists/productiontest/%s.txt --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s %s"%(execPython,executable,optionsDataMc[isMc],xsec,label,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
            #else:
            exec_me("python %s %s %s -a 4:%f --list 1:../lists/production13/%s.txt --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s %s"%(execPython,executable,optionsDataMc[isMc],xsec,label,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
        elif isMc in ['data','rereco']:
            exec_me("python %s %s %s -a 4:%f --list 1:../lists/production12a/%s.txt --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s  %s"%(execPython,executable,optionsDataMc[isMc],xsec,label,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
        else:            
            exec_me("python %s %s %s -a 4:%f --list 1:../lists/production13/%s.txt --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s %s"%(execPython,executable,optionsDataMc[isMc],xsec,label,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
            #exec_me("python %s %s %s -a 4:%f --list 1:../lists/production12/%s.txt --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s  %s"%(execPython,executable,optionsDataMc[isMc],xsec,label,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
