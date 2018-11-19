#!/usr/bin/env python
import sys, commands, os, fnmatch
from optparse import OptionParser

samplesDict = {}
samplesDict['JetHTrereco'] = {
    'JetHTRun2017B_17Nov2017_v1': 'rereco17',
    'JetHTRun2017C_17Nov2017_v1': 'rereco17',
    #'JetHTRun2017D_17Nov2017_v1': 'rereco17',
    'JetHTRun2017E_17Nov2017_v1': 'rereco17',
    #'JetHTRun2017F_17Nov2017_v1': 'rereco17',
}
samplesDict['JetHTrereco_13'] = {
    'JetHTRun2017B_17Nov2017_v1_noPF': 'rereco17',
    'JetHTRun2017C_17Nov2017_v1_noPF': 'rereco17',
    'JetHTRun2017D_17Nov2017_v1_noPF': 'rereco17',
    'JetHTRun2017E_17Nov2017_v1_noPF': 'rereco17',
    'JetHTRun2017F_17Nov2017_v1_noPF': 'rereco17',
}
samplesDict['JetHTprompt'] = {
    'JetHTRun2017B_PromptReco_v1_noPF' : 'prompt17',
    'JetHTRun2017B_PromptReco_v2_noPF' : 'prompt17',
    'JetHTRun2017C_PromptReco_v1_noPF' : 'prompt17',
    'JetHTRun2017C_PromptReco_v2_noPF' : 'prompt17',
    'JetHTRun2017C_PromptReco_v3_noPF' : 'prompt17',
    'JetHTRun2017D_PromptReco_v1_noPF' : 'prompt17',
    'JetHTRun2017E_PromptReco_v1_noPF' : 'prompt17',
    'JetHTRun2017F_PromptReco_v1_noPF' : 'prompt17',
}
samplesDict['SingleMuonrereco'] = { 
    'SingleMuonRun2017B_17Nov2017_v1': 'rereco17',
    'SingleMuonRun2017C_17Nov2017_v1': 'rereco17',
    'SingleMuonRun2017D_17Nov2017_v1': 'rereco17',
    'SingleMuonRun2017E_17Nov2017_v1': 'rereco17',
    'SingleMuonRun2017F_17Nov2017_v1': 'rereco17',
}
samplesDict['Hbb'] = {
    #'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_herwigpp': 'mc',
    #'ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'mc',
    'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8': 'mc',
    'GluGluHToBB_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8': 'mc',
    #'GluGluHToBB_M125_13TeV_powheg_pythia8': 'mc',
    'GluGluHToCC_M125_13TeV_powheg_pythia8': 'mc',
    #'ttHTobb_M125_TuneCP5_13TeV_powheg_pythia8': 'mc',
    #'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix': 'mc',
    #'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'mc',
    #'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'mc',
    #'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'mc',
    #'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'mc',
    }
samplesDict['QCD_8X'] = {
    # 'QCD_HT1000to1500_13TeV_8X': 'mc',
    # 'QCD_HT1000to1500_13TeV_ext_8X': 'mc',
    # 'QCD_HT100to200_13TeV_8X': 'mc',
    # 'QCD_HT1500to2000_13TeV_8X': 'mc',
    # 'QCD_HT1500to2000_13TeV_ext_8X': 'mc',
    # 'QCD_HT2000toInf_13TeV_8X': 'mc',
    # 'QCD_HT2000toInf_13TeV_ext_8X': 'mc',
    # 'QCD_HT200to300_13TeV_8X': 'mc',
    # 'QCD_HT200to300_13TeV_ext_8X': 'mc',
    # 'QCD_HT300to500_13TeV_8X': 'mc',
    # 'QCD_HT300to500_13TeV_ext_8X': 'mc',
    # 'QCD_HT500to700_13TeV_8X': 'mc',
    # 'QCD_HT500to700_13TeV_ext_8X': 'mc',
    'QCD_HT700to1000_13TeV_8X': 'mc',
    # 'QCD_HT700to1000_13TeV_ext_8X': 'mc',
    }
# 5
samplesDict['QCD'] = {       
    'QCD_HT1000to1500_TuneCP5_13TeV_madgraph_pythia8': 'mc',
    'QCD_HT100to200_TuneCP5_13TeV_madgraph_pythia8': 'mc',
    'QCD_HT1500to2000_TuneCP5_13TeV_madgraph_pythia8': 'mc',
    'QCD_HT2000toInf_TuneCP5_13TeV_madgraph_pythia8': 'mc',
    'QCD_HT200to300_TuneCP5_13TeV_madgraph_pythia8': 'mc',
    'QCD_HT300to500_TuneCP5_13TeV_madgraph_pythia8': 'mc',
    'QCD_HT500to700_TuneCP5_13TeV_madgraph_pythia8': 'mc',
    'QCD_HT700to1000_TuneCP5_13TeV_madgraph_pythia8': 'mc',
    }
samplesDict['SingleTop'] = {
    #'ST_s_channel_4f_leptonDecays_TuneCP5_13TeV_amcatnlo_pythia8_noPF': 'mc',
    'ST_t_channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8': 'mc',
    'ST_t_channel_top_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8': 'mc',
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8': 'mc',
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8': 'mc',
    }
samplesDict['Diboson'] = {
    'WW_TuneCP5_13TeV_pythia8': 'mc',
    'WZ_TuneCP5_13TeV_pythia8': 'mc',
    'ZZ_TuneCP5_13TeV_pythia8': 'mc',
    }
samplesDict['WLNu'] = {
    'WJetsToLNu_HT_100To200_TuneCP5_13TeV': 'mc',
    'WJetsToLNu_HT_200To400_TuneCP5_13TeV': 'mc',
    'WJetsToLNu_HT_400To600_TuneCP5_13TeV': 'mc',
    'WJetsToLNu_HT_600To800_TuneCP5_13TeV': 'mc',
    'WJetsToLNu_HT_800To1200_TuneCP5_13TeV': 'mc',
    'WJetsToLNu_HT_1200To2500_TuneCP5_13TeV': 'mc',
    #'W1JetsToLNu_TuneCP5_13TeV': 'mc',
    'W2JetsToLNu_TuneCP5_13TeV': 'mc',
    'W3JetsToLNu_TuneCP5_13TeV': 'mc',
    'W4JetsToLNu_TuneCP5_13TeV': 'mc',
    'WJetsToLNu_TuneCP5_13TeV': 'mc'
    }
# 4
samplesDict['VectorDiJet1Jet'] = {
    # 'VectorDiJet1Jet_madgraph_Mphi50Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi75Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi100Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi115Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi125Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi150Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi175Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi200Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi225Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi250Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi275Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi300Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi350Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi400Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi450Mchi1_13TeV': 'mc',
    # 'VectorDiJet1Jet_madgraph_Mphi500Mchi1_13TeV': 'mc',

    'VectorDiJet1Jet_madgraph_Mphi50Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi75Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi100Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi115Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi125Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi150Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi175Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi200Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi225Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi250Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi275Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi300Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi350Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi400Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi450Mchi3000_13TeV': 'mc',
    'VectorDiJet1Jet_madgraph_Mphi500Mchi3000_13TeV': 'mc',
}
samplesDict['VectorDiJet1Jet2016'] = {
    #'VectorDiJet1Jet_100_madgraph_2016': 'mc',
    #'VectorDiJet1Jet_125_madgraph_2017_noPF': 'mc',
    #'VectorDiJet1Jet_200_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_75_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_50_madgraph_2017_noPF' : 'mc',
    #'VectorDiJet1Jet_100_madgraph_2017_noPF' : 'mc',
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
samplesDict['ZprimeWW'] = {
    'ZprimeToWWToWlepWhad_narrow_M_1000_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_2000_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_4500_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_1200_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_2500_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_600_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_1400_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_3000_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_800_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_1600_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_3500_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_1800_TuneCP5_13TeV_madgraph': 'mc',
    'ZprimeToWWToWlepWhad_narrow_M_4000_TuneCP5_13TeV_madgraph': 'mc',
}

samplesDict['W'] = {
    #'WJetsToQQ_HT400to600_TuneCP5_13TeV_noPF': 'mc',
    #'WJetsToQQ_HT600to800_TuneCP5_13TeV_noPF': 'mc',
    #'WJetsToQQ_HT_800toInf_TuneCP5_13TeV_noPF': 'mc',
    'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV': 'mc',
    'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV': 'mc',
    'WJetsToQQ_HT_800toInf_qc19_3j_TuneCP5_13TeV': 'mc',
    }
samplesDict['DY'] = {
    #'ZJetsToQQ_HT400to600_TuneCP5_13TeV_noPF': 'mc',
    #'ZJetsToQQ_HT600to800_3j_TuneCP5_13TeV_noPF': 'mc',
    #'ZJetsToQQ_HT_800toInf_TuneCP5_13TeV_noPF': 'mc',
    'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV': 'mc',
    'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV': 'mc',
    'ZJetsToQQ_HT_800toInf_qc19_4j_TuneCP5_13TeV': 'mc',
    }
# 5
samplesDict['TT'] = {
    #'TTToHadronic_TuneCP5_13TeV_powheg_pythia8': 'mc',
    #'TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8': 'mc',
    #'TTToSemiLeptonic_WspTgt150_TuneCUETP8M2T4_13TeV_powheg_pythia8': 'mc',
    #'TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8': 'mc',
    'TTJets_TuneCP5_13TeV_amcatnloFXFX_pythia8': 'mc',
    }
samplesDict['TT_8X'] = {
    'TT_TuneCUETP8M2T4_13TeV_powheg_pythia8_8X': 'mc',
}
samplesDict['BulkGrav_8X'] = {
    'BulkGravTohhTohVVhbb_narrow_M_1000_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_1200_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_1400_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_1600_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_1800_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_2000_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_2500_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_3000_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_3500_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_4000_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_4500_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_600_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_650_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_700_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_800_13TeV_madgraph_8X': 'mc',
    'BulkGravTohhTohVVhbb_narrow_M_900_13TeV_madgraph_8X': 'mc',
}
samplesDict['BulkGravWW'] = {
    'BulkGravToWW_narrow_M_500_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_800_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_1000_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_1200_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_1400_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_1600_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_1800_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_2000_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_3000_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_3500_13TeV_madgraph':'mc',
    'BulkGravToWW_narrow_M_4500_13TeV_madgraph':'mc',
}

samplesDict['MC'] = dict(samplesDict['Hbb'].items() +
                         samplesDict['QCD'].items() +
                         samplesDict['SingleTop'].items() +
                         samplesDict['W'].items() +
                         samplesDict['DY'].items() +
                         samplesDict['TT'].items() +
                         samplesDict['Diboson'].items() +
                         samplesDict['VectorDiJet1Jet'].items())

samplesDict['Data'] = dict(samplesDict['JetHTrereco'].items() +
                           samplesDict['SingleMuonrereco'].items())

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
    parser.add_option("--monitor",default='',
                      help="Monitor mode (sub/resub/check directory of jobs)")
    parser.add_option('-s','--sample',dest="sample", default="All",
                      #choices=['All','Hbb','QCD','JetHT','SingleMuon','DMSpin0','TT','DY','W','Diboson','Triboson','SingleTop','VectorDiJet1Jet','VectorDiJet1Gamma','MC','Data'],
                      help="samples to produce")
    parser.add_option('-e','--executable',dest="executable", default="runZprime", 
                      help = "executable name")
    parser.add_option('-t','--tag',dest="tag", default = "zprimebits-v14.01", 
                      help = "tag, which is the same as folder") 
    parser.add_option('-b','--batch',dest="sub", default = False, 
                      help = "use condor or batch system")
    parser.add_option("--njobs-per-file",dest="njobs_per_file",type='int',default=1,
                      help="Split into n jobs per file, will automatically produce submission scripts")
    parser.add_option("--nfiles-per-job", dest="nfiles_per_job", type='int', default=1,
                      help="Split into n files per job, will automatically produce submission scripts")    
    (options,args) = parser.parse_args()

    monitorOption = ''
    if options.monitor is not '':
        monitorOption = '--monitor %s'%options.monitor
    
    jsonPrompt = "$PWD/../data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
    jsonRereco = "$PWD/../data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    jsonPrompt17 = "Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
    jsonRereco17 = "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"

    analysisDir = options.tag
    if 'runHWW' in options.executable:
        analysisDir += '-HWW'
    executable = options.executable
    samples = samplesDict[options.sample]

    if options.sub:
        eosOutDir = '/eos/cms/store/group/phys_exotica/dijet/dazsle'
        execPython = 'baconBatch.py'
        EOS = ''
        optionsDataMc = {
            'mc': "Output.root -a 6:subjob_i -a 7:%i -a 2:mc -a 3:none  -n 8000 -q 2nw4cores --njobs-per-file %d"%(options.njobs_per_file,options.njobs_per_file),
            'data': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonPrompt,options.njobs_per_file),        
            'rereco': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonRereco,options.njobs_per_file),
            }
    else:
        eosOutDir = '/store/user/lpcbacon/dazsle'
        execPython = 'baconCondor.py'
        EOS = ''#eos root://cmseos.fnal.gov'
        optionsDataMc = {
            'prompt17': "-a 4:Output.root -a 5:subjob_i -a 6:%i -a 2:data -a 3:%s -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,jsonPrompt17,options.njobs_per_file,options.nfiles_per_job),
            'rereco17': "-a 4:Output.root -a 5:subjob_i -a 6:%i -a 2:data -a 3:%s -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,jsonRereco17,options.njobs_per_file,options.nfiles_per_job),
            'mc': "-a 4:Output.root -a 5:subjob_i -a 6:%i -a 2:mc -a 3:none -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,options.njobs_per_file,options.nfiles_per_job),
            }

    exec_me('%s mkdir -p /eos/uscms/%s/%s'%(EOS,eosOutDir,analysisDir))  
    for label, isMc in samples.iteritems():
        exec_me('%s mkdir -p /eos/uscms/%s/%s/%s'%(EOS,eosOutDir,analysisDir,label))
        if isMc in ['prompt17','rereco17']:
            #listLabel = '../lists/production14/%s.txt'%label
            listLabel = '../lists/production13/%s.txt'%label
        elif isMc in ['data','rereco']:
            listLabel = '../lists/production12a/%s.txt'%label
        else:
            listLabel = '../lists/production14/%s.txt'%label

        exec_me("python %s %s %s --list 1:%s --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s %s"%(execPython,executable,optionsDataMc[isMc],listLabel,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
