#!/usr/bin/env python
import sys, commands, os, fnmatch
from optparse import OptionParser

samplesDict = {}
samplesDict['JetHTrereco_8X'] = {
    'JetHTRun2016B_07Aug17_ver1_v1_ddb8X': 'rereco16',
    'JetHTRun2016B_07Aug17_ver2_v1_ddb8X': 'rereco16',
    'JetHTRun2016C_07Aug17_v1_ddb8X': 'rereco16',
    'JetHTRun2016D_07Aug17_v1_ddb8X': 'rereco16',
    'JetHTRun2016E_07Aug17_v1_ddb8X': 'rereco16',
    'JetHTRun2016F_07Aug17_v1_ddb8X': 'rereco16',
    'JetHTRun2016G_07Aug17_v1_ddb8X': 'rereco16',
    'JetHTRun2016H_07Aug17_v1_ddb8X': 'rereco16',
    }
samplesDict['JetHTrereco_8X_withPF'] = {
    'JetHTRun2016B_07Aug17_ver1_v1_8X_withPF': 'rereco16',
    'JetHTRun2016B_07Aug17_ver2_v1_8X_withPF': 'rereco16',
    'JetHTRun2016C_07Aug17_v1_8X_withPF': 'rereco16',
    'JetHTRun2016D_07Aug17_v1_8X_withPF': 'rereco16',
    'JetHTRun2016E_07Aug17_v1_8X_withPF': 'rereco16',
    'JetHTRun2016F_07Aug17_v1_8X_withPF': 'rereco16',
    'JetHTRun2016G_07Aug17_v1_8X_withPF': 'rereco16',
    'JetHTRun2016H_07Aug17_v1_8X_withPF': 'rereco16',
    }
samplesDict['JetHTrereco_9X'] = {
    #'JetHTRun2017B_17Nov2017_v1': 'rereco17',
    #'JetHTRun2017C_17Nov2017_v1': 'rereco17',
    #'JetHTRun2017D_17Nov2017_v1': 'rereco17',
    #'JetHTRun2017E_17Nov2017_v1': 'rereco17',
    'JetHTRun2017F_17Nov2017_v1': 'rereco17',
    #'JetHTRun2017F_17Nov2017_v1_try2': 'rereco17',
    #'JetHTRun2017F_17Nov2017_v1_try2_missingjson': 'rereco17', 
    }
samplesDict['JetHTrereco_9X_13'] = {
    'JetHTRun2017B_17Nov2017_v1_noPF': 'rereco17',
    'JetHTRun2017C_17Nov2017_v1_noPF': 'rereco17',
    'JetHTRun2017D_17Nov2017_v1_noPF': 'rereco17',
    'JetHTRun2017E_17Nov2017_v1_noPF': 'rereco17',
    'JetHTRun2017F_17Nov2017_v1_noPF': 'rereco17',
    }
samplesDict['SingleMuonrereco_9X'] = { 
    'SingleMuonRun2017B_17Nov2017_v1': 'rereco17',
    'SingleMuonRun2017C_17Nov2017_v1': 'rereco17',
    'SingleMuonRun2017D_17Nov2017_v1': 'rereco17',
    'SingleMuonRun2017E_17Nov2017_v1': 'rereco17',
    'SingleMuonRun2017F_17Nov2017_v1': 'rereco17',
}
samplesDict['JetHTprompt_10X'] = {
    # 'JetHTRun2018A_17Sep2018_pilot_v1': 'prompt18',
    # 'JetHTRun2018A_17Sep2018_v1': 'prompt18',
    # 'JetHTRun2018A_22May2018_v1': 'prompt18',
    # 'JetHTRun2018A_PromptReco_v1': 'prompt18', #
    # 'JetHTRun2018A_PromptReco_v2': 'prompt18',
    # 'JetHTRun2018A_PromptReco_v3': 'prompt18', #
    # 'JetHTRun2018B_17Sep2018_v1': 'prompt18',
    # 'JetHTRun2018B_26Sep2018_HEM_v1': 'prompt18',
    # 'JetHTRun2018B_26Sep2018_HEMmitigation_v1': 'prompt18',
    # 'JetHTRun2018B_26Sep2018_v1': 'prompt18',
    # 'JetHTRun2018B_PromptReco_v1': 'prompt18',
    # 'JetHTRun2018B_PromptReco_v2': 'prompt18',
    # 'JetHTRun2018C_PromptReco_v1': 'prompt18',
    # 'JetHTRun2018C_PromptReco_v2': 'prompt18',
    # 'JetHTRun2018C_PromptReco_v3': 'prompt18',
    # 'JetHTRun2018D_PromptReco_v1': 'prompt18',
    'JetHTRun2018D_PromptReco_v2': 'prompt18', #
    #'JetHTRun2018E_PromptReco_v1': 'prompt18',
}                                                                                                                                                                 
samplesDict['SingleMuonprompt_10X'] = {
    'SingleMuonRun2018A_PromptReco_v1': 'prompt18', #
    'SingleMuonRun2018A_PromptReco_v2': 'prompt18',
    'SingleMuonRun2018A_PromptReco_v3': 'prompt18', #
    'SingleMuonRun2018B_PromptReco_v1': 'prompt18',
    'SingleMuonRun2018B_PromptReco_v2': 'prompt18',
    'SingleMuonRun2018C_PromptReco_v1': 'prompt18',
    'SingleMuonRun2018C_PromptReco_v2': 'prompt18', #
    'SingleMuonRun2018C_PromptReco_v3': 'prompt18', #
    'SingleMuonRun2018D_PromptReco_v2': 'prompt18', #
    'SingleMuonRun2018A_22May2018_v1': 'prompt18', #
    'SingleMuonRun2018A_06Jun2018_v1': 'prompt18',
    'SingleMuonRun2018A_17Sep2018_v2': 'prompt18', #
    'SingleMuonRun2018B_17Sep2018_v1': 'prompt18',
    'SingleMuonRun2018C_17Sep2018_v1': 'prompt18',
    }
samplesDict['Hcc_9X'] = {
    'GluGluHToCC_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8': 'mc',
    'VBFHToCC_M-125_13TeV_powheg_pythia8': 'mc',
    'WminusH_HToCC_WToLNu_M125_13TeV_powheg_pythia8': 'mc',
    'WminusH_HToCC_WToQQ_M125_13TeV_powheg_pythia8': 'mc',
    'WplusH_HToCC_WToLNu_M125_13TeV_powheg_pythia8': 'mc',
    'WplusH_HToCC_WToQQ_M125_13TeV_powheg_pythia8': 'mc',
    'ZH_HToCC_ZToLL_M125_13TeV_powheg_pythia8': 'mc',
    'ZH_HToCC_ZToQQ_M125_13TeV_powheg_pythia8': 'mc',
    'ggZH_HToCC_ZToQQ_M125_13TeV_powheg_pythia8': 'mc',
    }
samplesDict['Hbb_9X'] = {
   'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'mc',
   'ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'mc',
   'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8': 'mc',
   'GluGluHToBB_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8': 'mc',
   'GluGluHToBB_M125_13TeV_powheg_pythia8': 'mc',
   'ttHTobb_M125_TuneCP5_13TeV_powheg_pythia8': 'mc',
   'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix': 'mc',
   'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'mc',
   'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'mc',
   'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'mc',
   'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'mc',
   }
samplesDict['Hbb_10X'] = {
   #'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'mc',
   'ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_10X': 'mc',
   'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8_10X': 'mc',
   'GluGluHToBB_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8_10X': 'mc',
   'GluGluHToBB_M125_13TeV_powheg_pythia8_10X': 'mc',
   'ttHTobb_M125_TuneCP5_13TeV_powheg_pythia8_10X': 'mc',
   'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_10X': 'mc',
   'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_10X': 'mc',
   'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_10X': 'mc',
   'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_10X': 'mc',
   'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_10X': 'mc',
   }
samplesDict['QCD_9X'] = {       
    'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8': 'mc',
    'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8': 'mc',
    'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8': 'mc',
    'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8': 'mc',
    'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8': 'mc',
    }
samplesDict['QCD_10X'] = {
    'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_10X': 'mc',
    'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_10X': 'mc',
    'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_10X': 'mc',
    'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_10X': 'mc',
    'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_10X': 'mc',
    }
samplesDict['SingleTop_9X'] = {
    'ST_t_channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8': 'mc',
    'ST_t_channel_top_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8': 'mc',
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8': 'mc',
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8': 'mc',
    'ST_s_channel_4f_leptonDecays_TuneCP5_13TeV_amcatnlo_pythia8': 'mc',
    }
samplesDict['SingleTop_10X'] = {
    'ST_t_channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV_powheg_madspin_pythia8_10X': 'mc',
    'ST_t_channel_top_4f_InclusiveDecays_TuneCP5_13TeV_powheg_madspin_pythia8_10X': 'mc',
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8_10X': 'mc',
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8_10X': 'mc',
    'ST_s_channel_4f_hadronicDecays_TuneCP5_13TeV_madgraph_pythia8_10X': 'mc',
    #'ST_s_channel_4f_leptonDecays_TuneCP5_13TeV_amcatnlo_pythia8': 'mc',
    }
samplesDict['W_9X'] = {
    #'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV': 'mc',
    #'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV': 'mc',
    'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV': 'mc',
    }
samplesDict['W_10X'] = {
    #'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV_PS_10X': 'mc',
    #'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV_PS_10X': 'mc',
    'WJetsToQQ_HT_800toInf_qc19_3j_TuneCP5_13TeV_PS_10X': 'mc',
}
samplesDict['DY_9X'] = {
    'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV': 'mc',
    'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV': 'mc',
    'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV': 'mc',
    }
samplesDict['DY_10X'] = {
    'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV_PS_10X': 'mc',
    'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV_PS_10X': 'mc',
    'ZJetsToQQ_HT_800toInf_qc19_4j_TuneCP5_13TeV_PS_10X': 'mc',
}
samplesDict['TT_9X'] = {
    'TTToHadronic_TuneCP5_13TeV_powheg_pythia8': 'mc',
    'TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8': 'mc',
    'TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8': 'mc',
    }
samplesDict['TT_10X'] = {
    #'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_10X': 'mc',
    #'TTToHadronic_TuneCP5_13TeV-powheg-pythia8_10X': 'mc',
    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_PS_10X': 'mc'
}
samplesDict['VV_9X'] = {
    'WW_TuneCP5_13TeV_pythia8': 'mc',
    'WZ_TuneCP5_13TeV_pythia8': 'mc',
    'ZZ_TuneCP5_13TeV_pythia8': 'mc',
    }
samplesDict['VV_10X'] = {
    'WW_TuneCP5_13TeV_pythia8_10X': 'mc',
    'WZ_TuneCP5_13TeV_pythia8_10X': 'mc',
    'ZZ_TuneCP5_13TeV_pythia8_10X': 'mc',
    }
samplesDict['WLNu_9X'] = {
    # 'WJetsToLNu_HT_100To200_TuneCP5_13TeV': 'mc',
    # 'WJetsToLNu_HT_200To400_TuneCP5_13TeV': 'mc',
    # 'WJetsToLNu_HT_400To600_TuneCP5_13TeV': 'mc',
    # 'WJetsToLNu_HT_600To800_TuneCP5_13TeV': 'mc',
    #'WJetsToLNu_HT_800To1200_TuneCP5_13TeV': 'mc',
    # 'WJetsToLNu_HT_1200To2500_TuneCP5_13TeV': 'mc',
    #'WJetsToLNu_HT_2500ToInf_TuneCP5_13TeV': 'mc',
    }
samplesDict['WLNu_10X'] = {
    # 'WJetsToLNu_HT_100To200_TuneCP5_13TeV_10X': 'mc',
    # 'WJetsToLNu_HT_200To400_TuneCP5_13TeV_10X': 'mc',
    # 'WJetsToLNu_HT_400To600_TuneCP5_13TeV_10X': 'mc',
    # 'WJetsToLNu_HT_600To800_TuneCP5_13TeV_10X': 'mc',
    'WJetsToLNu_HT_800To1200_TuneCP5_13TeV_10X': 'mc',
    #'WJetsToLNu_HT_1200To2500_TuneCP5_13TeV_10X': 'mc',
    #'WJetsToLNu_HT_2500ToInf_TuneCP5_13TeV_10X': 'mc',
}
samplesDict['DYLL_9X'] = {
    # 'DYJetsToLL_M_50_HT_100to200_TuneCP5_13TeV': 'mc',
    # 'DYJetsToLL_M_50_HT_200to400_TuneCP5_13TeV': 'mc',
    # 'DYJetsToLL_M_50_HT_400to600_TuneCP5_13TeV': 'mc',
    # 'DYJetsToLL_M_50_HT_600to800_TuneCP5_13TeV': 'mc',
    # 'DYJetsToLL_M_50_HT_800to1200_TuneCP5_13TeV': 'mc',
    'DYJetsToLL_M_50_HT_1200to2500_TuneCP5_13TeV': 'mc',
    'DYJetsToLL_M_50_HT_2500toInf_TuneCP5_13TeV': 'mc',
    }
samplesDict['DYLL_10X'] = {
    'DYJetsToLL_M_50_TuneCP5_13TeV_10X': 'mc',
    #'DY3JetsToLL_M_50_TuneCP5_13TeV_10X': 'mc',
    }

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
    parser.add_option('-s','--sample',dest="sample", default="",
                      help="samples to produce")
    parser.add_option('-e','--executable',dest="executable", default="runZprime", 
                      help = "executable name")
    parser.add_option('-t','--tag',dest="tag", default = "zprimebits-v15.01", 
                      help = "tag, which is the same as folder") 
    parser.add_option('--production',dest="production", default = "15",
                      help="bacon production") 
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
    
    jsonPrompt16 = "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
    jsonRereco16 = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    jsonPrompt17 = "Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
    jsonRereco17 = "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
    jsonPrompt18 = "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"

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
            'data': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonPrompt16,options.njobs_per_file),        
            'rereco': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonRereco16,options.njobs_per_file),
            }
    else:
        eosOutDir = '/store/user/lpcbacon/dazsle'
        if 'runHWW' in options.executable:
            eosOutDir = '/store/user/cmantill/baconbits/'
        execPython = 'baconCondor.py'
        EOS = ''#eos root://cmseos.fnal.gov'
        optionsDataMc = {
            'rereco16': "-a 4:Output.root -a 5:subjob_i -a 6:%i -a 2:data -a 3:%s -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,jsonRereco16,options.njobs_per_file,options.nfiles_per_job),
            'prompt17': "-a 4:Output.root -a 5:subjob_i -a 6:%i -a 2:data -a 3:%s -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,jsonPrompt17,options.njobs_per_file,options.nfiles_per_job),
            'rereco17': "-a 4:Output.root -a 5:subjob_i -a 6:%i -a 2:data -a 3:%s -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,jsonRereco17,options.njobs_per_file,options.nfiles_per_job),
            'prompt18': "-a 4:Output.root -a 5:subjob_i -a 6:%i -a 2:data -a 3:%s -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,jsonPrompt18,options.njobs_per_file,options.nfiles_per_job),
            'mc': "-a 4:Output.root -a 5:subjob_i -a 6:%i -a 2:mc -a 3:none -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,options.njobs_per_file,options.nfiles_per_job),
            }

    exec_me('%s mkdir -p /eos/uscms/%s/%s'%(EOS,eosOutDir,analysisDir))  
    for label, isMc in samples.iteritems():
        if '_10X' in options.sample:
            labelOut = label.replace('_10X','')
            if '_PS' in labelOut:
                labelOut = labelOut.replace('_PS','')
        else:
            labelOut = label
        exec_me('%s mkdir -p /eos/uscms/%s/%s/%s'%(EOS,eosOutDir,analysisDir,labelOut))
        listLabel = '../lists/production%s/%s.txt'%(options.production,label)
        exec_me("python %s %s %s --list 1:%s --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s %s"%(execPython,executable,optionsDataMc[isMc],listLabel,analysisDir,label,isMc,eosOutDir,analysisDir,labelOut,monitorOption),options.dryRun)
