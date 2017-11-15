#!/usr/bin/env python

# baconCondor.py #############################################################################
# Python driver for Bacon Analyzer executable in Condor
# Original Author N.Wardle (CERN) 

# ------------------------------------------------------------------------------------

import ROOT as r

import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from numpy import arange
from itertools import product
from BaconAna.Utils.makeFilelist import *

default_args = []
EOS = ''

# Options
parser = OptionParser()
parser = OptionParser(usage="usage: %prog analyzer outputfile [options] \nrun with --help to get list of options")
parser.add_option("-l", "--list", default='', 
                  help="Pick up files from a particular list of files")
parser.add_option("-o", "--outdir", default='bacon',
                  help="output for analyzer. This will always be the output for job scripts.")
parser.add_option("-e", "--eosoutdir", default='', 
                  help="eos output directory for analyzer files.")
parser.add_option("-a", "--args", dest="args", default=[], action="append",
                  help="Pass executable args n:arg OR named arguments name:arg. Multiple args can be passed with <val1,val2...> or lists of integers with [min,max,stepsize]")
parser.add_option("-v", "--verbose", dest="verbose", default=False, action="store_true", 
                  help="Spit out more info")
parser.add_option("", "--passSumEntries", dest="passSumEntries", default="",
                  help="x:treename Get Entries in TTree treename and pass to argument x")
parser.add_option("", "--blacklist", dest="blacklist", default=[], action="append",
                  help="Add blacklist file types (search for this string in files and ignore them")

# Make condor submission scripts options
parser.add_option("-n", "--njobs", dest="njobs", type='int', default=-1,
                  help="Split into n jobs, will automatically produce submission scripts")
parser.add_option("--njobs-per-file", dest="njobs_per_file", type='int', default=1,
                  help="Split into n jobs per file, will automatically produce submission scripts")
parser.add_option("-q", "--queue", default='1nh', 
                  help="submission queue")
parser.add_option("--dryRun", default=False, action="store_true", 
                  help="Do nothing, just create jobs if requested")

# Monitor options (submit,check,resubmit failed)  -- just pass outodir as usual but this time pass --monitor sub --monitor check or --monitor resub
parser.add_option("--monitor", default='', help="Monitor mode (sub/resub/check directory of jobs)")

cwd = os.getcwd()
(options, args) = parser.parse_args()
if len(args) < 2 and not options.monitor: sys.exit('Error -- must specify ANALYZER and OUTPUTNAME')
njobs = options.njobs if options.njobs > 0 else 1
njobs_per_file = options.njobs_per_file
cmssw = "CMSSW_9_2_12"

def write_job(exec_line, out, analyzer, i, n, j, eosout=''):
    cwd = os.getcwd()
    analyzer_short = analyzer.split("/")[-1]

    exec_line = exec_line.replace(analyzer, analyzer_short)
    sub_file = open('%s/sub_%s_job%d_subjob%d.sh' % (out, analyzer_short, i, j), 'w')
    sub_file.write('#!/bin/bash\n')
    sub_file.write('# Job Number %d, running over %d files, subjob %d \n' % (i, n, j))
    sub_file.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    sub_file.write('pwd\n')
    sub_file.write('tar -xf %s.tgz\n'% (cmssw))
    sub_file.write('rm %s.tgz\n'% (cmssw))
    sub_file.write('export SCRAM_ARCH=slc6_amd64_gcc630\n')
    sub_file.write('mkdir -p %s/src\n'% (cmssw))
    sub_file.write('cd CMSSW_9_2_12/src\n')
    sub_file.write('scram b ProjectRename\n')
    sub_file.write('eval `scramv1 runtime -sh`\n')
    sub_file.write('cp ../../data.tgz .\n')
    sub_file.write('mkdir -p ${PWD}/BaconAnalyzer/Analyzer/\n')
    sub_file.write('tar -xvzf data.tgz -C ${PWD}/BaconAnalyzer/Analyzer/\n')
    sub_file.write("export TWD=${PWD}/%s_job%d_subjob%d\n" % (analyzer_short, i, j))
    sub_file.write("mkdir -p $TWD\n")
    sub_file.write("cd $TWD\n")
    sub_file.write("cp ../../../runZprime .\n")
    sub_file.write('%s\n' % exec_line)
    sub_file.write('cd ..\n')
    sub_file.write('rm -rf $TWD\n')
    sub_file.close()
    os.system('chmod +x %s' % os.path.abspath(sub_file.name))


def submit_jobs(lofjobs):
    for sub_file in lofjobs:
        os.system('rm -f %s.stdout' % sub_file)
        os.system('rm -f %s.stderr' % sub_file)
        os.system('rm -f %s.log' % sub_file)
        os.system('rm -f %s.jdl'% sub_file)
        condor_file = open('%s.jdl' % sub_file, 'w')
        condor_file.write('universe = vanilla\n')
        condor_file.write('Executable = %s\n'% sub_file)
        condor_file.write('Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )\n')
        condor_file.write('request_disk = 3000000\n')
        condor_file.write('request_memory = 5000\n')
        condor_file.write('Should_Transfer_Files = YES\n')
        condor_file.write('WhenToTransferOutput = ON_EXIT\n')
        condor_file.write('Transfer_Input_Files = /uscms_data/d3/cmantill/bacon/baconbits/CMSSW_9_2_12/src/BaconAnalyzer/Analyzer/production/cmsset_default.sh, /uscms_data/d3/cmantill/bacon/baconbits/CMSSW_9_2_12.tgz, /uscms_data/d3/cmantill/bacon/baconbits/CMSSW_9_2_12/bin/slc6_amd64_gcc630/runZprime, /uscms_data/d3/cmantill/bacon/baconbits/CMSSW_9_2_12/src/BaconAnalyzer/Analyzer/data.tgz\n')
        condor_file.write('use_x509userproxy = true\n')
        condor_file.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
        condor_file.write('Output = %s.stdout\n' % os.path.abspath(condor_file.name))
        condor_file.write('Error = %s.stderr\n' % os.path.abspath(condor_file.name))
        condor_file.write('Log = %s.log\n' % os.path.abspath(condor_file.name))
        condor_file.write('Queue 1\n')
        condor_file.close()
        os.system('chmod +x %s'% os.path.abspath(condor_file.name))
        os.system('condor_submit %s'%(os.path.abspath(condor_file.name)))

if options.monitor:
    if options.monitor not in ['sub']: sys.exit('Error -- Unknown monitor mode %s' % options.monitor)
    dir = options.outdir

    if options.monitor == 'sub':
        # pick up job scripts in output directory (ends in .sh)
        lofjobs = []
        for root, dirs, files in os.walk(dir):
            for file in fnmatch.filter(files, '*.sh'):
                lofjobs.append('%s/%s' % (os.path.abspath(root), file))
        print 'Submitting %d jobs from directory %s' % (len(lofjobs), dir)
        submit_jobs(lofjobs)


    sys.exit('Finished Monitor -- %s' % options.monitor)

def parse_to_dict(l_list):
    if len(l_list) < 1: return {}
    ret = {}
    nkey = 0
    for item in l_list:
        vargs = item.split(':')  # should put a try here
        ni = vargs[0]
        varg = vargs[1:]
        varg = ":".join(varg)
        if not '-' in ni:
            ni = int(ni)
            nkey += 1
        if not "[" in item and "<" not in item:
            ret[(ni)] = ['', [varg]]
        else:
            if "[" in varg:
                varg = varg.replace("[", "")
                varg = varg.replace("]", "")
                min, max, step = varg.split(",")
                ret[(ni)] = ['', arange(int(min), int(max), int(step))]
            elif "<" in varg:
                varg = varg.replace("<", "")
                varg = varg.replace(">", "")
                largs = varg.split(",")
                ret[(ni)] = ['', largs]

    iskey = 0
    for kr in ret.keys():
        if type(kr) == type(''):
            ll = ret.pop(kr)

            ll[1] = [kr + ' ' + str(l) for l in ll[1]]
            ret[nkey + iskey] = ll
            iskey += 1
    return ret

def getArgsJob(interationsobject, job_id, njobs):
    injobs = []
    #   nf = 0
    ifile = 0
    for ff in iterationsobject:
        if (njobs > 0) and (ifile % njobs != job_id):
            injobs.append((ff, False))
        else:
            injobs.append((ff, True))
        ifile += 1

    return injobs


# -- MAIN
print options.outdir
os.system('mkdir -p %s' % (options.outdir))

mindeces = []
analyzer = args[0]
outfile = args[1]
print options.args
analyzer_args = parse_to_dict(options.args)
print analyzer_args
if options.passSumEntries:
    pos, treenam = options.passSumEntries.split(":")
    numEntries = 0
    if options.list:
        with open(options.list.split(":")[1], 'r') as mylist:
            files = [(myfile.replace('\n', ''), True) for myfile in mylist.readlines()]
        # print files
        for fi in files:
            tf = r.TFile.Open(fi[0])
            try:
                tf.IsOpen()
            except:
                continue
            hf = tf.Get(treenam)  # first try to see if its a histogram
            hist = r.TH1F("tmpA", "tmpA", 1, -1000, 1000)
            try:
                getattr(tf, treenam).Draw("Info.nPUmean>>tmpA", "GenEvtInfo.weight")
                numEntries += hist.Integral()
                print "Adding :", hist.Integral(), " -- ", numEntries
            except:
                try:
                    numEntries += int(getattr(tf, treenam).GetEntries())
                    print "Adding Old :", numEntries
                except:
                    continue
    else:
        numEntries = -1;
    analyzer_args[int(pos)] = ['', ['{:e}'.format(float(numEntries))]]
else:
    if options.list:
        with open(options.list.split(":")[1], 'r') as mylist:
            files = [(myfile.replace('\n', ''), True) for myfile in mylist.readlines()]
            
exec_line = '%s' % analyzer

if options.list:
    filepos, options.list = options.list.split(':')
    analyzer_args[int(filepos)] = ['', "fileinput"]
    
# NEED TO ITERATE OF MAP OF ARGS, FORGET DEFAULT ARGGS I THINK, forec them set!!!!!
sortedkeys = analyzer_args.keys()
if len(sortedkeys): sortedkeys.sort()

for key in sortedkeys:
    #  if arg_i in analyzer_args.keys():
    arg = analyzer_args[key][1]
    if arg == 'fileinput':
        exec_line += ' fileinput '
    elif len(arg) > 1:
        mindeces.append(key)
        exec_line += ' MULTARG_%d ' % key
    else:
        exec_line += ' %s ' % arg[0]

# check that from max to 0 all arguments are accounted for (could always add defaults above) !

for arg_c in range(1, max(analyzer_args.keys())):
    if arg_c not in analyzer_args.keys(): sys.exit("ERROR -- missing argument %d" % arg_c)

print 'running executable -- (default call) \n\t%s' % exec_line

if not options.dryRun and njobs > 0:
    print 'Writing %d Submission Scripts to %s (submit after with --monitor sub)' % (njobs, options.outdir)

if len(files) < njobs:
    njobs = len(files)
for job_i in range(njobs):
    for subjob_i in range(njobs_per_file):
        ################################ WHY does this need to be recreate?
        # This must be the sorted set of keys from the dictionary to build the iterations
        listoflists = [analyzer_args[k][1] for k in sortedkeys]
        # itertools section, make object containing all args to be considered
        # i.e it iterates over all combinations of arguments in the args list
        iterationsobject = product(*listoflists)
        ################################
        if options.list:
            with open(options.list, 'r') as mylist:
                allfiles = [(myfile.replace('\n', ''), True) for myfile in mylist.readlines()]
                #print allfiles
                files = [allfiles[job_i]]
        else:
            files = getArgsJob(iterationsobject, job_i, njobs)  # use itertools to split up any arglists into jobs
        # else: files=[]
        job_exec = ''

        nfiles_i = 0
        for fil_i, fil in enumerate(files):
            if not fil[1]: continue
            if options.list:
                exec_line_i = exec_line.replace('subjob_i','%d'%subjob_i)
                exec_line_i = exec_line_i.replace('fileinput', " " + fil[0] + " ")
            else:
                exec_line_i = exec_line.replace('subjob_i','%d'%subjob_i)
                for i, m in enumerate(fil[0]):  # no defaults so guarantee (make the check) that all of the args are there)
                    exec_line_i = exec_line_i.replace(" MULTARG_%d " % i, " " + str(m) + " ")  # LIST  OVER iterated arguments and produce and replace MULTIARG_i with arguemnt at i in list ?
            if options.eosoutdir:
                if njobs_per_file > 1:                    
                    job_exec += exec_line_i + '; xrdcp -s %s root://cmseos.fnal.gov/%s/%s_job%d_file%d_subjob%d.root; ' % (
                    outfile, options.eosoutdir, outfile, job_i, fil_i, subjob_i)
                else:
                    job_exec += exec_line_i + '; xrdcp -s %s root://cmseos.fnal.gov/%s/%s_job%d_file%d.root; ' % (
                    outfile, options.eosoutdir, outfile, job_i, fil_i)
            else:
                if njobs_per_file > 1:
                    job_exec += exec_line_i + '; mv %s %s/%s_job%d_file%d_subjob%d.root; ' % (
                    outfile, options.outdir, outfile, job_i, fil_i, subjob_i)
                else:
                    job_exec += exec_line_i + '; mv %s %s/%s_job%d_file%d.root; ' % (
                    outfile, options.outdir, outfile, job_i, fil_i)
            nfiles_i += 1
        if options.verbose: print "VERB -- job exec line --> ", job_exec

        if options.dryRun:
            print 'job %d/%d -> ' % (job_i + 1, njobs), job_exec
        elif options.njobs > 0:
            write_job(job_exec, options.outdir, analyzer, job_i, nfiles_i, subjob_i, options.eosoutdir)
        else:
            print "Running: ", job_exec
            os.system(job_exec)
