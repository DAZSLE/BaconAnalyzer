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
parser.add_option("--dat", dest="dat", default=False, action="store_true",
                  help="output is txt instead of root")
# Make condor submission scripts options
parser.add_option("-n", "--njobs", dest="njobs", type='int', default=-1,
                  help="Split into n jobs, will automatically produce submission scripts")
parser.add_option("--njobs-per-file", dest="njobs_per_file", type='int', default=1,
                  help="Split into n jobs per file, will automatically produce submission scripts")
parser.add_option("--nfiles-per-job", dest="nfiles_per_job", type='int', default=1,
                  help="Split into n files per job, will automatically produce submission scripts")
parser.add_option("--dry-run", dest="dryRun", default=False, action="store_true", 
                  help="Do nothing, just create jobs if requested")

# Monitor options (submit,check,resubmit failed)  -- just pass outodir as usual but this time pass --monitor sub --monitor check or --monitor resub
parser.add_option("--monitor", default='', help="Monitor mode (sub/resub/check directory of jobs)")
parser.add_option("--submit", default='', dest="submit", help="submit")

cwd = os.getcwd()
(options, args) = parser.parse_args()
if len(args) < 1 and not options.monitor: sys.exit('Error -- must specify ANALYZER')
njobs = options.njobs if options.njobs > 0 else 1
njobs_per_file = options.njobs_per_file
nfiles_per_job = options.nfiles_per_job
cmssw = os.getenv('CMSSW_VERSION', 'CMSSW_10_2_6')
cmssw_base = os.getenv('CMSSW_BASE', 'CMSSW_10_2_6')

# submit jobs by looping over job scripts in output dir
if options.monitor:
    if options.monitor not in ['sub']: sys.exit('Error -- Unknown monitor mode %s' % options.monitor)
    odir = options.outdir
    if options.monitor == 'sub':
        lofjobs = []
        for root, dirs, files in os.walk(odir):
            for file in fnmatch.filter(files, '*.jdl'):
                lofjobs.append('%s/%s' % (os.path.abspath(root), file))
        for sub_file in lofjobs:
            cwd = os.getcwd()
            os.chdir(odir)
            print 'Submitting from directory %s' % (odir)
            os.system("condor_submit %s"%sub_file)
            os.chdir(cwd)
    sys.exit('Finished Monitor -- %s' % options.monitor)

# write job
def write_job(exec_line, out, analyzer):
    cwd = os.getcwd()
    analyzer_short = analyzer.split("/")[-1]
    exec_line = exec_line.replace(analyzer, analyzer_short)
    sub_file = open('%s/sub_%s_%s.sh' % (out, analyzer_short,out.split("/")[-1]), 'w')
    sub_file.write('#!/bin/bash\n')
    sub_file.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    sub_file.write('pwd\n')
    sub_file.write('ls\n')
    sub_file.write('export SCRAM_ARCH=slc6_amd64_gcc630\n')
    sub_file.write('scramv1 project CMSSW %s\n'%cmssw)
    sub_file.write('tar -xzf %s.tgz\n'% (cmssw))
    sub_file.write('ls %s/bin\n'%cmssw)
    sub_file.write('rm %s.tgz\n'% (cmssw))
    sub_file.write('cd %s/src\n'%(cmssw))
    sub_file.write('eval `scramv1 runtime -sh`\n')
    sub_file.write('echo "CMSSW: "$CMSSW_BASE \n')
    sub_file.write('cp ../../data.tgz .\n')
    sub_file.write('cp ../../baconSubmit.py .\n')
    sub_file.write('cp ../../*.txt .\n')
    sub_file.write('mkdir -p ${PWD}/BaconAnalyzer/Analyzer/\n')
    sub_file.write('tar -xvzf data.tgz -C ${PWD}/BaconAnalyzer/Analyzer/\n')
    sub_file.write('echo %s\n' % exec_line) 
    sub_file.write('%s\n' % exec_line)
    sub_file.write('cd ..\n')
    sub_file.write('rm -rf $TWD\n')
    sub_file.close()
    os.system('chmod +x %s' % os.path.abspath(sub_file.name))

# write condor submission script
def write_condor(odir, exe='runjob.sh', arguments = [], files = [],nqueue=1):
    job_name = odir+'/'+exe.replace('.sh','.jdl')
    out = 'universe = vanilla\n'
    out += 'Executable = %s\n'%exe
    out += 'Should_Transfer_Files = YES\n'
    out += 'WhenToTransferOutput = ON_EXIT\n'
    out += 'Transfer_Input_Files = %s\n'%(','.join(files))
    out += 'request_memory = 2.5GB\n'
    out += 'Output = %s_$(Cluster)_$(Process).stdout\n'%(job_name)
    out += 'Error  = %s_$(Cluster)_$(Process).stderr\n'%(job_name)
    out += 'Log    = %s_$(Cluster)_$(Process).log\n'   %(job_name)
    out += 'Arguments = %s $(Process)\n'%(' '.join(arguments))
    out += 'Queue %i\n'%nqueue
    with open(job_name, 'w') as f:
        f.write(out)
    os.system('chmod +x %s'% os.path.abspath(job_name))

# parse arguments to dictionary
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

if __name__ == "__main__":
    os.system('mkdir -p %s' % (options.outdir))

    analyzer = args[0]
    analyzer_args = parse_to_dict(options.args)
    
    with open(options.list.split(":")[1], 'r') as mylist:
        files = [(myfile.replace('\n', ''), True) for myfile in mylist.readlines()]
    filepos, options.list = options.list.split(':')
    analyzer_args[int(filepos)] = ['', "fileinput"]
    
    # iterate over args and write exec line
    exec_line = analyzer
    sortedkeys = analyzer_args.keys()
    if len(sortedkeys): sortedkeys.sort()
    for key in sortedkeys:
        arg = analyzer_args[key][1]
        if arg == 'fileinput':
            exec_line += ' fileinput '
        elif len(arg) > 1:
            exec_line += ' MULTARG_%d ' % key
        else:
            exec_line += ' %s ' % arg[0]

    # check that from max to 0 all arguments are accounted for (could always add defaults above) !
    for arg_c in range(1, max(analyzer_args.keys())):
        if arg_c not in analyzer_args.keys(): sys.exit("ERROR -- missing argument %d" % arg_c)

    # check njobs
    if len(files) < njobs:
        njobs = len(files)
    if nfiles_per_job > 1:
        njobs = njobs//nfiles_per_job
        if njobs==0 and len(files) > 0:
            njobs = 1
            options.nfiles_per_job = len(files)

    print 'running executable -- (default call) \n\t%s' % exec_line
    job_exec = 'python baconSubmit.py --eosoutdir %s --list %s --njobs-per-file %s --nfiles-per-job %s %s '%(options.eosoutdir,
                                                                                                             options.list.split('/')[-1],
                                                                                                             options.njobs_per_file,
                                                                                                             options.nfiles_per_job,
                                                                                                             exec_line)
    job_exec+= '--index ${1}'

    if not options.dryRun and njobs > 0:
        print 'Writing 1 Submission Script for %s to %s (submit after with --monitor sub)' % (njobs, options.outdir)
        print 'with arg ',job_exec

    print 'njobs ',njobs
    if options.njobs > 0:
        filesTransfer  = ['CMSSW/src/BaconAnalyzer/Analyzer/production/cmsset_default.sh', 
                          'CMSSW/src/BaconAnalyzer/Analyzer/production/baconSubmit.py',
                          'CMSSW/src/BaconAnalyzer/Analyzer/production/%s'%options.list,
                          'CMSSW.tgz', 
                          'CMSSW/src/BaconAnalyzer/Analyzer/data.tgz']
        #filesTransfer = [x.replace('CMSSW',cmssw_base) for x in filesTransfer]
        # this is a temp fix - 24/06/2019
        # so that all guys submit from my dir
        crisdir = '/uscms_data/d3/cmantill/bacon/baconbits/CMSSW_10_2_6'
        filesTransfer = [x.replace('CMSSW',crisdir) for x in filesTransfer]
        arguments = []
        write_job(job_exec, options.outdir, analyzer)
        print 'odir  ',options.outdir
        write_condor(options.outdir,'sub_%s_%s.sh'%(analyzer.split("/")[-1],options.outdir.split("/")[-1]), arguments, filesTransfer,njobs)
    else:
        print "Running: ", job_exec
        os.system(job_exec)

