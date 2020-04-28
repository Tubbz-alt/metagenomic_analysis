#!/usr/bin/python3

import os
import subprocess
import argparse
import psutil
import time
import copy
import sys
import select
import tty
import re
import shutil

tasks = {}

s_counts = "stability.trim.contigs.good.unique.good.filter.precluster.unique.pick.pick.pick.opti_mcc.0.03.cons.taxonomy"
i_counts = "stability.trim.contigs.good.unique.pick.unique.precluster.opti_mcc.0.03.cons.taxonomy"

s_cmd = lambda x: '\n'.join(["make.file(inputdir=., type=gz, prefix=stability)",
               "make.contigs(file=stability.files, processors=1)",
               "screen.seqs(fasta=current, maxambig=0, maxlength=" + str(x['maxlen']) +")",
               "unique.seqs()",
               "count.seqs(name=current, group=current)",
               "align.seqs(fasta=current, reference=/qfs/projects/minos/reference/silva.nr_v128.align, processors=1)",
               "screen.seqs(fasta=current, count=current, start=" + str(x['start']) + ", end=" + str(x['end']) + ", maxhomop=8)",
               "filter.seqs(fasta=current, vertical=T, trump=.)",
               "pre.cluster(fasta=current, count=current, diffs=2)",
               "unique.seqs(fasta=current, count=current)",
               "chimera.vsearch(fasta=current, count=current, dereplicate=t)",
               "remove.seqs(fasta=current, accnos=current)",
               "classify.seqs(fasta=current, count=current, reference=/qfs/projects/minos/reference/trainset16_022016.pds/trainset16_022016.pds.fasta, taxonomy=/qfs/projects/minos/reference/trainset16_022016.pds/trainset16_022016.pds.tax, cutoff=80, processors=1)",
               "remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown)",
               "remove.groups(count=current, fasta=current, taxonomy=current, groups=Mock)",
               "dist.seqs(fasta=current, cutoff=0.03, processors=1)",
               "cluster(column=current, count=current, cutoff=0.03)",
               "make.shared(list=current, count=current, label=0.03)",
               "classify.otu(list=current, count=current, taxonomy=current, label=0.03)",
               "phylotype(taxonomy=current)"])

i_cmd = lambda x: '\n'.join(["make.file(inputdir=., type=gz, prefix=stability)",
               "make.contigs(file=stability.files, processors=1)",
               "screen.seqs(fasta=current, maxambig=0, maxlength=" + str(x['maxlen']) +")",
               "unique.seqs(fasta=current)",
               "count.seqs(name=current, group=current)",
               "chimera.vsearch(fasta=current, count=current, dereplicate=t)",
               "remove.seqs(fasta=current, accnos=current)",
               "unique.seqs(fasta=current)",                             
               "pre.cluster(fasta=current, count=current, diffs=2)",
               "pairwise.seqs(fasta=current)",
               "classify.seqs(fasta=current, count=current, reference=/qfs/projects/minos/reference/UNITE/UNITE_public_mothur_02.02.2019.fasta, taxonomy=/qfs/projects/minos/reference/UNITE/UNITE_public_mothur_02.02.2019_taxonomy.txt, cutoff=80, processors=1)",
               "cluster(column=current, count=current, cutoff=0.03)",
               "make.shared(list=current, count=current, label=0.03)",
               "classify.otu(list=current, count=current, taxonomy=current, label=0.03)",
               "phylotype(taxonomy=current)"])

s_params = s_cmd({"maxlen": 500, "start": 6388, "end": 25318})
i_params = i_cmd({"maxlen": 500, "start": 1046, "end": 43116})

def PreProcessDir(pdir):
    for i in os.listdir(pdir):
        cfile = os.path.abspath(os.path.join(pdir,i))
        if os.path.isfile(cfile):
            samp = '-'.join(i.split('-')[0:2])
            sampdir = os.path.abspath(os.path.join(pdir,samp))
            if not os.path.isdir(sampdir):
                os.makedirs(sampdir)
            os.rename(cfile, os.path.abspath(os.path.join(sampdir,i.replace('-','_'))))

def ProcessDir(pdir, pparams, clean, tax_counts):
    for i in os.listdir(pdir):
        cdir = os.path.abspath(os.path.join(pdir,i))
        if os.path.isdir(cdir):
            if clean:
                for j in os.listdir(cdir):
                    if "stability" in j.lower() or "mothur" in j.lower():
                        os.remove(os.path.abspath(os.path.join(cdir,j)))
            if not os.path.exists(os.path.join(cdir,tax_counts)):
                with open(os.path.join(cdir, 'mothur.batch'), 'w') as f:
                    f.write(pparams)
                tasks[cdir] = {'fnx' : (lambda x=cdir: subprocess.Popen(["mothur", "mothur.batch"], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, cwd=x)), 'status' : 0}

def CombineFiles(pdir, tax_counts):
    rmpar = re.compile("\(([^\)]+)\)")
    taxdict = {}
    samps = []
    for i in os.listdir(pdir):
        cdir = os.path.abspath(os.path.join(pdir,i))
        if os.path.isdir(cdir):
            cnts = os.path.join(cdir,tax_counts)
            if os.path.exists(cnts):
                samps += [i]
                with open(cnts) as f:
                    for j in f.readlines()[1:]:
                        line = j.split('\t')
                        tax = rmpar.sub('', line[-1]).replace('\n','').replace('\r','')
                        if tax not in taxdict.keys():
                            taxdict[tax] = {}
                        if i not in taxdict[tax].keys():
                            taxdict[tax][i] = 0
                        taxdict[tax][i] += int(line[1])
            else:
                print("WARNING!: No taxonomy counts for sample: " + i)

    hkeys = sorted(taxdict.keys())
    with open(os.path.join(pdir,'taxonomy_counts.txt'), 'w') as f:
        f.write('Sample_ID\t' + '\t'.join(hkeys) + '\n')
        for i in samps:
            f.write(i)
            for j in hkeys:
                if i in taxdict[j]:
                    f.write('\t' + str(taxdict[j][i]))
                else:
                    f.write('\t0')
            f.write('\n')

def RetreiveFiles(pdir, tax_counts):
    results_dir = os.path.join(pdir,"results")
    os.mkdir(results_dir)
    for i in os.listdir(pdir):
        cdir = os.path.abspath(os.path.join(pdir,i))
        if os.path.isdir(cdir):
            cnts = os.path.join(cdir,tax_counts)
            if os.path.exists(cnts):
                shutil.copyfile(cnts, os.path.join(results_dir, i + ".taxonomy"))
            else:
                print("WARNING!: No taxonomy counts for sample: " + i)
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Mothur wrapper for microbiome analysis.')
    parser.add_argument('--dir', type = str, required = True)
    parser.add_argument('--loc', type = str, required = True)
    parser.add_argument('--cpu', type = int, required = True)
    parser.add_argument('--cln', default = False, action='store_false')
    parser.add_argument('--com', type = int, required = True)
    parser.add_argument('--pre', default = False, action='store_false')
    args = parser.parse_args()

    if args.pre:
        print("Beginning PreProcess")
        PreProcessDir(args.dir)
        quit()
    
    if args.com == 2:
        print("Combining Files")
        if args.loc.lower() == "16s":
            CombineFiles(args.dir, s_counts)
        elif args.loc.lower() == "its":
            CombineFiles(args.dir, i_counts)
        quit()

    if args.com == 3:
        print("Retreiving Files")
        if args.loc.lower() == "16s":
            RetreiveFiles(args.dir, s_counts)
        elif args.loc.lower() == "its":
            RetreiveFiles(args.dir, i_counts)
        quit()    

    print("Preparing Batch")
    if args.loc.lower() == "16s":
        ProcessDir(args.dir, s_params, args.cln, s_counts)
    elif args.loc.lower() == "its":
        ProcessDir(args.dir, i_params, args.cln, i_counts)
    else:
        print("Invalid location!")

    done = False
    buf = ""
    if sys.stdin.isatty():
        tty.setraw(sys.stdin.fileno())
    while done == False:
        que = [i for i, j in tasks.items() if j['status'] == 0]
        run = [i for i, j in tasks.items() if j['status'] == 1]
        don = [i for i, j in tasks.items() if j['status'] == 2]
        err = [i for i, j in tasks.items() if j['status'] == 3]

        cpu_avail = args.cpu - len(run)

        for i in que:
            if cpu_avail > 0:
                tasks[i]['ret'] = tasks[i]['fnx']()
                tasks[i]['status'] = 1
                cpu_avail -= 1
        for i in run:
            tasks[i]['ret'].poll()
            if not psutil.pid_exists(tasks[i]['ret'].pid):
                tasks[i]['status'] = 2

        if len(don) + len(err) == len(tasks):
            done = True


        if sys.stdin.isatty():
            os.system("clear")
            sys.stdout.write('{0}\n\r'.format("Queued: " + str(len(que)) + " Running: " + str(len(run)) + " Done: " + str(len(don)) + " Error: " + str(len(err))))
            sys.stdout.write(">" + buf)
            sys.stdout.flush()
            if select.select([sys.stdin], [], [], 0)[0]:
                s = sys.stdin.read(1)
                if s == '\x03':
                    for i in run:
                        tasks[i]['ret'].kill()
                    quit()
                elif s == '\b':
                    buf = buf[0:-1]
                elif s == '\n':
                    sys.stdout.write("\r\n")
                    try:
                        exec(buf)
                    except Exception as e:
                        sys.stderr.write(str(e))
                        sys.stderr.flush()
                    sys.stdin.read(1)
                    buf = ""
                else:
                    buf += s
            
        time.sleep(0.05)

    if args.com > 1:
        if args.loc.lower() == "16s":
            CombineFiles(args.dir, s_counts)
        elif args.loc.lower() == "its":
            CombineFiles(args.dir, i_counts)

    print("All Done!")
