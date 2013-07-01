import os
import sys
import cPickle
import pylab
from subprocess import call
import numpy

def ind(modes):
    print "File\tMode\tAlpha\tAUC\tDelta\t#MutEst\t#MutReal"
    with open(os.path.join('..','16S','key_to_Bench.txt')) as f:
        key = eval(f.readline())
    d_key = {x[0]:x[1:] for x in key}
    with open(os.path.join('..','16S','d_clusters_trim.cPickle'),'rb') as f:
        data = cPickle.load(f)
    mutsfiles = {}
    for f_name in os.listdir(os.path.join('..','16S','Mut')):
        pos = f_name.find("id")
        idnum = int(f_name[:pos])
        mutsfiles[idnum] = os.path.join('..','16S','Mut',f_name)
    for f_name in os.listdir(os.path.join('..','16S','Bench_iso')):
        sys.stderr.write("%s\n"%f_name)
        if f_name.endswith('mut'):
            continue
        id = int(f_name.split('.')[0])
        seq = d_key[id][0]
        fold = d_key[id][-1]
        mut = data[seq]['mut'][fold]
        mut_file = mutsfiles[id]
        try:
            os.mkdir('Figs_Bench')
        except OSError:
            pass
        for m in modes:
            cmd = ['python',
                   'createROC.py',
                   '-f', '%s' % os.path.join('..','16S','Bench_iso',f_name),
                   '-o', '%s' % os.path.join('Figs_Bench','%sid_%sfold_%s.pdf' % (id, fold,m)),
                   '-k', '%s' % mut_file,
                   '-m', '%s' % m]
            #print cmd
            call(cmd)
    
def batch(modes):
    print "File\tMode\tAlpha\tAUC\tDelta\t#MutEst\t#MutReal"
    with open(os.path.join('..','16S','key_to_Bench.txt')) as f:
        key = eval(f.readline())
    d_key = {x[0]:x[1:] for x in key}
    with open(os.path.join('..','16S','d_clusters_trim.cPickle'),'rb') as f:
        data = cPickle.load(f)

    l_f_name = os.listdir(os.path.join('..','16S','Bench_iso'))
    l_id = set([x for x in l_f_name if not x.endswith('mut')])
    d_fold = {}
    for id in l_id:
        fold = d_key[int(id.split('.')[0])][-1]
        if fold not in d_fold:
            d_fold[fold] = [id]
        else:
            d_fold[fold].append(id)

    for fold in d_fold:
        for m in modes:
            to_do = d_fold[fold]

            cmd = ['python',
                   'createROC.py',
                   '-f']
            flist = []
            for x in to_do:
                flist.append(os.path.join('..','16S','Bench_iso',x))
            cmd.extend([",".join(flist)])

            cmd.extend(["-k"])
            flist = []
            for x in to_do:
                id = int(x[:x.find(".txt")])
                flist.append(os.path.join('..','16S','Mut',"%sid_%sfold.mut"%(id,fold)))
            cmd.extend([",".join(flist)])
            cmd.extend(['-o', '%s' % os.path.join('Figs_Bench','%sfold-%s.pdf' % (fold,m)),
                   '-m', '%s' % m])
            call(cmd)

def getTimes(path):
    fileHandler = open(path,'r')
    lines = [line[:] for line in fileHandler.readlines()]
    seq = lines[0][:-1]
    struct = lines[1][:-1]
    nbMut = int(lines[2][:-1])
    profiles = eval(lines[3][:-1])
    times = eval(lines[4])
    print len(seq),nbMut,times
    return times

def createTimeFigMultiAl():
    dataPoints = {}
    pylab.figure()
    for yamaska_name in ("TimeBench_yamaska.dat",
                 "TimeBench_yamaska_75.dat",
                         "TimeBench_yamaska_100.dat"):
        for l in open(yamaska_name,"r"):
            data = l[:-1].split()
            print data
            try:
                if len(data)==2:
                    f_name,time = data[0],float(data[1])
                    args = f_name.split("-")
                    length = int(args[1][len("length"):])
                    if (length not in dataPoints):
                        dataPoints[length] = []
                    dataPoints[length].append(time)
            except Exception, e:
                print e
                pass
        x = dataPoints.keys()
        x.sort()
        y = [numpy.mean(dataPoints[v]) for v in x]    
        err = [numpy.std(dataPoints[v]) for v in x]    
        pylab.errorbar(x, y, yerr=err,linewidth=1.5)
    pylab.title("")
    pylab.xlabel('Length (nts)',fontsize=16)
    pylab.ylabel('Time (sec.)',fontsize=16)
    pylab.savefig(os.path.join("Figs_Bench","TimeBenchmark.pdf"), format='pdf',bbox_inches='tight')

def createTimeFig():
    dataPoints = {}
    for l in open("TimeBench_yamaska.dat","r"):
        data = l[:-1].split()
        print data
        try:
            if len(data)==2:
                f_name,time = data[0],float(data[1])
                args = f_name.split("-")
                length = int(args[1][len("length"):])
                if (length not in dataPoints):
                    dataPoints[length] = []
                dataPoints[length].append(time)
        except Exception, e:
            print e
            pass
    x = dataPoints.keys()
    x.sort()
    y = [numpy.mean(dataPoints[v]) for v in x]    
    err = [numpy.std(dataPoints[v]) for v in x]    
    pylab.figure()
    pylab.errorbar(x, y, yerr=err,linewidth=1.5)
    pylab.title("")
    pylab.xlabel('Length (nts)',fontsize=16)
    pylab.ylabel('Time (sec.)',fontsize=16)
    pylab.savefig(os.path.join("Figs_Bench","TimeBenchmark.pdf"), format='pdf',bbox_inches='tight')

if __name__ == '__main__':
    modes = ['a','b']
    batch_flag = False
    time_flag = False
    for i, opt in enumerate(sys.argv):
        if opt == '-m':
            modes = [sys.argv[i + 1]]
        elif opt == '-b':
            batch_flag = True
        elif opt == '-t':
            time_flag = True

                  
    if batch_flag:
        batch(modes)
    elif time_flag:
        createTimeFigMultiAl()
    else:
        ind(modes)

