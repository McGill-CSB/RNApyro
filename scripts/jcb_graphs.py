import os
import sys
import cPickle
from subprocess import call

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


if __name__ == '__main__':
    modes = ['a','b']
    batch_flag = False
    for i, opt in enumerate(sys.argv):
        if opt == '-m':
            modes = [sys.argv[i + 1]]
        elif opt == '-b':
            batch_flag = True

                  
    if batch_flag:
        batch(modes)
    else:
        ind(modes)
