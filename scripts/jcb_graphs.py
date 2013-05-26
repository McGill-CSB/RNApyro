import os
import sys
import cPickle
from subprocess import call

def ind(m='a'):
    with open(os.path.join('..','16S','key_to_Bench.txt')) as f:
        key = eval(f.readline())
    d_key = {x[0]:x[1:] for x in key}
    with open(os.path.join('..','16S','d_clusters_trim.cPickle')) as f:
        data = cPickle.load(f)
    

    for f_name in os.listdir(os.path.join('..','16S','Bench_iso')):
        if f_name.endswith('mut'):
            continue
        id = int(f_name.split('.')[0])
        seq = d_key[id][0]
        fold = d_key[id][-1]
        mut = data[seq]['mut'][fold]
        try:
            os.mkdir('Figs_Bench')
        except OSError:
            pass
        cmd = ['python',
               'createROC.py',
               '-f', '%s' % os.path.join('..','16S','Bench_iso',f_name),
               '-o', '%s' % os.path.join('Figs_Bench','%sid_%sfold_%s.pdf' % (id, fold,m)),
               '-m', '%s' % m]
        print cmd
        call(cmd)

def batch(m='a'):
    with open(os.path.join('..','16S','key_to_Bench.txt')) as f:
        key = eval(f.readline())
    d_key = {x[0]:x[1:] for x in key}
    with open(os.path.join('..','16S','d_clusters_trim.cPickle')) as f:
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
        to_do = d_fold[fold]

        cmd = ['python',
               'createROC.py',
               '-f']
        for x in to_do:
            cmd.append(os.path.join('..','16S','Bench_iso',x))

        cmd.extend(['-o', '%s' % os.path.join('Figs_Bench','%sfold.pdf' % (fold)),
               '-m', '%s' % m])
        print cmd
        continue
        call(cmd)


if __name__ == '__main__':
    m = 'a'
    batch_flag = False
    for i, opt in enumerate(sys.argv):
        if opt == '-m':
            m = sys.argv[i + 1]
        elif opt == '-b':
            batch_flag = True

                  
    if batch_flag:
        batch(m)
    else:
        ind(m)
