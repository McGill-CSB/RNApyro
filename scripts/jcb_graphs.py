import os
import sys
import cPickle
from subprocess import call

def main(m='a'):
    with open(os.path.join('..','16S','key_to_Bench.txt')) as f:
        key = eval(f.readline())
    d_key = {x[0]:x[1:] for x in key}
    with open(os.path.join('..','16S','d_clusters_trim.cPickle')) as f:
        data = cPickle.load(f)


    for f_name in os.listdir(os.path.join('..','16S','Bench_iso')):
        id = int(f_name.split('.')[0])
        seq = d_key[id][0]
        fold = d_key[id][-1]
        mut = data[seq]['mut'][fold]
        with open(os.path.join('Figs_Bench','%sid_%sfold.mut' % (id, fold)), 'w') as f:
            f.write(mut)
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

if __name__ == '__main__':
    m = 'a'
    for i, opt in enumerate(sys.argv):
        if opt == '-m':
            m = sys.argv[i + 1]
    main(m)
