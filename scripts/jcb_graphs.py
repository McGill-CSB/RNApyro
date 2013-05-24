import os
import sys
from subprocess import call

def main(m='a'):
    with open(os.path.join('..','16S','key_to_Bench.txt')) as f:
        key = eval(f.readline())
    d_key = {x[0]:x[1:] for x in key}

    for f_name in os.listdir(os.path.join('..','16S','Bench')):
        id = int(f_name.split('.')[0])
        fold = d_key[id][-1]
        try:
            os.mkdir('Figs_Bench')
        except OSError:
            pass
        cmd = ['python',
               'createROC.py',
               '-f', '%s' % os.path.join('..','16S','Bench',f_name),
               '-o', '%s' % os.path.join('Figs_Bench','%sid_%sfold.pdf' % (id, fold)),
               '-m', '%s' % m]
        print cmd
        call(cmd)

if __name__ == '__main__':
    m = 'a'
    for i, opt in enumerate(sys.argv):
        if opt == '-m':
            m = sys.argv[i + 1]
    main(m)
