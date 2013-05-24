import os
import sys
import cPickle
#import art_main
from shutil import rmtree
from tempfile import NamedTemporaryFile as NTF, mkdtemp
from subprocess import check_output
from multiprocessing import Pool
#from parse_stockholm import get_data, l_basepairs
from pprint import pprint
from random import choice, shuffle
from numpy import mean

PROCS = 18
DATA = 'RF00177_seed.stockholm.txt'
NAME = 'RF00177'
IUPACBASES = {
  'A':['A'],
  'C':['C'],
  'G':['G'],
  'U':['U'],
  'N':['A','C','G','U'],
  'R':['A','G'],
  'Y':['C','U'],
  'S':['G','C'],
  'W':['A','U'],
  'K':['G','U'],
  'M':['A','C'],
  'B':['C','G','U'],
  'D':['A','G','U'],
  'H':['A','C','U'],
  'V':['A','C','G']
  }

d_s_trim = {}

def ss_from_l_bp(l_bp,s):
    l = set(x[0] for x in l_bp)
    r = set(x[1] for x in l_bp)
    ss = []
    for i in range(len(s)):
        if s[i] == '.':
            continue
        if i in l:
            ss.append('(')
        elif i in r:
            ss.append(')')
        else:
            ss.append('.')
    return ''.join(ss)

def remove_seq_gap_ss(l_s,ss):
    l_ss = l_basepairs(ss)
    l_s_1 = []
    for s in l_s:
        for x,y in l_ss:
            if s[x] == '.' or s[y] == '.':
                break
        else:
            l_s_1.append(s)

    return l_s_1

def cluster_by_seq_same_gaps(l_s):
    l_c = []
    s_done = set()
    for i,s in enumerate(l_s):
        if s in s_done:
            continue
        l_c.append([s])
        s_done.add(s)
        for s2 in l_s[i+1:]:
            if (all(s2[j] == '.' for j,nt in enumerate(s) 
                     if nt == '.') and 
                    all(s[j] == '.' for j,nt in enumerate(s2)
                     if nt == '.')):
                l_c[-1].append(s2)
                s_done.add(s2)
    return l_c

def rand_iupac(x):
    return choice(IUPACBASES[x])

def cluster_remove_ref_gaps(l_s,ss):
    l_ss = l_basepairs(ss)
    d_s = {}
    for s in l_s:
        to_trim = set(i for i,x in enumerate(s) if x == '.')
        s_trim = ''.join(rand_iupac(x) for i,x in enumerate(s) if i not in to_trim)
        d_s[s_trim] = {}
        l_ss_trim = [x for x in l_ss if not (x[0] in to_trim or
                                             x[1] in to_trim)]
        d_s[s_trim]['ss'] = ss_from_l_bp(l_ss_trim,s)
        d_s[s_trim]['seq'] = []
        for s2 in l_s:
            if s2 == s:
                continue
            d_s[s_trim]['seq'].append(
                ''.join(rand_iupac(x) if x != '.' else '.' for i,x in enumerate(s2) if i not in to_trim)
                )
    return d_s

def get_RNAPyro_out(seq,ref_seq,struct,m,alpha):
    tmp_dir = mkdtemp(dir='.')
    tmp_in = NTF(dir=tmp_dir,delete=False)
    tmp_in.write('%s\n' % seq)
    tmp_in.write('\n'.join(ref_seq))
    tmp_in.write('\n%s' % struct)
    tmp_in.close()

    
    cmd = ['python','/scratch/Vreinh/Applications/RNApyro/src/RNAPyro.py',
           '-f', tmp_in.name, '-m',str(m),'-a',str(alpha), '-p', '15']
    print ' '.join(cmd)

    try:
        out = check_output(cmd)
    except:
        rmtree(tmp_dir)
        print "Error calling RNAPyro", cmd
        sys.exit(1)

    rmtree(tmp_dir)
    out = out.split('\n')
    probs = []
    for x in out:
        if x:
            x = x.split()
            print x
            if len(x) == 5:
                probs.append(x[1:])
            elif len(x) == 1:
                time = x[0]
    return probs, time

def do_benchmarks(d_s_trim):
    for seq in d_s_trim:
        seq = seq.strip()
        mut = art_main.main(seq.replace('U','T'),75,10)
        print seq
        print mut
        continue
        ref_seq = d_s_trim[seq]['seq'][:5]
        struct = d_s_trim[seq]['ss']
        m = int(len(seq)*0.1)
        alpha = 0.5
        prob, time = get_RNAPyro_out(seq,ref_seq,struct,m,alpha)

def l_basepairs(ss):
    l_bp = []
    l = []
    for i, x in enumerate(ss):
        if x == '(':
            l.append(i)
        elif x == ')':
            l_bp.append((l.pop(), i))
    return l_bp

def correct_ss(ss, seq):

    l_bp = l_basepairs(ss)
    new_l_bp = []
    for l,r in l_bp:
        if (seq[l],seq[r]) in (('A', 'U'), ('U', 'A'),
                             ('G', 'C'), ('C', 'G'),
                             ('G', 'U'), ('U', 'G')):
            new_l_bp.append((l,r))

    return ss_from_l_bp(new_l_bp, seq)


def multiproc_benchmark_slave(work):

    id, seq, cover = work

    if cover == 15: m = int(0.0120471730053*len(seq)) # 2*error * len(seq) [first precomputed]
    if cover == 10: m = int(0.0181836135877*len(seq))
    if cover == 3: m = int(0.137103600045*len(seq))
    if cover == 5: m = int(0.0479629838618*len(seq))
    if cover == 7: m = int(0.0265207087236*len(seq))

    #data_s_trim[(id.pop(), seq, ref_seq, struct, m, alpha)]]

    ref_seq = d_s_trim[seq]['seq']
    mut_seq = d_s_trim[seq]['mut'][cover]
    ss = d_s_trim[seq]['ss']

    #ss == correct_ss(ss, seq)
    
    l_alpha = [0., 0.5, 0.8, 1.]
    shuffle(l_alpha)
    d_prob = {x:[] for x in l_alpha}
    d_time = {x:[] for x in l_alpha}
    for alpha in l_alpha:
      prob, time = get_RNAPyro_out(mut_seq,ref_seq,ss,m,alpha)
      d_prob[alpha] = prob
      d_time[alpha] = time

    with open(os.path.join('Bench', '%s.txt' % id), 'w') as f:
      out = [seq, ss, m, d_prob, d_time]
      out = [str(x) for x in out]
      f.write('\n'.join(out))
    
def multiproc_benchmark():
    l_s = d_s_trim.keys() #93
    l_cover = [3, 5, 10, 15]
    l_id = range(len(l_s) * len(l_cover))
    to_do = [(l_id.pop(), s, cover) for s in l_s
                    for cover in l_cover]
    with open('key_to_Bench.txt', 'w') as f:
      f.write(str(to_do))
    shuffle(to_do)

    pool = Pool(PROCS)
    pool.map(multiproc_benchmark_slave,to_do)

def get_mutant(s,cov):
    try:
        mut = art_main.main(s.replace('U','T'),75,cov)
    except:
        return get_mutant(s,cov)
    return mut

def get_dic_mut(s):
    d = {}
    for cov in (3,5,7,10,15):
        d[cov] = get_mutant(s,cov)
    return d

def main():
    global d_s_trim
    """
    d = get_data(DATA)
    d = d[NAME]
    ss = d['ss']
    l_s = [d['sequences'][x]['seq'] for x in d['sequences']] 
    d_s_trim = cluster_remove_ref_gaps(l_s,ss)
    d_tmp = {}
    for s in d_s_trim:
        d_tmp[s] = get_dic_mut(s)

    good = 0
    bad = 0
    for s in d_tmp:
        for m in d_tmp[s]:
            if not d_tmp[s][m] :
                bad += 1 
            else:
                good += 1
    print good, bad
    for s in d_tmp:
        d_s_trim[s]['mut'] = d_tmp[s]
    with open('d_clusters_trim.cPickle', 'w') as f:
        cPickle.dump(d_s_trim,f,-1)
    """
    with open('d_clusters_trim.cPickle') as f:
        d_s_trim = cPickle.load(f)

    multiproc_benchmark()


if __name__ == '__main__':
    main()
