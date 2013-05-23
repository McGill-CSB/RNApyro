import os
import sys
import shlex
from random import choice
from subprocess import call
from tempfile import mkdtemp, NamedTemporaryFile as NTF 
from shutil import rmtree
from pprint import pprint
from itertools import combinations


def get_raw_data(seq,l_reads,fold):
    dir_temp = mkdtemp(dir='.')
    cwd = os.getcwd()
    os.chdir(dir_temp)

    temp_file = NTF(dir='.',delete=False,suffix='.fasta')
    temp_file.file.write('>seq\n%s' % seq)
    temp_file.file.close()

    try:
        cmd = 'art_illumina -i %s -l %s -f %s -o out' % (
                temp_file.name,
                l_reads,
                fold)
        cmd  = shlex.split(cmd)
        call(cmd)
    except:
        os.chdir(cwd)
        rmtree(dir_temp)
        print "Error with command:"
        print '\t', cmd
        sys.exit(1)

    raw_reads = [x.strip() for x in open('out.fq')]
    raw_align = [x.strip() for x in open('out.aln')]

    os.chdir(cwd)
    rmtree(dir_temp)

    return raw_reads, raw_align

def reverse_sens(seq):
    dic = {'A':'U','U':'A','C':'G','G':'C'}
    return ''.join(dic[x] for x in reversed(seq))

def clean_align(raw_align,l_seq):
    l = l_seq
    d_clean = {}

    for i,x in enumerate(raw_align):
        if x.startswith('>'):
            _, name, start, orientation = x.split()
            start = int(start)
            d_clean[name] = {}

            if orientation == '+':
                d_clean[name]['ref'] = raw_align[i+1].replace('T','U')
                d_clean[name]['seq'] = raw_align[i+2].replace('T','U')
                d_clean[name]['start'] = start
            else:
                d_clean[name]['ref'] = reverse_sens(raw_align[i+1].replace('T','U'))
                d_clean[name]['seq'] = reverse_sens(raw_align[i+2].replace('T','U'))
                d_clean[name]['start'] = l - start - len(d_clean[name]['seq'])

    return d_clean

def check_align(wt,ref,seq,start):
    wt = wt[start:start+len(ref)]
    if wt != ref:
        return False
    elif '.' in ref+seq:
        return False
    return True

def merge_reads(d_clean_align,d_nt_probs):
    l_names_sorted = sorted(d_clean_align.keys(),
                            key=lambda x: d_clean_align[x]['start'])

    first_pos = d_clean_align[l_names_sorted[0]]['start']
    seq = d_clean_align[l_names_sorted[0]]['seq']
    d_seq_probs = {'start':first_pos,
                   'probs':{i+first_pos:{x:[d_nt_probs[l_names_sorted[0]][i]]} for i,x in enumerate(seq)}
                  }
    
    for name in l_names_sorted[1:]:
        seq = d_clean_align[name]['seq']
        start = d_clean_align[name]['start']

        for i,x in enumerate(seq):
            i = i + start
            if i in d_seq_probs['probs']:
                if x in d_seq_probs['probs'][i]:
                    d_seq_probs['probs'][i][x].append(d_nt_probs[name][i-start])
                else:
                    d_seq_probs['probs'][i][x] = [d_nt_probs[name][i-start]]
            else:
                d_seq_probs['probs'][i] = {x:[d_nt_probs[name][i-start]]}


    return d_seq_probs

def get_nt_probs(raw_reads):
    d_nt_probs = {}

    for i,x in enumerate(raw_reads):
        if x.startswith('@'):
            d_nt_probs[x[1:]] = [float((ord(z))-33)/100  for z in raw_reads[i+3]]

    return d_nt_probs

def list_ss_bp(ss):
    d_lbp = {}
    d_rbp = {}
    left = []
    for i,x in enumerate(ss):
        if x == '(':
            left.append(i)
        elif x == ')':
            l = left.pop()
            r = i
            d_lbp[l] = r
            d_rbp[r] = l
    return d_lbp,d_rbp

def get_max_vote(d_merge_probs,ss):
    seq = []
    for pos in sorted(d_merge_probs['probs'].keys()):
        t = tuple((x,len(d_merge_probs['probs'][pos][x])) 
                  for x in d_merge_probs['probs'][pos])
        nt = sorted(t, key=lambda z:z[1])[-1][0] 
        if pos > len(seq):
            seq.extend(['-']*(pos-len(seq)))
        elif pos < len(seq):
             print "problem get_max_vote order probs"
             sys.exit(1)
        seq.append(nt)
    return ''.join(seq)


def get_most_reads_nt(d):
    listed_keys = sorted(d.keys(), key=lambda x:len(d[x]))
    return listed_keys[-1]


def get_most_reads_seq(d_merge_probs, seq):
    new_seq = []
    for i,x in enumerate(seq):
        if i in d_merge_probs['probs']:
            new_seq.append(get_most_reads_nt(d_merge_probs['probs'][i]))
        else:
            new_seq.append(choice('ACGU'))
    return ''.join(new_seq)

def main(seq,l_reads,fold):
    seq = seq.upper().replace('U','T')
    if any(x not in ('ACGT') for x in seq):
        print help()
        sys.exit(0)
    raw_reads, raw_align = get_raw_data(seq,l_reads,fold) 
    seq = seq.replace('T','U')
    l_seq = len(seq)
    d_clean_align = clean_align(raw_align,l_seq)
    d_nt_probs = get_nt_probs(raw_reads)
    d_merge_probs = merge_reads(d_clean_align,d_nt_probs)
    new_seq = get_most_reads_seq(d_merge_probs, seq)
    print new_seq


def help():
    print """ 
        mandatory arguments:
            -s <string AGC(U|T)> Sequence (DNA or RNA)

            -l <int> length of reads
                re-calibrated:
                    36
                    44
                    50
                    75
            -f <int> average number of fold


        ATTENTION the first few and last nuclotides might
        not be sequenced

        OUTPUT:
            penultimate line: input sequence (as RNA)
            ultimate line   : simulated merge reads (already 
                aligned, majority call, can be incomplete!!!!)
          """

if __name__ == '__main__':
    opts = sys.argv
    
    seq = ''
    l_reads = 0
    fold = 0

    for i,x in enumerate(opts):
        if x == '-l':
            l_reads = int(opts[i+1])
        elif x == '-s':
            seq = opts[i+1].upper().replace('U','T')
        elif x == '-f':
            fold = int(opts[i+1])

    if not seq or not l_reads or not fold:
        print seq, l_reads,fold
        help()
        sys.exit(1)

    main(seq,l_reads,fold)
    
