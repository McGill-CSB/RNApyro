import os
import sys
import shlex
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

def clean_align(raw_align,seq):
    l = len(seq)
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
                   'probs':[{x:[d_nt_probs[l_names_sorted[0]][i]] for i,x in enumerate(seq)}]
                  }

    for name in l_names_sorted[1:]:
        seq = d_clean_align[name]['seq']
        start = d_clean_align[name]['start']

        for i,x in enumerate(seq):
            i = i + start
            if i < len(d_seq_probs['probs']):
                if x in d_seq_probs['probs'][i]:
                    d_seq_probs['probs'][i][x].append(d_nt_probs[name][i-start])
                else:
                    d_seq_probs['probs'][i][x] = [d_nt_probs[name][i-start]]
            else:
                d_seq_probs['probs'].append({x:[d_nt_probs[name][i-start]]})


    return d_seq_probs

def get_nt_probs(raw_reads):
    d_nt_probs = {}

    for i,x in enumerate(raw_reads):
        if x.startswith('@'):
            d_nt_probs[x[1:]] = [float((ord(z))-33)/100  for z in raw_reads[i+3]]

    return d_nt_probs

def main(seq,l_reads,fold):
    raw_reads, raw_align = get_raw_data(seq,l_reads,fold) 
    d_clean_align = clean_align(raw_align,seq)
    d_nt_probs = get_nt_probs(raw_reads)
    d_merge_probs = merge_reads(d_clean_align,d_nt_probs)


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
            fold = opts[i+1]

    main(seq,l_reads,fold)
    
