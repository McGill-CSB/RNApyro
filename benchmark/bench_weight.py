import sys
import os
import subprocess
import cPickle

def  fix_bad_substitution(struct,l_rnas,data_2):
  new_data_set = []
  for i,rna in enumerate(l_rnas):
    rnainv = data_2[i]
    new_rna = []
    for j,x in enumerate(rna):
      if struct[j] == '.':
        new_rna.append(rnainv[j])
      else:
        new_rna.append(x)
    new_data_set.append(''.join(new_rna))
  return new_data_set

def benchmark(list_sec_struc,weights,nb,output_cPickle=None):
  """f
  Will output a dictionary with keys (sec_struct,gc_content_target)
  and value a tuple (res,data_dump)
    res = list of sequences with +- 10% of GC target
    data_dump = list updating 'C' weight in profile,
                list of GC contents for sequences in given round

  If a 'output_cPickle' arguement is given, will save the dic in it
  """
  d_all_bench = {}
  tmp_struct_file = '.sec_struct.tmp'
  content_dump = '.content.tmp'
  present_loc = os.path.abspath('.')
  cmd = 'python ../src/RNAPyroProfile.py -d "%s" -a 1 -m 15 -no_profile -s_gc %s %s -gc_sec_struct -gc_data "%s"'
  for struct in list_sec_struct:
    results = []
    with open(tmp_struct_file,'w') as f:
      f.write(struct)
    for w in weights:
      results.append([])
      full_path_struct = os.path.join(present_loc,tmp_struct_file)
      full_path_content = os.path.join(present_loc,content_dump)
      to_do = cmd % (full_path_struct,
                   w,nb,
                   full_path_content)
      res = subprocess.check_output(to_do,shell=True)
      data_dump = [x.strip() for x in open(full_path_content)]
      print struct,w,len(struct),len(data_dump)/2
      d_all_bench[struct,w] = (res,data_dump)
  if output_cPickle:
    with open(output_cPickle, 'wb') as f:
      cPickle.dump(d_all_bench,f,-1)
  os.remove(tmp_struct_file)
  return d_all_bench

def format_rna(struct,rna):
  new_rna = []
  for i,x in enumerate(rna):
    if not struct[i] == '.':
      new_rna.append(x.lower())
    else:
      new_rna.append(x)
  return ''.join(new_rna)

def rnainverse(struct,l_rnas):
  l_rnas = l_rnas.split()
  tmp_rnainverse = 'tmp_file_rnainverse.tmp'
  with open(tmp_rnainverse,'w') as f:
    for r in l_rnas[:-1]:
      r = format_rna(struct,r)
      f.write('%s\n%s\n' % (struct,r)) 
    f.write('%s\n%s' % (struct,l_rnas[-1])) 
  to_do = 'RNAinverse < %s' % tmp_rnainverse
  data = subprocess.Popen(to_do,shell=True,stdout=subprocess.PIPE)
  data = data.communicate()[0]
  data = [x for x in data.split('\n') if x]
  data_2 = []
  for x in data:
    y = x.split()[0]
    if all(z in 'ACGU' for z in y.upper()):
      data_2.append(y.upper())
  os.remove(tmp_rnainverse)
  fix_bad_substitution(struct,l_rnas,data_2)
  return '\n'.join(data_2)

def get_perfect_fit(data_file,output_file):
  tmp_rnafold = 'tmp_file_rnafold.tmp'
  with open(data_file) as f:
    d_results = cPickle.load(f)
  d_good = {}
  to_do = 'RNAfold < %s' % tmp_rnafold 
  for struct,weight in d_results:
    print struct,weight
    inversed_rnas = rnainverse(struct,d_results[struct,weight][0])
    with open(tmp_rnafold, 'w') as f:
      f.write(inversed_rnas)
    data = subprocess.Popen(to_do,shell=True,stdout=subprocess.PIPE)
    data = data.communicate()[0]
    data = [x for x in data.split('\n') if x]
    i = 0
    this_batch_good = []
    while i < len(data):
      if not all(x in 'ACGU' for x in data[i]):
        i += 1
        continue
      if data[i+1].split()[0] == struct:
        this_batch_good.append(data[i])
      i += 2
    d_good[struct,weight] = this_batch_good
  with open(output_file,'wb') as f:
      cPickle.dump(d_good,f,-1)
  os.remove(tmp_rnafold)
  return d_good


if __name__ == '__main__':
  output_file = 'dump_50K_40.cPickle'
  output_good_sequences = 'good_rnastrand_rnainverse_100_highT.cPickle'
  DATA = '50K.tmp'
  list_sec_struct = [x.strip() for x in open(
    os.path.join('..','data',DATA)
    )]
  nb = 50000
  #weights = [0.1,0.3,0.5,0.7,0.9]
  weights=[0.4]
  benchmark(list_sec_struct,weights,nb,output_file)
  #get_perfect_fit(output_file,output_good_sequences)
