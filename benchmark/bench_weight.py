import sys
import os
import subprocess
import cPickle

def benchmark(list_sec_struc,weights,nb,output_cPickle=None):
  """f
  Will output a dictionary with keys (sec_struct,gc_content_target)
  and value a tuple (res,data_dump)
    res = list of sequences with +- 10% of GC target
    data_dump = list alterning 'C' weight in profile,
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

def get_perfect_fit(data_file,output_file):
  tmp_rnafold = 'tmp_file_rnafold.tmp'
  with open(data_file) as f:
    d_results = cPickle.load(f)
  d_good = {}
  to_do = 'RNAfold < %s' % tmp_rnafold 
  for struct,weight in d_results:
    print struct,weight
    with open(tmp_rnafold, 'w') as f:
      f.write(d_results[struct,weight][0])
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
  output_file = 'dump_rnastrand_dataset_filtered_nodup_10000.cPickle'
  output_good_sequences = 'good_sequences_rnastrand_dataset_10000.cPickle'
  list_sec_struct = [x.strip() for x in open(
    os.path.join('..','data','rnastrand_dataset_filtered_nodup.txt')
    )]
  nb = 10000
  weights = [0.1,0.3,0.5,0.7,0.9]
  benchmark(list_sec_struct,weights,nb,output_file)
  get_perfect_fit(output_file,output_good_sequences)
