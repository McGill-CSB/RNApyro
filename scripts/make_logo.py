import os,sys
import random
import weblogolib as wl



def logo_from_frequency(data, output_file):
  #must have integer data only
  array_data = wl.array(data)
  data = wl.LogoData.from_counts('ACGU', array_data)
  options = wl.LogoOptions()
  options.color_scheme = wl.colorscheme.nucleotide
  options.unit_name = 'probability'
  options.fineprint = ''
  options.creator_text = ''

  format = wl.LogoFormat(data, options)
  with open(output_file, 'w') as f:
    wl.png_formatter(data, format, f)

def data_table_to_integers(data):
  #must first transforme into integers
  for i, x in enumerate(data):
    data[i] = [int(1000*y) for y in x]
    to_distribute = 1000 - sum(data[i])
    for j in range(to_distribute):
      insert = random.choice(range(4))
      data[i][insert] += 1
  return data

def read_file(file_path):
  data = [x.strip().split()[1:] for x in open(file_path)]
  for i,x in enumerate(data):
    data[i] = [float(y) for y in x]
  return data

def help():
  print """The arguments are <input_file> <output_file>, where <input_file> comes
  from RNAPyroEx"""

if __name__ == '__main__':
  options = sys.argv
  if len(options) < 3:
    help()
    sys.exit(1)
  file = options[1]
  output = options[2]
  data = data_table_to_integers(read_file(file))
  logo_from_frequency(data, output)
