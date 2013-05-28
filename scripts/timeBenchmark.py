import os,sys,os.path;
from  math import log;
import random
from subprocess import call
import time

def parseSecStr(secstr):
  p = []
  result = []
  for i in range(len(secstr)):
    c = secstr[i]
    if c =="(":
      p.append(i)
    elif c==")":
      j = p.pop()
      result.append((j,i))
  return result
      
BASES = ["A","C","G","U"]
BASE_PAIRS = [("A","U"),
  ("U","A"),
  ("U","G"),
  ("G","C"),
  ("G","U"),
  ("C","G")]

def genMSA(secstr,numseqs):
  result = []
  bps = parseSecStr(secstr)
  for i in range(numseqs):
    seq = [random.choice(BASES) for j in secstr]
    for (x,y) in bps:
      a,b = random.choice(BASE_PAIRS)
      seq[x]=a
      seq[y]=b
    result.append("".join(seq))
  return result

def genRNAPyroInputFiles(basepath,fromlength=50,tolength=500,lengthincr=25,numseqs=25,alignmentlength=50):
  for i in range(fromlength,tolength,25):
    sys.stderr.write("*%s"%(i))
    os.system("grgfreqs --draw RNALoopsOut.ggd %s %s > out.txt 2> dummy"%(i,numseqs))
    j = 1
    for l in open("out.txt"):
      data = l[:-1].replace("H","(").replace("a","(").replace("b",")").replace("m",".").replace("u",".").replace("M",".").split()
      secstr = "".join(data)
      seqs = genMSA(secstr,alignmentlength)
      fname  = "%s-length%s-num%s.dat"%(basepath,i,j)
      print fname
      outfile = open(fname,"w")
      for s in seqs:
        outfile.write(s+"\n")
      outfile.write("".join(data)+"\n")
      outfile.close()    
      j += 1



if __name__=="__main__":
#  if not os.path.exists("TimeBench"):
#    os.mkdir("TimeBench")
#  genRNAPyroInputFiles(os.path.join("TimeBench","MSA"),50,500,25,50,50)
  exec_name = os.path.join('..','src','RNAPyro_iso.py')
  nbmuts = 10
  alpha = .5

  dataPoints = {}
  for f_name in os.listdir(os.path.join("TimeBench")):
    args = f_name.split("-")
    length = int(args[1][len("length"):])
    if (length not in dataPoints):
      dataPoints[length] = []
    dataPoints[length].append(f_name)
  maxSize = max([len(x) for x in dataPoints.values()])

  outfile = open("TimeBenchTst.dat","w")
  lengths = dataPoints.keys()
  lengths.sort()
  for i in range(maxSize):
    for length in lengths:
      if len(dataPoints[length])>i:
        f_name = dataPoints[length][i]
        outfile.write("%s\t"%(f_name))
        cmd = ['python',
                       exec_name,
                       '-f', '%s' % os.path.join("TimeBench",f_name),
                       '-m', '%s' % nbmuts,
                       '-a', '%s' % alpha
                       ]
        tic = time.time()
        call(cmd)
        tac = time.time() - tic
        print "%s\t%s"%(f_name,tac)
        outfile.write("%s\n"%(tac))
        outfile.flush()
