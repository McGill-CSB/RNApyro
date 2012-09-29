import os,sys
import itertools
import math
from mpmath import *

sys.setrecursionlimit(10000)

BASES = ['A','C','G','U']
BOLTZMANN = 0.0019872041
T = 310.15

#Just generate all possible combination with maxint
STACKING_ENERGY = {k:sys.maxint for k in itertools.product(
                    BASES, repeat=4)}
#Adjust with turner04
#The order of the nucleotides is from 5' -> 3'
STACKING_ENERGY.update({('A', 'A', 'U', 'U'):-0.9,
                        ('A', 'C', 'G', 'U'):-2.2,
                        ('A', 'G', 'C', 'U'):-2.1,
                        ('A', 'G', 'U', 'U'):-0.6,
                        ('A', 'U', 'A', 'U'):-1.1,
                        ('A', 'U', 'G', 'U'):-1.4,
                        ('C', 'A', 'U', 'G'):-2.1,
                        ('C', 'G', 'U', 'G'):-1.4,
                        ('C', 'U', 'A', 'G'):-2.1,
                        ('C', 'U', 'G', 'G'):-2.1,
                        ('G', 'A', 'U', 'C'):-2.4,
                        ('G', 'C', 'G', 'C'):-3.4,
                        ('G', 'G', 'C', 'C'):-3.3,
                        ('G', 'G', 'U', 'C'):-1.5,
                        ('G', 'U', 'A', 'C'):-2.2,
                        ('G', 'U', 'G', 'C'):-2.5,
                        ('G', 'A', 'U', 'U'):-1.3,
                        ('G', 'C', 'G', 'U'):-2.5,
                        ('G', 'G', 'C', 'U'):-2.1,
                        ('G', 'G', 'U', 'U'):-0.5,
                        ('G', 'U', 'A', 'U'):-1.4,
                        ('G', 'U', 'G', 'U'):1.3,
                        ('U', 'A', 'U', 'A'):-1.3,
                        ('U', 'C', 'G', 'A'):-2.4,
                        ('U', 'G', 'C', 'A'):-2.1,
                        ('U', 'G', 'U', 'A'):-1.0,
                        ('U', 'U', 'A', 'A'):-0.9,
                        ('U', 'U', 'G', 'A'):-1.3,
                        ('U', 'A', 'U', 'G'):-1.0,
                        ('U', 'C', 'G', 'G'):-1.5,
                        ('U', 'G', 'C', 'G'):-1.4,
                        ('U', 'G', 'U', 'G'):0.3,
                        ('U', 'U', 'A', 'G'):-0.6,
                        ('U', 'U', 'G', 'G'):-0.5,
                        ('C', 'C', 'G', 'G'):-3.3,
                        ('C', 'G', 'C', 'G'):-2.4})

ISO = {((k1,k2),(k3,k4)):sys.maxint for (k1,k2,k3,k4) in 
                             itertools.product(BASES, repeat=4)}
ISO.update({(('A', 'U'), ('A', 'U')): 0.0,
            (('A', 'U'), ('C', 'G')): 0.34,
            (('A', 'U'), ('G', 'C')): 0.21,
            (('A', 'U'), ('G', 'U')): 2.11,
            (('A', 'U'), ('U', 'A')): 0.31,
            (('A', 'U'), ('U', 'G')): 2.4,
            (('C', 'G'), ('A', 'U')): 0.34,
            (('C', 'G'), ('C', 'G')): 0.0,
            (('C', 'G'), ('G', 'C')): 0.26,
            (('C', 'G'), ('G', 'U')): 2.39,
            (('C', 'G'), ('U', 'A')): 0.21,
            (('C', 'G'), ('U', 'G')): 2.14,
            (('G', 'C'), ('A', 'U')): 0.21,
            (('G', 'C'), ('C', 'G')): 0.26,
            (('G', 'C'), ('G', 'C')): 0.0,
            (('G', 'C'), ('G', 'U')): 2.14,
            (('G', 'C'), ('U', 'A')): 0.34,
            (('G', 'C'), ('U', 'G')): 2.39,
            (('G', 'U'), ('A', 'U')): 2.11,
            (('G', 'U'), ('C', 'G')): 2.39,
            (('G', 'U'), ('G', 'C')): 2.14,
            (('G', 'U'), ('G', 'U')): 0.0,
            (('G', 'U'), ('U', 'A')): 2.4,
            (('G', 'U'), ('U', 'G')): 4.48,
            (('U', 'A'), ('A', 'U')): 0.31,
            (('U', 'A'), ('C', 'G')): 0.21,
            (('U', 'A'), ('G', 'C')): 0.34,
            (('U', 'A'), ('G', 'U')): 2.4,
            (('U', 'A'), ('U', 'A')): 0.0,
            (('U', 'A'), ('U', 'G')): 2.11,
            (('U', 'G'), ('A', 'U')): 2.4,
            (('U', 'G'), ('C', 'G')): 2.14,
            (('U', 'G'), ('G', 'C')): 2.39,
            (('U', 'G'), ('G', 'U')): 4.48,
            (('U', 'G'), ('U', 'A')): 2.11,
            (('U', 'G'), ('U', 'G')): 0.0})

class memoize(object):
    """Generically memoizes a function results."""
    cache = {}
    fun = None
    
    def __init__(self, f):
        self.fun = f
    
    def __call__(self,seq,struct,*args):
        nargs = (args)
        if nargs in self.cache:
            return self.cache[nargs]
        else:
            val = self.fun(seq,struct,*args)
            self.cache[nargs] = val
            return val
    def resetCache(self):
        self.cache = {}

def delta(seq,i,c):
  if i<0 or i>=len(seq):
    return 0
  if c==seq[i]:
    return 0
  else:
    return 1

def energy((a,b),(a2,b2),alpha=1.0):
  #stacking energy of base pair (a,b) around base pair (a2,b2)
  E = STACKING_ENERGY[a,a2,b2,b]
  return  math.exp(-(alpha*E)/(BOLTZMANN*T))

def isostericity(seq,(i,j),(a,b), alpha=1.0):
  iso = ISO[(seq[i],seq[j]),(a,b)]
  return  math.exp(-((1-alpha)*iso)/(BOLTZMANN*T))

@memoize
def forward(seq,struct,(i,j),(a,b),m, alpha=1.0):
  #alpha gives the weight energy vs isostericity
  result = 0.
  if m<0: return mpf(0)
  if i > j :
    if m==0:
      result=1.
    else:
      result=0.
  else:
    k = struct[i]
    if k==-1:
      for a2 in BASES:
        d = delta(seq,i,a2)
        if d <= m:
          result += forward(seq,struct,
                            (i+1,j),
                            (a2,b),
                            m-d)
    elif i < k <= j: #If k > j we return 0
      for a2 in BASES:
        for b2 in BASES:
          d = delta(seq,i,a2)+delta(seq,k,b2)
          #if not stacked outside (or border, then no stack possible)
          if i==0 or j==len(seq)-1 or not (j==k and struct[i-1]==j+1):
            for m2 in range(m-d+1):
              result += forward(seq,struct,
                                  (i+1,k-1),
                                  (a2,b2),
                                  m2)*forward(seq,struct,
                                  (k+1,j),
                                  (b2,b),
                                  m-m2-d)*isostericity(seq,
                                                       (i,k),
                                                       (a2,b2),
                                                       alpha)

          #if stack, we add energy
          else :
            result += forward(seq,struct,
                              (i+1,k-1),
                              (a2,b2),
                              m-d)*energy((a,b),
                                          (a2,b2),
                                          alpha)*isostericity(seq,
                                                              (i,k),
                                                              (a2,b2),
                                                              alpha)
  return mpf(result)

@memoize
def backward(seq,struct,(i,j),(a,b),m,alpha=1.0):
  result = mpf(0.)
  if m<0: return mpf(0)
  if i<0:
    result=forward(seq,struct,
                        (j,len(seq)-1),
                        ('X','X'),
                         m)
  else:
    k = struct[i]
    if k==-1:
      for a2 in BASES:
        d = delta(seq,i,a2)
        result += backward(seq,struct,
                        (i-1,j),
                        (a2,b),
                         m-d)
    #BP to the left
    elif k<i:
      for a2 in BASES:
        for b2 in BASES:
          d = delta(seq,i,b2) + delta(seq,k,a2)
          for m2 in range(m-d+1):
            result += backward(seq,struct,
                                  (k-1,j),
                                  (a2,b),
                                  m-m2-d)*forward(seq,struct,
                                                  (k+1,i-1),
                                                  (a2,b2),
                                                  m2)*isostericity(seq,
                                                                   (k,i),
                                                                   (a2,b2),
                                                                   alpha)
            #print "  L",i,j,m,"B:",m-m2-d,b,"F:(",(k+1,i-1),")",m2,f
    #BP to the right
    elif k>=j:
      for a2 in BASES:
        for b2 in BASES:
          d = delta(seq,i,a2) + delta(seq,k,b2)
          if not (j==k and struct[i+1]==j-1):
            for m2 in range(m-d+1):
              result += backward(seq,struct,
                                 (i-1,k+1),
                                 (a2,b2),
                                 m-m2-d)*forward(seq,struct,
                                                 (j,k-1),
                                                 (b,b2),
                                                 m2)*isostericity(seq,
                                                                  (i,k),
                                                                  (a2,b2),
                                                                  alpha)
          else:
            result +=  backward(seq,struct,
                                 (i-1,k+1),
                                 (a2,b2),
                                 m-d)*energy((a2,b2),
                                             (a,b),
                                             alpha)*isostericity(seq,
                                                                 (i,k),
                                                                 (a2,b2),
                                                                 alpha)
  return mpf(result)

def parseStruct(dbn):
  p = []
  result = [-1 for c in dbn]
  for i in range(len(dbn)):
    c = dbn[i]
    if c=='(':
      p.append(i)
    elif c==')':
      j = p.pop()
      result[j] = i
      result[i] = j
  return result


def testSingleSequence(seq,struct,m):
  forward.resetCache()
  backward.resetCache()
  n = len(seq)
  print "  Forward: \t",forward(seq,struct,(0,n-1),('A','G'),m)
  res = 0
  i = n-1
  for j in BASES:
    d = delta(seq,i,j)
    res += backward(seq,struct,(i-1,i+1),(j,j),m-d)
  print "  Backward of nuc %s:\t" % i, res

def test():
  seq = "UAUAUAUGUAAUACAACAAACAAUAAAAGGUGAUGCGAAUGACAAAAAUUUGGUGUACAAGCUGAUUUUUUCUCAGCUUUGUCAUAUGCCGGCAAACGUCUGGUGAAGGAUUUUAGUGUCAAAGGACAGAUGGUUAGCUUACAGGAGGAAGGCGAGGAGUAUACGCCGCGUCAGGGACUAUAUGCUAACUUAAGUUCUAGGGAGACGGAAUUUCGUUACGUCGCUCAUUUAUCAAUAUUAUCAUCCUCACUACAUUGGUAAAUAAGCGAAUAAAGUAUCGUAUCAUUUAUGGCCACUAUUUAAACAAGUUACGGGCGCUCAUUUAUUCGCAAGUUUGAAUAUCUUGAGAAGUAAGUAGUAAGAUAAAUAAUUCGCCUACUUGUUUGAAAAUGCCUCACCGUUUUCCGAAUGUUGUAUGUUUAUUCAGAAACAUCCAGACAUGGUCCGGCCCAUCAGAUGAGUGGCAAGACAAGCUCAUCCGAAAAGAAAAACCUCGUGUGACGAAAUCG"
  dbn = ".............................((.((((((..........(((((((....((((((.......))))))(((((...(((((.(...(((((((((..((((....(((((..(((((....(((((((.((.((......((((.((.(....).)).)))).....)).)).)))))))...))))).(((((((((.....((.(((.(((((.((((((((((.................)))))))))).))))).....))).))..........(((.....(((((((((((...(((((....((((..(....)..))))(((((((..(......)..))))))).......))))))))))))))))...)))...)))))))))......((.((((((......)))))).))))))).))))....))))))))).)))))).)))))...)).)))))..........)))))).))......."

  seq = seq*4
  dbn = dbn*4
  print len(seq)
  struct = parseStruct(dbn)

  m = 10

  testSingleSequence(seq,struct,m)



if __name__ == "__main__":
  test()

