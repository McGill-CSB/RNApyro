import os,sys
import itertools
import math

sys.setrecursionlimit(100000)

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
                        ('C', 'C', 'G', 'G'):-3.3,
                        ('C', 'G', 'C', 'G'):-2.4,
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
                        ('U', 'U', 'G', 'G'):-0.5})

ISO = {((k1,k2),(k3,k4)):sys.maxint for (k1,k2,k3,k4) in 
                             itertools.product(BASES, repeat=4)}
ISO.update({(('A', 'U'), ('A', 'U')): 0.0,
            (('A', 'U'), ('C', 'G')): 0.34,
            (('A', 'U'), ('G', 'C')): 0.21,
            (('A', 'U'), ('G', 'U')): 2.11,
            (('A', 'U'), ('U', 'A')): 0.31,
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

def energy(seq,struct,(a,b),(a2,b2),(i,j),alpha=1.0):
  #stacking energy, so if not stacked or i == 0 or j+1 out of len, return 1
  return 1
  if i == 0 or j == len(seq)-1 or struct[i-1] != j+1 :
    return 1
  E = STACKING_ENERGY[a,a2,b2,b]
  iso = ISO[(seq[i-1],seq[j+1]),(a,b)]+ISO[(seq[i],seq[j]),(a2,b2)]
  return  math.exp(-((alpha*E)+(1-alpha)*iso)/(BOLTZMANN*T))

@memoize
def forward(seq,struct,(i,j),(a,b),m, alpha=1.0):
  result = 0.
  if m<0: return 0
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
    elif k>i:
      for a2 in BASES:
        for b2 in BASES:
          d = delta(seq,i,a2)+delta(seq,k,b2)
          #if 'k' in middle
          if k < j:
            for m2 in range(m-d+1):
              result += forward(seq,struct,
                                  (i+1,k-1),
                                  (a2,b2),
                                  m2)*forward(seq,struct,
                                  (k+1,j),
                                  (b2,b),
                                  m-m2-d)   
          else :
            result += forward(seq,struct,
                              (i+1,k-1),
                              (a2,b2),
                              m-d)*energy(seq,struct,
                                         (a,b),
                                         (a2,b2),
                                         (i,j),
                                         alpha)
  return result

def getNext(struct,j):
  for s in range(j,len(struct)):
    if struct[s]!=-1 and struct[s]<j:
      return s
  return len(struct)

def getPrevious(struct,j):
  for s in range(j,-1,-1):
    if struct[s]!=-1 and struct[s]>j:
      return s
  return -1

def getAlts(seq,i):
  if i<0 or i>= len(seq):
    return [BASES[0]]
  else:
    return BASES

BP_CONTRIBS = {
  }

"""
def backward(seq,struct,(i,j),(a,b),m):
  result = 0.
  if i<0:
    if m==0:
      result=1.
    elif m>0:
      result=0.
  else:
    k = struct[i]
    if k==j:
      ni = getPrevious(struct,i-1)
      nj = getNext(struct,k+1)
      for a2 in getAlts(seq,ni):
        for b2 in getAlts(seq,nj):
          d = delta(seq,ni,a2) + delta(seq,nj,b2)
          for m2 in range(m-d+1):
            for m3 in range(m-m2-d+1):
              b = backward(seq,struct,
                                  (ni,nj),
                                  (a2,b2),
                                  m-m2-m3-d)
              t1 = forward(seq,struct,
                                  (ni+1,i-1),
                                  (a2,a),
                                  m2)
              t2 = forward(seq,struct,
                                  (j+1,nj-1),
                                  (b,b2),
                                  m3)
              BF = 1.0
              if (i,j)==(ni+1,nj-1):
                
              result += b*t1*t2
    else:
      print "Error! i and j should be base-paired. (i,j,part(i))=",(i,j,k)
      return -sys.maxint
  #print "B",i,j,m, result
  return result
"""

@memoize
def backward(seq,struct,(i,j),(a,b),m):
  result = 0.
  if m<0: return 0
  if i-1<0:
    result=forward(seq,struct,
                        (j+1,len(seq)-1),
                        (b,'X'),
                         m)
  else:
    k = struct[i-1]
    if k==-1:
      for a2 in getAlts(seq,i-1):
        d = delta(seq,i-1,a2)
        b = backward(seq,struct,
                        (i-1,j),
                        (a2,b),
                         m-d)
        result += b
    #BP to the left
    elif k<i-1:
      for a2 in getAlts(seq,k):
        for b2 in getAlts(seq,i-1):
          d = delta(seq,i-1,b2) + delta(seq,k,a2)
          for m2 in range(m-d+1):
            b = backward(seq,struct,
                                  (k,j),
                                  (a2,b),
                                  m-m2-d)
            f = forward(seq,struct,
                                  (k+1,i-2),
                                  (a2,b2),
                                  m2)
            result += b*f
            #print "  L",i,j,m,"B:",m-m2-d,b,"F:(",(k+1,i-1),")",m2,f
    #BP to the right
    else:
      for a2 in getAlts(seq,i-1):
        for b2 in getAlts(seq,k):
          d = delta(seq,i-1,a2) + delta(seq,k,b2)
          for m2 in range(m-d+1):
              b = backward(seq,struct,
                                  (i-1,k),
                                  (a2,b2),
                                  m-m2-d)
              f = forward(seq,struct,
                                  (j+1,k-1),
                                  (b,b2),
                                  m2)
              result += b*f
  #print "B",i,j,m, result
  return result

def validate_struct(seq,struct):
  val_bp = (('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G'))
  for i, first in enumerate(seq):
    if struct[i] == -1 or struct[i] < i:
      continue
    if (first, seq[struct[i]]) not in val_bp:
      return False
  return True


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

def expandStruct(struct):
  result = ["." for c in struct]
  for i in range(len(struct)):
    if struct[i]>i:
      result[i]="("
      result[struct[i]]=")"
  return "".join(result)

def allSecStr(n):
  if n==0:
    return [""]
  else:
    result = []
    s1 = allSecStr(n-1)
    for a in s1:
      result.append("."+a)
    for i in range(0,n-2+1):
      s1 = allSecStr(n-2-i)
      s2 = allSecStr(i)
      for a in s1:
        for b in s2:
          result.append("("+a+")"+b)
    return result

def testSingleSequence(seq,struct,m):
  forward.resetCache()
  backward.resetCache()
  n = len(seq)
  print "  Forward: \t",forward(seq,struct,(0,n-1),('A','G'),m)
  """
  print "  Backward: "
  for x in range(n):
    tot = 0.
    if struct[x]>x:
      (i,j) = (x,struct[x])
      for a2 in getAlts(seq,i):
        for b2 in getAlts(seq,j):
          d = delta(seq,i,a2) + delta(seq,j,b2)
          for m2 in range(0,m-d+1):
                tot += backward(seq,struct,(i,j),(a2,b2),m-m2-d)*forward(
                  seq,struct,(i+1,j-1),(a2,b2),m2)
      print "    ",(i,j),"P->\t",tot
      #if tot!= forward(seq,struct,(0,n-1),('A','G'),m):
      #  sys.exit()
    elif struct[x]==-1:
      (i,j) = (x,x)
      for a2 in getAlts(seq,i):
        d = delta(seq,i,a2)
        tot += backward(seq,struct,(i,j),(a2,a2),m-d)
      print "    ",(i,j),"U->\t",tot
      #if tot!= forward(seq,struct,(0,n-1),('A','G'),m):
      #  sys.exit()
      """
  


def test(n,m):
  seq = "".join(["G" for i in range(n)])
  secStrs = allSecStr(n)
  for dbn in secStrs:
    struct = parseStruct(dbn)
    print dbn
    testSingleSequence(seq,struct,m)
      


if __name__ == "__main__":
  #try:
   # n = int(sys.argv[1])
   # m = int(sys.argv[2])
    #test(n,m)
  #except ValueError, e:
    #seq = sys.argv[1]
    #dbn = sys.argv[2]
    #m = int(sys.argv[3])
  seq = "GCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACUGCUAGUCUGCGAUCGCAUCGACU"
#  dbn = "((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))((((...)))((...)))((.))"
  dbn = "...((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((((.(((((((((((((((..(((.....)))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))).)))))))))))))))..))))"
  struct = parseStruct(dbn)

  #print getBPs(seq,struct,m)
  if not validate_struct(seq, struct):
    print 'The secondary structure contains non watson-crick base pairs'
  m =  15
  testSingleSequence(seq,struct,m)





