'RNApyro 2.0, for python3'
###############################################################################
#Copyright (c) 2013, Vladimir Reinharz, Yann Ponty & Jerome Waldispuhl        #
#All rights reserved.                                                         #
#                                                                             #
#Redistribution and use in source and binary forms, with or without           #
#modification, are permitted provided that the following conditions are met:  #
#* Redistributions of source code must retain the above copyright             #
#notice, this list of conditions and the following disclaimer.                #
#* Redistributions in binary form must reproduce the above copyright          #
#notice, this list of conditions and the following disclaimer in the          #
#documentation and/or other materials provided with the distribution.         #
#* Neither the name of the <organization> nor the                             #
#names of its contributors may be used to endorse or promote products         #
#derived from this software without specific prior written permission.        #
#                                                                             #
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  #
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    #
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   #
#ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY       #
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   #
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  #
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT   #
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS#
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE                  #
###############################################################################
import os,sys
import itertools
import random
import math
import time
from decimal import Decimal

if sys.version_info.major != 3:
    print("You need to run this software with python 3")
    return 1


def MPMATH_MISSING():
  print("""The module `mpmath` was not found. This might impeed the
  processing of long RNAs.
  More information can be found on http://code.google.com/p/mpmath/
  We also recommand installation of `gmpy` (http://code.google.com/p/gmpy/),
  automatically leveraged by `mpmath` to increase the speed of computations""")

try: #For infinite precision
  from mpmath import mpf
except ImportError:
  MPMATH_MISSING()
  def mpf(n):
    return n

sys.setrecursionlimit(10000)

BASES = ['A','C','G','U']
BOLTZMANN = Decimal('0.0019872041')
T = Decimal('310.15')

#control STACKING_ENERGY with get(key,STACKING_ENERGY_MALUS)
global STACKING_ENERGY_MALUS
STACKING_ENERGY_MALUS = Decimal(sys.maxsize)
#turner04
#The order of the nucleotides is from 5' -> 3'
STACKING_ENERGY  = {('A', 'A', 'U', 'U'):Decimal('-0.9'),
                    ('A', 'C', 'G', 'U'):Decimal('-2.2'),
                    ('A', 'G', 'C', 'U'):Decimal('-2.1'),
                    ('A', 'G', 'U', 'U'):Decimal('-0.6'),
                    ('A', 'U', 'A', 'U'):Decimal('-1.1'),
                    ('A', 'U', 'G', 'U'):Decimal('-1.4'),
                    ('C', 'A', 'U', 'G'):Decimal('-2.1'),
                    ('C', 'G', 'U', 'G'):Decimal('-1.4'),
                    ('C', 'U', 'A', 'G'):Decimal('-2.1'),
                    ('C', 'U', 'G', 'G'):Decimal('-2.1'),
                    ('G', 'A', 'U', 'C'):Decimal('-2.4'),
                    ('G', 'C', 'G', 'C'):Decimal('-3.4'),
                    ('G', 'G', 'C', 'C'):Decimal('-3.3'),
                    ('G', 'G', 'U', 'C'):Decimal('-1.5'),
                    ('G', 'U', 'A', 'C'):Decimal('-2.2'),
                    ('G', 'U', 'G', 'C'):Decimal('-2.5'),
                    ('G', 'A', 'U', 'U'):Decimal('-1.3'),
                    ('G', 'C', 'G', 'U'):Decimal('-2.5'),
                    ('G', 'G', 'C', 'U'):Decimal('-2.1'),
                    ('G', 'G', 'U', 'U'):Decimal('-0.5'),
                    ('G', 'U', 'A', 'U'):Decimal('-1.4'),
                    ('G', 'U', 'G', 'U'):Decimal('1.3'),
                    ('U', 'A', 'U', 'A'):Decimal('-1.3'),
                    ('U', 'C', 'G', 'A'):Decimal('-2.4'),
                    ('U', 'G', 'C', 'A'):Decimal('-2.1'),
                    ('U', 'G', 'U', 'A'):Decimal('-1.0'),
                    ('U', 'U', 'A', 'A'):Decimal('-0.9'),
                    ('U', 'U', 'G', 'A'):Decimal('-1.3'),
                    ('U', 'A', 'U', 'G'):Decimal('-1.0'),
                    ('U', 'C', 'G', 'G'):Decimal('-1.5'),
                    ('U', 'G', 'C', 'G'):Decimal('-1.4'),
                    ('U', 'G', 'U', 'G'):Decimal('0.3'),
                    ('U', 'U', 'A', 'G'):Decimal('-0.6'),
                    ('U', 'U', 'G', 'G'):Decimal('-0.5'),
                    ('C', 'C', 'G', 'G'):Decimal('-3.3'),
                    ('C', 'G', 'C', 'G'):Decimal('-2.4')}

#All isostericity values. GG->anything or anything->GG is 10
ISO = {(('A', 'C'), ('C', 'A')):Decimal('4.93'),
        (('G', 'C'), ('A', 'G')):Decimal('3.5'),
        (('G', 'C'), ('G', 'C')):Decimal('0.0'),
        (('U', 'G'), ('U', 'G')):Decimal('0.0'),
        (('A', 'C'), ('G', 'U')):Decimal('0.8'),
        (('A', 'G'), ('U', 'U')):Decimal('8.18'),
        (('G', 'A'), ('A', 'G')):Decimal('2.25'),
        (('A', 'A'), ('A', 'C')):Decimal('4.58'),
        (('A', 'C'), ('U', 'G')):Decimal('4.76'),
        (('U', 'G'), ('C', 'G')):Decimal('2.14'),
        (('U', 'G'), ('G', 'G')):Decimal('10.0'),
        (('U', 'C'), ('U', 'U')):Decimal('4.36'),
        (('G', 'G'), ('U', 'A')):Decimal('10.0'),
        (('C', 'G'), ('G', 'C')):Decimal('0.26'),
        (('A', 'U'), ('C', 'A')):Decimal('2.75'),
        (('A', 'U'), ('U', 'U')):Decimal('3.8'),
        (('U', 'U'), ('G', 'C')):Decimal('3.94'),
        (('U', 'G'), ('G', 'C')):Decimal('2.39'),
        (('C', 'A'), ('C', 'A')):Decimal('0.0'),
        (('G', 'A'), ('U', 'G')):Decimal('3.8'),
        (('A', 'G'), ('G', 'C')):Decimal('4.38'),
        (('A', 'A'), ('C', 'G')):Decimal('3.44'),
        (('G', 'G'), ('C', 'U')):Decimal('10.0'),
        (('G', 'C'), ('U', 'A')):Decimal('0.34'),
        (('U', 'A'), ('A', 'A')):Decimal('4.66'),
        (('U', 'U'), ('U', 'C')):Decimal('6.46'),
        (('U', 'U'), ('C', 'G')):Decimal('3.8'),
        (('C', 'C'), ('U', 'C')):Decimal('8.25'),
        (('G', 'A'), ('C', 'U')):Decimal('0.0'),
        (('C', 'U'), ('G', 'U')):Decimal('5.06'),
        (('C', 'G'), ('U', 'C')):Decimal('3.44'),
        (('U', 'C'), ('G', 'U')):Decimal('2.89'),
        (('U', 'U'), ('U', 'G')):Decimal('2.89'),
        (('G', 'A'), ('U', 'U')):Decimal('5.97'),
        (('G', 'U'), ('A', 'G')):Decimal('4.44'),
        (('A', 'C'), ('G', 'A')):Decimal('5.05'),
        (('G', 'C'), ('U', 'G')):Decimal('2.39'),
        (('G', 'G'), ('A', 'A')):Decimal('10.0'),
        (('U', 'U'), ('C', 'A')):Decimal('2.39'),
        (('G', 'G'), ('C', 'A')):Decimal('10.0'),
        (('A', 'A'), ('C', 'U')):Decimal('1.53'),
        (('G', 'U'), ('C', 'C')):Decimal('6.25'),
        (('C', 'C'), ('U', 'U')):Decimal('2.37'),
        (('C', 'U'), ('U', 'C')):Decimal('7.97'),
        (('C', 'C'), ('G', 'G')):Decimal('10.0'),
        (('G', 'G'), ('A', 'U')):Decimal('10.0'),
        (('C', 'G'), ('A', 'U')):Decimal('0.34'),
        (('C', 'U'), ('A', 'U')):Decimal('5.3'),
        (('U', 'A'), ('C', 'A')):Decimal('2.47'),
        (('U', 'U'), ('A', 'U')):Decimal('3.8'),
        (('G', 'G'), ('U', 'U')):Decimal('10.0'),
        (('U', 'A'), ('C', 'C')):Decimal('5.3'),
        (('G', 'G'), ('A', 'G')):Decimal('10.0'),
        (('A', 'G'), ('A', 'U')):Decimal('4.52'),
        (('A', 'A'), ('A', 'G')):Decimal('2.33'),
        (('U', 'C'), ('U', 'A')):Decimal('3.8'),
        (('G', 'U'), ('U', 'U')):Decimal('5.27'),
        (('A', 'A'), ('C', 'A')):Decimal('5.3'),
        (('U', 'G'), ('C', 'U')):Decimal('3.8'),
        (('G', 'U'), ('G', 'U')):Decimal('0.0'),
        (('G', 'A'), ('A', 'U')):Decimal('3.57'),
        (('U', 'A'), ('A', 'U')):Decimal('0.31'),
        (('U', 'A'), ('U', 'C')):Decimal('3.57'),
        (('U', 'G'), ('A', 'G')):Decimal('4.33'),
        (('C', 'G'), ('C', 'A')):Decimal('2.55'),
        (('U', 'A'), ('U', 'A')):Decimal('0.0'),
        (('G', 'G'), ('G', 'A')):Decimal('10.0'),
        (('U', 'U'), ('C', 'U')):Decimal('5.97'),
        (('A', 'C'), ('A', 'G')):Decimal('5.14'),
        (('U', 'C'), ('A', 'A')):Decimal('6.91'),
        (('A', 'C'), ('A', 'U')):Decimal('2.47'),
        (('G', 'C'), ('A', 'C')):Decimal('2.55'),
        (('C', 'C'), ('A', 'C')):Decimal('5.91'),
        (('A', 'U'), ('A', 'U')):Decimal('0.0'),
        (('C', 'C'), ('C', 'C')):Decimal('0.0'),
        (('U', 'C'), ('A', 'C')):Decimal('2.39'),
        (('G', 'A'), ('C', 'C')):Decimal('7.97'),
        (('A', 'A'), ('G', 'C')):Decimal('3.39'),
        (('U', 'C'), ('A', 'G')):Decimal('6.96'),
        (('G', 'C'), ('G', 'A')):Decimal('3.49'),
        (('U', 'U'), ('U', 'A')):Decimal('3.63'),
        (('G', 'U'), ('G', 'G')):Decimal('10.0'),
        (('C', 'G'), ('U', 'G')):Decimal('2.14'),
        (('A', 'G'), ('U', 'A')):Decimal('4.66'),
        (('A', 'A'), ('U', 'A')):Decimal('3.57'),
        (('U', 'G'), ('U', 'A')):Decimal('2.11'),
        (('C', 'A'), ('U', 'C')):Decimal('5.3'),
        (('G', 'A'), ('C', 'G')):Decimal('3.39'),
        (('U', 'A'), ('C', 'G')):Decimal('0.21'),
        (('A', 'U'), ('C', 'U')):Decimal('3.57'),
        (('G', 'C'), ('U', 'U')):Decimal('3.94'),
        (('U', 'U'), ('C', 'C')):Decimal('2.37'),
        (('G', 'U'), ('A', 'A')):Decimal('4.1'),
        (('A', 'C'), ('C', 'U')):Decimal('5.3'),
        (('A', 'G'), ('C', 'C')):Decimal('9.77'),
        (('U', 'A'), ('G', 'C')):Decimal('0.34'),
        (('A', 'U'), ('A', 'C')):Decimal('2.47'),
        (('U', 'G'), ('G', 'A')):Decimal('4.44'),
        (('G', 'U'), ('A', 'U')):Decimal('2.11'),
        (('A', 'U'), ('A', 'A')):Decimal('4.52'),
        (('C', 'U'), ('C', 'C')):Decimal('3.02'),
        (('C', 'U'), ('U', 'A')):Decimal('5.39'),
        (('A', 'C'), ('U', 'C')):Decimal('4.58'),
        (('C', 'G'), ('G', 'A')):Decimal('3.5'),
        (('C', 'A'), ('G', 'G')):Decimal('10.0'),
        (('G', 'A'), ('G', 'G')):Decimal('10.0'),
        (('C', 'G'), ('U', 'U')):Decimal('3.8'),
        (('U', 'A'), ('A', 'C')):Decimal('2.75'),
        (('G', 'A'), ('C', 'A')):Decimal('4.58'),
        (('G', 'A'), ('A', 'C')):Decimal('5.3'),
        (('A', 'C'), ('C', 'G')):Decimal('2.78'),
        (('U', 'C'), ('G', 'G')):Decimal('10.0'),
        (('A', 'A'), ('A', 'U')):Decimal('3.5'),
        (('C', 'C'), ('C', 'U')):Decimal('7.97'),
        (('A', 'G'), ('U', 'C')):Decimal('2.71'),
        (('A', 'C'), ('A', 'C')):Decimal('0.0'),
        (('U', 'A'), ('G', 'A')):Decimal('3.67'),
        (('C', 'C'), ('G', 'U')):Decimal('6.25'),
        (('A', 'U'), ('U', 'G')):Decimal('2.4'),
        (('U', 'C'), ('C', 'G')):Decimal('3.94'),
        (('U', 'C'), ('U', 'C')):Decimal('5.97'),
        (('A', 'G'), ('A', 'A')):Decimal('0.0'),
        (('A', 'A'), ('U', 'G')):Decimal('4.59'),
        (('C', 'C'), ('A', 'U')):Decimal('5.39'),
        (('A', 'G'), ('C', 'U')):Decimal('3.77'),
        (('C', 'C'), ('U', 'G')):Decimal('5.06'),
        (('G', 'U'), ('C', 'G')):Decimal('2.39'),
        (('G', 'C'), ('A', 'A')):Decimal('4.38'),
        (('U', 'C'), ('U', 'G')):Decimal('5.27'),
        (('U', 'G'), ('C', 'A')):Decimal('0.8'),
        (('A', 'A'), ('G', 'U')):Decimal('3.8'),
        (('C', 'U'), ('U', 'U')):Decimal('4.31'),
        (('C', 'C'), ('G', 'C')):Decimal('5.56'),
        (('G', 'C'), ('G', 'U')):Decimal('2.14'),
        (('U', 'U'), ('U', 'U')):Decimal('0.0'),
        (('U', 'C'), ('G', 'A')):Decimal('6.91'),
        (('C', 'U'), ('C', 'U')):Decimal('8.25'),
        (('C', 'G'), ('A', 'C')):Decimal('2.78'),
        (('G', 'U'), ('U', 'C')):Decimal('3.8'),
        (('A', 'G'), ('C', 'G')):Decimal('4.5'),
        (('G', 'G'), ('U', 'C')):Decimal('10.0'),
        (('U', 'G'), ('U', 'C')):Decimal('4.59'),
        (('C', 'A'), ('C', 'U')):Decimal('4.58'),
        (('G', 'A'), ('U', 'C')):Decimal('1.53'),
        (('G', 'G'), ('C', 'G')):Decimal('10.0'),
        (('A', 'U'), ('G', 'U')):Decimal('2.11'),
        (('G', 'U'), ('C', 'A')):Decimal('4.76'),
        (('A', 'C'), ('U', 'U')):Decimal('5.21'),
        (('A', 'G'), ('C', 'A')):Decimal('6.7'),
        (('C', 'A'), ('U', 'A')):Decimal('2.47'),
        (('A', 'G'), ('U', 'G')):Decimal('6.07'),
        (('C', 'G'), ('U', 'A')):Decimal('0.21'),
        (('U', 'G'), ('A', 'U')):Decimal('2.4'),
        (('G', 'G'), ('C', 'C')):Decimal('10.0'),
        (('C', 'G'), ('A', 'A')):Decimal('4.5'),
        (('C', 'U'), ('C', 'A')):Decimal('5.91'),
        (('G', 'U'), ('U', 'G')):Decimal('4.48'),
        (('G', 'C'), ('C', 'A')):Decimal('2.78'),
        (('A', 'A'), ('U', 'C')):Decimal('0.0'),
        (('G', 'G'), ('G', 'C')):Decimal('10.0'),
        (('A', 'U'), ('G', 'C')):Decimal('0.21'),
        (('C', 'A'), ('C', 'G')):Decimal('2.55'),
        (('U', 'C'), ('C', 'A')):Decimal('5.21'),
        (('U', 'G'), ('U', 'U')):Decimal('2.89'),
        (('G', 'G'), ('G', 'U')):Decimal('10.0'),
        (('C', 'A'), ('A', 'G')):Decimal('5.05'),
        (('A', 'C'), ('U', 'A')):Decimal('2.75'),
        (('A', 'U'), ('C', 'G')):Decimal('0.34'),
        (('U', 'G'), ('C', 'C')):Decimal('5.06'),
        (('A', 'C'), ('G', 'C')):Decimal('2.55'),
        (('G', 'U'), ('A', 'C')):Decimal('0.8'),
        (('C', 'A'), ('C', 'C')):Decimal('4.49'),
        (('A', 'U'), ('G', 'G')):Decimal('10.0'),
        (('A', 'A'), ('U', 'U')):Decimal('6.46'),
        (('G', 'G'), ('G', 'G')):Decimal('0.0'),
        (('A', 'C'), ('G', 'G')):Decimal('10.0'),
        (('A', 'C'), ('A', 'A')):Decimal('4.8'),
        (('A', 'A'), ('G', 'G')):Decimal('10.0'),
        (('G', 'C'), ('U', 'C')):Decimal('3.39'),
        (('U', 'C'), ('G', 'C')):Decimal('3.8'),
        (('U', 'A'), ('G', 'U')):Decimal('2.4'),
        (('G', 'A'), ('G', 'C')):Decimal('3.44'),
        (('G', 'G'), ('A', 'C')):Decimal('10.0'),
        (('A', 'A'), ('A', 'A')):Decimal('2.71'),
        (('C', 'C'), ('A', 'G')):Decimal('8.82'),
        (('C', 'G'), ('C', 'C')):Decimal('5.49'),
        (('U', 'A'), ('U', 'G')):Decimal('2.11'),
        (('U', 'G'), ('A', 'C')):Decimal('4.76'),
        (('C', 'U'), ('U', 'G')):Decimal('6.25'),
        (('C', 'C'), ('C', 'G')):Decimal('5.49'),
        (('G', 'U'), ('G', 'A')):Decimal('4.33'),
        (('G', 'C'), ('C', 'G')):Decimal('0.26'),
        (('A', 'G'), ('A', 'G')):Decimal('2.41'),
        (('A', 'A'), ('C', 'C')):Decimal('8.25'),
        (('A', 'G'), ('G', 'A')):Decimal('2.18'),
        (('C', 'G'), ('G', 'U')):Decimal('2.39'),
        (('A', 'C'), ('C', 'C')):Decimal('5.91'),
        (('A', 'G'), ('G', 'G')):Decimal('10.0'),
        (('U', 'U'), ('G', 'A')):Decimal('6.96'),
        (('U', 'C'), ('C', 'U')):Decimal('6.46'),
        (('A', 'U'), ('U', 'A')):Decimal('0.31'),
        (('G', 'A'), ('G', 'A')):Decimal('2.33'),
        (('G', 'C'), ('C', 'U')):Decimal('3.44'),
        (('G', 'A'), ('U', 'A')):Decimal('3.5'),
        (('U', 'G'), ('A', 'A')):Decimal('6.07'),
        (('C', 'U'), ('A', 'G')):Decimal('8.86'),
        (('U', 'U'), ('G', 'U')):Decimal('5.27'),
        (('A', 'U'), ('A', 'G')):Decimal('3.67'),
        (('G', 'C'), ('C', 'C')):Decimal('5.56'),
        (('A', 'U'), ('C', 'C')):Decimal('5.39'),
        (('U', 'C'), ('A', 'U')):Decimal('3.63'),
        (('C', 'A'), ('G', 'C')):Decimal('2.78,')
        (('C', 'G'), ('C', 'G')):Decimal('0.0,')
        (('C', 'C'), ('A', 'A')):Decimal('9.77,')
        (('G', 'C'), ('A', 'U')):Decimal('0.21,')
        (('A', 'U'), ('G', 'A')):Decimal('3.67,')
        (('U', 'U'), ('A', 'A')):Decimal('8.18,')
        (('U', 'C'), ('C', 'C')):Decimal('4.31,')
        (('C', 'U'), ('G', 'A')):Decimal('8.82,')
        (('C', 'A'), ('U', 'G')):Decimal('0.8,')
        (('C', 'U'), ('A', 'A')):Decimal('9.05,')
        (('C', 'A'), ('A', 'A')):Decimal('6.7,')
        (('C', 'C'), ('U', 'A')):Decimal('5.3,')
        (('G', 'U'), ('C', 'U')):Decimal('4.59,')
        (('C', 'A'), ('A', 'C')):Decimal('4.93,')
        (('C', 'A'), ('G', 'U')):Decimal('4.76,')
        (('G', 'C'), ('G', 'G')):Decimal('10.0,')
        (('U', 'A'), ('G', 'G')):Decimal('10.0,')
        (('G', 'A'), ('A', 'A')):Decimal('3.77,')
        (('U', 'U'), ('A', 'G')):Decimal('6.91,')
        (('C', 'G'), ('C', 'U')):Decimal('3.39,')
        (('A', 'G'), ('A', 'C')):Decimal('4.8,')
        (('U', 'U'), ('A', 'C')):Decimal('5.21,')
        (('U', 'A'), ('U', 'U')):Decimal('3.63,')
        (('G', 'U'), ('G', 'C')):Decimal('2.14,')
        (('C', 'U'), ('A', 'C')):Decimal('4.49,')
        (('C', 'A'), ('G', 'A')):Decimal('5.14,')
        (('C', 'G'), ('G', 'G')):Decimal('10.0,')
        (('C', 'C'), ('G', 'A')):Decimal('8.86,')
        (('G', 'A'), ('G', 'U')):Decimal('4.59,')
        (('U', 'A'), ('A', 'G')):Decimal('3.67,')
        (('U', 'U'), ('G', 'G')):Decimal('10.0,')
        (('U', 'A'), ('C', 'U')):Decimal('3.5,')
        (('C', 'C'), ('C', 'A')):Decimal('4.49,')
        (('C', 'U'), ('G', 'C')):Decimal('5.49,')
        (('C', 'U'), ('C', 'G')):Decimal('5.56,')
        (('A', 'A'), ('G', 'A')):Decimal('2.25,')
        (('C', 'U'), ('G', 'G')):Decimal('10.0,')
        (('C', 'A'), ('U', 'U')):Decimal('2.39,')
        (('A', 'G'), ('G', 'U')):Decimal('4.1,')
        (('G', 'U'), ('U', 'A')):Decimal('2.4,')
        (('U', 'G'), ('G', 'U')):Decimal('4.48,')
        (('C', 'A'), ('A', 'U')):Decimal('2.75,')
        (('C', 'G'), ('A', 'G')):Decimal('3.49,')
        (('G', 'G'), ('U', 'G')):Decimal('10.0,')
        (('A', 'U'), ('U', 'C')):Decimal('3.5')}


class memoize(dict):
  """Generically memoizes a function results."""
  fun = None
  
  def __init__(self, f):
      self.fun = f
  
  def __call__(self,seq,ref_seq,struct,*args):
      nargs = (args)
      if nargs in self:
          return self[nargs]
      else:
          val = mpf(self.fun(seq,ref_seq,struct,*args))
          self[nargs] = val
          return val
  def resetCache(self):
      self.clear()

class memoize_iso(dict):
  """Generically memoizes a function results."""
  fun = None
  
  def __init__(self, f):
      self.fun = f
  
  def __call__(self,*args):
      nargs = [args[0]]
      nargs.extend(args[2:])
      nargs = tuple(nargs)
      if nargs in self:
          return self[nargs]
      else:
          val = mpf(self.fun(*args))
          self[nargs] = val
          return val
  def resetCache(self):
      self.clear()


def delta(seq,i,c):
  if i<0 or i>=len(seq):
    return 0
  if c==seq[i]:
    return 0
  else:
    return 1

def energy(xxx_todo_changeme6, xxx_todo_changeme7,alpha):
  #stacking energy of base pair (a,b) around base pair (a2,b2)
  (a,b) = xxx_todo_changeme6
  (a2,b2) = xxx_todo_changeme7
  E = STACKING_ENERGY.get((a,a2,b2,b), MAX_NB)
  return  math.exp(-(alpha*E)/(BOLTZMANN*T))

def isGapChar(c):
  return c=="." or c=="-"

@memoize_iso
def isostericity(seq,ref_seq, xxx_todo_changeme, xxx_todo_changeme1, alpha):
  (i,j) = xxx_todo_changeme
  (a,b) = xxx_todo_changeme1
  if not ref_seq:
    return 1.
  #isostericity of going from original base pair to (a,b)
  iso_delta = 0.
  nb_correct = 0
  for ref in ref_seq:
    if not (isGapChar(ref[i]) or isGapChar(ref[j])):
      iso_delta += ISO[(ref[i],ref[j]),(a,b)]-ISO[(ref[i],ref[j]),(seq[i],seq[j])]
      nb_correct += 1
  if nb_correct>0:
    iso = mpf(iso_delta)/nb_correct
  else:
    iso = mpf(0.)
  return  math.exp(-((1-alpha)*iso)/(BOLTZMANN*T))

@memoize
def forward(seq,ref_seq,struct, xxx_todo_changeme2, xxx_todo_changeme3,m, alpha):
  #alpha gives the weight energy vs isostericity
  (i,j) = xxx_todo_changeme2
  (a,b) = xxx_todo_changeme3
  result = 0.
  if m<0 or m>j-i+1: return 0
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
          result += forward(seq,ref_seq,struct,
                            (i+1,j),
                            (a2,b),
                            m-d,alpha)
    elif i < k <= j: #If k > j we return 0
      for a2 in BASES:
        for b2 in BASES:
          d = delta(seq,i,a2)+delta(seq,k,b2)
          #if not stacked outside (or border, then no stack possible)
          if i==0 or j==len(seq)-1 or not (j==k and struct[i-1]==j+1):
            for m2 in range(m-d+1):
              result += forward(seq,ref_seq,struct,
                                  (i+1,k-1),
                                  (a2,b2),
                                  m2,alpha)*forward(seq,ref_seq,struct,
                                  (k+1,j),
                                  (b2,b),
                                  m-m2-d,alpha)*isostericity(seq,ref_seq,
                                                       (i,k),
                                                       (a2,b2),
                                                       alpha)

          #if stack, we add energy
          else :
            result += forward(seq,ref_seq,struct,
                              (i+1,k-1),
                              (a2,b2),
                              m-d,alpha)*energy((a,b),
                                          (a2,b2),
                                          alpha)*isostericity(seq,ref_seq,
                                                              (i,k),
                                                              (a2,b2),
                                                              alpha)
  return result

def random_weighted_sampling(l_samples):
  tot = sum(x[1] for x in l_samples)
  if tot == 0:
    print(l_samples)
    return random.choice([x[0] for x in l_samples])
  scaled_weights = [x[1]/tot for x in l_samples]
  rand_nb = random.random()
  accumulation = 0
  for i,x in enumerate(scaled_weights):
    accumulation += x
    if accumulation > rand_nb:
      return l_samples[i][0]
  return l_samples[-1][0]

def backtrack(seq,ref_seq,struct, xxx_todo_changeme8, xxx_todo_changeme9,m, alpha):
  #alpha gives the weight energy vs isostericity
  (i,j) = xxx_todo_changeme8
  (a,b) = xxx_todo_changeme9
  result = 0.
  max_seq = ''
  if m==0: return seq[i:j+1]
  if i > j : return ''
  else:
    k = struct[i]
    if k==-1:
      l_samples = []
      for a2 in BASES:
        d = delta(seq,i,a2)
        result = forward(seq,ref_seq,struct,
                          (i+1,j),
                          (a2,b),
                          m-d,alpha)
        l_samples.append((a2,result))
      a2 = random_weighted_sampling(l_samples)
      d = delta(seq,i,a2)
      max_seq = a2 + backtrack(seq,ref_seq,struct,(i+1,j),(a2,b),m-d,alpha)
    elif i < k <= j: #If k > j we return 0
      l_samples=[]
      for a2 in BASES:
        for b2 in BASES:
          d = delta(seq,i,a2)+delta(seq,k,b2)
          #if not stacked outside (or border, then no stack possible)
          if i==0 or j==len(seq)-1 or not (j==k and struct[i-1]==j+1):
            for m2 in range(m-d+1):
              result = forward(seq,ref_seq,struct,
                                  (i+1,k-1),
                                  (a2,b2),
                                  m2,alpha)*forward(seq,ref_seq,struct,
                                  (k+1,j),
                                  (b2,b),
                                  m-m2-d,alpha)*isostericity(seq,ref_seq,
                                                       (i,k),
                                                       (a2,b2),
                                                       alpha)
              l_samples.append(((a2,b2,m2),result))
          #if stack, we add energy
          else :
            result = forward(seq,ref_seq,struct,
                              (i+1,k-1),
                              (a2,b2),
                              m-d,alpha)*energy((a,b),
                                          (a2,b2),
                                          alpha)*isostericity(seq,ref_seq,
                                                              (i,k),
                                                              (a2,b2),
                                                              alpha)
            l_samples.append(((a2,b2,m-d),result))
      a2,b2,m2 = random_weighted_sampling(l_samples)
      d = delta(seq,i,a2)+delta(seq,k,b2)
      best_1 = backtrack(seq,ref_seq,struct,(i+1,k-1),(a2,b2),m2,alpha)
      best_2 = backtrack(seq,ref_seq,struct,(k+1,j),(b2,b),m-m2-d,alpha)
      max_seq = a2+best_1+b2+best_2
  return max_seq

@memoize
def backward(seq,ref_seq,struct, xxx_todo_changeme4, xxx_todo_changeme5,m,alpha):
  (i,j) = xxx_todo_changeme4
  (a,b) = xxx_todo_changeme5
  result = 0.
  if m<0: return 0
  if i<0:
    result=forward(seq,ref_seq,struct,
                        (j,len(seq)-1),
                        ('X','X'),
                         m,alpha)
  else:
    k = struct[i]
    if k==-1:
      for a2 in BASES:
        d = delta(seq,i,a2)
        result += backward(seq,ref_seq,struct,
                        (i-1,j),
                        (a2,b),
                         m-d,alpha)
    #BP to the left
    elif k<i:
      for a2 in BASES:
        for b2 in BASES:
          d = delta(seq,i,b2) + delta(seq,k,a2)
          for m2 in range(m-d+1):
            result += backward(seq,ref_seq,struct,
                                  (k-1,j),
                                  (a2,b),
                                  m-m2-d,alpha)*forward(seq,ref_seq,struct,
                                                  (k+1,i-1),
                                                  (a2,b2),
                                                  m2,alpha)*isostericity(seq,ref_seq,
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
              result += backward(seq,ref_seq,struct,
                                 (i-1,k+1),
                                 (a2,b2),
                                 m-m2-d,alpha)*forward(seq,ref_seq,struct,
                                                 (j,k-1),
                                                 (b,b2),
                                                 m2,alpha)*isostericity(seq,ref_seq,
                                                                  (i,k),
                                                                  (a2,b2),
                                                                  alpha)
          else:
            result +=  backward(seq,ref_seq,struct,
                                 (i-1,k+1),
                                 (a2,b2),
                                 m-d,alpha)*energy((a2,b2),
                                             (a,b),
                                             alpha)*isostericity(seq,ref_seq,
                                                                 (i,k),
                                                                 (a2,b2),
                                                                 alpha)
  return result

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

def product_given_i_m(seq,ref_seq,struct,i,a,m,alpha):
  """Will compute the sum of boltzmann weights of structures 
  with 'm' mutations from 'seq' where the 'i-th' nucleotide is 'a'.
  """
  n = len(seq)
  tot = forward(seq,ref_seq,struct,(0,n-1),('X', 'X'),m,alpha)
  k = struct[i]
  result = mpf(0)
  if k == -1:
    d = delta(seq,i,a)
    result += backward(seq,ref_seq,struct,(i-1,i+1),(a,a),m-d,alpha)
  elif k < i:
    for c in BASES:
      d = delta(seq,i,a) + delta(seq,k,c)
      for m2 in range(m-d+1):
        f = forward(seq,ref_seq,struct,(k+1,i-1),(c,a),m-d-m2,alpha) 
        b = backward(seq,ref_seq,struct,(k-1,i+1),(c,a),m2,alpha)
        iso = isostericity(seq,ref_seq,(k,i),(c,a), alpha)
        result += f*b*iso
  else:
    for c in BASES:
      d = delta(seq,i,a) + delta(seq,k,c)
      for m2 in range(m-d+1):
        f = forward(seq,ref_seq,struct,(i+1,k-1),(a,c),m-d-m2,alpha) 
        b = backward(seq,ref_seq,struct,(i-1,k+1),(a,c),m2,alpha)
        iso = isostericity(seq,ref_seq,(i,k),(a,c), alpha)
        result += f*b*iso
  return result

def probability_given_i_m(seq,ref_seq,struct,i,a,m,alpha):
  """Will compute the probability that the 'i-th' nucleotide
  is 'a' over all sequences at 'm' mutations from seq
  """
  n = len(seq)
  tot = forward(seq,ref_seq,struct,(0,n-1),('X', 'X'),m,alpha)
  result = product_given_i_m(seq,ref_seq,struct,i,a,m,alpha)
  if tot == 0:
    print("""The total partition function is 0, you might want to increase
    the number of mutations allowed""")
    sys.exit(1)
  return result/tot

def probability_given_i_most_m(seq,ref_seq,struct,i,a,m,alpha):
  """Will compute the probability that the 'i-th' nucleotide
  is 'a' over all sequences between 0 and 'm' mutations from seq
  """
  n = len(seq)
  tot = 0
  result = 0
  for m2 in range(m+1):
    tot += forward(seq,ref_seq,struct,(0,n-1),('X', 'X'),m2,alpha)
    result += product_given_i_m(seq,ref_seq,struct,i,a,m2,alpha)
  if tot == 0:
    print("""The total partition function is 0, you might want to increase
    the number of mutations allowed""")
    sys.exit(1)

  return result/tot

def testSingleSequence(seq,ref_seq,struct,m,alpha):
  forward.resetCache()
  backward.resetCache()
  n = len(seq)

  i = 29
  #print "Given i=%s and m=%s the probability of sequence is:" %(i,m)
  #print "\t", probability_given_i_most_m(seq,ref_seq,struct,i,'C',m,alpha=1.0)
  print("  Forward: \t",forward(seq,ref_seq,struct,(0,n-1),('A','G'),m,alpha))
  i = n-1
  res = 0
  for j in BASES:
    d = delta(seq,i,j)
    res += backward(seq,ref_seq,struct,(i-1,i+1),(j,j),m-d,alpha)
  print("  Backward of nuc %s:\t" % i, res)

def test():
  
  seq = "UAUAUAUGUAAUACAACAAACAAUAAAAGGUG"
  dbn = "................................"
  ref_seq = seq
  alpha = 0.2 

  ref_seq = ref_seq*1
  seq = seq*1
  dbn = dbn*1
  struct = parseStruct(dbn)


  m = 10

  testSingleSequence(seq,ref_seq,struct,m,alpha)

def parse_fasta(file_name):
  #first seq is reference, other MSE, one struct only
  seq = []
  struct = []
  with open(file_name) as f:
    for line in f:
      line = line.strip()
      if not line:
        continue
      if all(x in 'AUGC.' for x in line):
        seq.append(line)
        continue
      if all(x in '(.)' for x in line):
        struct.append(line)
        continue
  return seq[1:],seq[0], parseStruct(struct[0])

def all_probabilities(seq,ref_seq, stuct, m, alpha):
  n = len(seq)
  results = []
  for i in range(n):
    results.append([])
    for a in BASES:
      results[-1].append(
        probability_given_i_most_m(seq,ref_seq,struct,i,a,m,alpha))
  return results

def display_all_probabilities(results):
  for i, x in enumerate(results):
    print(i, '\t'.join(str(y) for y in x))

def verify_file_name(file_name):
  if not os.path.isfile(file_name):
    help(file_name=True)
    sys.exit(1)
  return parse_fasta(file_name) 
     
def verify_mutants(m):
  try:
    m = int(m)
  except ValueError:
    help(mutants=True)
    sys.exit(1)
  if m < 0:
    help(mutants=True)
    sys.exit(1)
  return m

def verify_alpha(a):
  try:
    a = float(a)
  except ValueError:
    help(alpha=True)
    sys.exit(1)
  if a < 0 or a > 1:
    help(alpha=True)
    sys.exit(1)
  return a

def verify_penality(z):
  global STACKING_ENERGY
  try:
    z = float(opts[4])
    for x in STACKING_ENERGY:
      if STACKING_ENERGY[x] == sys.maxsize:
        STACKING_ENERGY[x] = z
  except ValueError:
    print(help(penality=True))
    sys.exit(1)
  return None

def verify_nb_backtrack(b):
  try:
    b = int(b)
  except ValueError:
    help(nb_backtrack=True)
    sys.exit(1)
  return b

def help(file_name=False,
         mutants=False,
         alpha=False,
         penality=False,
         nb_backtrack=False):

  print("""Required:
      -f <file_path>  Path to the file containing the reference
                      sequence, the MSE (optional) and the secondary
                      structure. The FIRST sequence in the file
                      will be used as reference
      -m <int>        The numbers of mutants required
      -a <float>    The weight given to match the MSE or the secondary
                    structure, between [0,1]. A weight of 1 will only
                    take into account the secondary structure, 0 only
                    the MSA
    Optional:
      -p <float>    The penality for a non canonical base pair (default
                    sys.maxint)
      -b <int>      Backtrack stochastic, number of sequences to output
      -no_profile   Do no output the resultant profile
    """)




  if file_name:
    print("""The location of the file does not exist, please enter
    a valid path\n""")
 
  if mutants:
    print("""The number of mutants should be an integer bigger than 0\n""")

  if alpha:
    print("""The weight to the MSE/SecStruct match should be a 
    float between 0 and 1\n""")

  if penality:
    print("""The penality argument should be a float\n""")

  if nb_backtrack:
    print("""The number of backtracked sequence should be an int\n""")

if __name__ == "__main__":

  f_no_profile = False
  f_nb_backtrack = 0

  opts = sys.argv
  for i,opt in enumerate(opts):
    if opt == '-f':
      ref_seq,seq,struct =  verify_file_name(opts[i+1])
    elif opt == '-m':
      m = verify_mutants(opts[i+1])
    elif opt == '-a':
      alpha = verify_alpha(opts[i+1])
    elif opt == '-p':
      verify_penality(opts[i+1])
    elif opt == '-b':
      f_nb_backtrack = verify_nb_backtrack(opts[i+1])
    elif opt == '-no_profile':
      f_no_profile = True
    elif opt == '-h':
        help()
        sys.exit(0)


  if len(seq) < m:
    print("The number of mutants is bigger than the length of the seq.")
    help()
    sys.exit(1)

  n = len(seq)
  if not f_no_profile:
    tic = time.time()
    results = all_probabilities(seq,ref_seq,struct,m,alpha)
    display_all_probabilities(results)
    print(time.time() - tic)
  else:
    forward(seq,ref_seq,struct,(0,n-1),('X', 'X'),m,alpha)

  for i in range(f_nb_backtrack):
    print(backtrack(seq,ref_seq,struct,(0,n-1),('X','X'),m, alpha))
