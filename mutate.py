# mutate.py
# BIMM 185, Spring 2014
# Max Shen A10082759
#
# python mutate.py %identity len1 len2 %mismatch
#
# Randomly generates a DNA sequence of length len1. Copies it into a 2nd sequence.
# Both sequences are mutated until their % identity is as specified in the parameter.
# %mismatch determines the number of mutations that are mismatches as opposed to indels.
# The two sequences are then surrounded on both sides by 'len2' random nucleotides
#
# Writes to stdout by default in proper fasta format

import sys
import copy
import string
import datetime
import random
import math
import numpy as np
from collections import defaultdict

def main():
  if len(sys.argv) != 5:
    print 'Usage: python mutate.py <pct identity> <len 1> <len2> <pct mismatch>'
    sys.exit(0)

  identity = float(sys.argv[1])
  len1 = int(sys.argv[2])
  len2 = int(sys.argv[3])
  pctmm = float(sys.argv[4])

  mutate(identity, len1, len2, pctmm)

def mutate(identity, len1, len2, pctmm):
  seq1 = randomNT(len1)
  seq2 = copy.copy(seq1)
  seqs = [list(seq1), list(seq2)]

  numMutations = int(round(len1*(1-identity)))
  # print 'Inducing', numMutations, 'mutations:', 1 - float(numMutations)/float(len1), identity

  for i in range(numMutations):
    seq = seqs[np.random.randint(0, high = 2)]
    pos = np.random.randint(0, high = len(seq))
    if pctmm >= np.random.random():
      seq = makeMM(seq, pos)
    else:
      seq = makeID(seq, pos)

  seq1 = ''.join(seqs[0])
  seq2 = ''.join(seqs[1])

  # print seq1, '\n', seq2

  seq1 = randomNT(len2) + seq1 + randomNT(len2)
  seq2 = randomNT(len2) + seq2 + randomNT(len2)

  print '> seq1\n', seq1, '\n> seq2\n', seq2


def makeMM(seq, pos):
  lib = ['A', 'C', 'T', 'G']
  old = copy.copy(seq[pos])
  while seq[pos] == old:
    seq[pos] = lib[np.random.randint(0, high = 4)]
  return seq

def makeID(seq, pos):
  if np.random.randint(0, high = 2) == 0:
    seq.insert(pos, seq[pos])
  else:
    del seq[pos]
  return seq

def randomNT(len):
  nts = ''
  lib = ['A', 'C', 'T', 'G']
  for i in range(len):
    nts += lib[np.random.randint(0, high = 4)]
  return nts

# main
if __name__ == '__main__':
  start = datetime.datetime.now()
  # print 'Start:', start
  main()
  end = datetime.datetime.now()
  # print '\nEnd:', end, '\nTotal:', end - start
