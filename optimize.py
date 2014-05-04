import sys
import string
import datetime
import random
import numpy as np
import locAL

from collections import defaultdict

def main():
  if len(sys.argv) != 3:
    print 'Usage: python optimize <input file> <input file 2>'
    sys.exit(0)

  # Read in seqs
  with open(sys.argv[1]) as f:
    lines = f.readlines()
    if len(lines) == 0 or lines[0][0] != '>':
      print 'Error: Bad file'
      sys.exit(0)
    seq1 = ''.join(lines[1:])

  with open(sys.argv[2]) as f:
    lines = f.readlines()
    if len(lines) == 0 or lines[0][0] != '>':
      print 'Error: Bad file'
      sys.exit(0)
    seq2 = ''.join(lines[1:])

  optimize(seq1, seq2)

def optimize(seq1, seq2):
  ms = 1
  mms = -2
  go = -2
  ge = -1

  results = defaultdict(list)
  for i in np.arange(-0.6, -4.5, -0.1):
    mms = i
    for j in np.arange(0, -4.5, -0.1):
      go = j
      for k in np.arange(0, -3.5, -0.1):
        ge = k
        # print ms, mms, go, ge
        stats = locAL.external(seq1, seq2, ms, mms, go, ge)
        if stats not in results:
          results[stats] = (ms, mms, go, ge)
        else:
          results[stats].append((ms, mms, go, ge))

  for i in results:
    print float(i[1])/float(i[0]), i, results[i]



# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start