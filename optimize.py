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

  global ms
  global mms
  global go
  global ge
  global mms_lowbound
  global mms_highbound
  global go_highbound
  global ge_highbound
  ms = 1
  mms = -2
  go = -2
  ge = -1
  mms_lowbound = -4.5
  mms_highbound = -0.6
  go_highbound = 0
  ge_highbound = 0


  optimize(seq1, seq2)

def optimize(seq1, seq2):
  seed = [-1, -1, -1]
  iterate = True
  while iterate:
    explore(seq1, seq2, seed, .1)

    iterate = False

  return

  # Brute force 
  results = defaultdict(list)
  for i in np.arange(-0.6, -4.5, -0.2):
    mms = i
    for j in np.arange(0, -4.5, -0.2):
      go = j
      for k in np.arange(0, -3, -0.2):
        ge = k
        # print ms, mms, go, ge
        stats = locAL.external(seq1, seq2, ms, mms, go, ge)
        if stats not in results:
          results[stats] = [(ms, mms, go, ge)]
        else:
          results[stats].append((ms, mms, go, ge))

  for i in results:
    print float(i[1])/float(i[0]), i, ':',
    for j in results[i]:
      print '({0}, {1}, {2}, {3})'.format(round(j[0], 2), round(j[1], 2), round(j[2], 2), round(j[3], 2)),
    print '\n',

# Explores seed in all directions in 3D Space (27 possibilities) with 
def explore(seq1, seq2, seed, step):
  for i in [-step, 0, step]:
    for j in [-step, 0, step]:
      for k in [-step, 0, step]:
        mms_t = seed[0] + i
        go_t = seed[1] + j
        ge_t = seed[2] + k
        if mms_lowbound < mms_t < mms_highbound and go_t < go_highbound and ge_t < ge_highbound:
          print mms_t, go_t, ge_t
          stats = locAL.external(seq1, seq2, ms, mms_t, go_t, ge_t)
          print stats, float(stats[1]*100)/float(stats[0])


# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
