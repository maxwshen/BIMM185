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
  global go_lowbound
  global ge_highbound
  global ge_lowbound
  ms = 1
  mms = -2
  go = -2
  ge = -1
  mms_lowbound = -4.5
  mms_highbound = -0.6
  go_highbound = 0
  go_lowbound = -4.5
  ge_highbound = 0
  ge_lowbound = -3


  optimize(seq1, seq2)

def optimize(seq1, seq2):
  seeds = []
  traversed = set()
  gridsize = 0.5
  numiterations = 0
  for i in np.arange(mms_lowbound + gridsize/2, mms_highbound - gridsize/2, gridsize):
    for j in np.arange(go_lowbound + gridsize/2, go_highbound - gridsize/2, gridsize):
      for k in np.arange(ge_lowbound + gridsize/2, ge_highbound - gridsize/2, gridsize):
        seeds.append([i, j, k])
  print len(seeds), 'starting seeds at distance', gridsize
  seeds = list(np.random.permutation(seeds))

  best = defaultdict(list)
  best[0] = []
  while len(seeds) > 0:
    # print numiterations
    seed = seeds[0]
    seeds = seeds[1:]
    if len(seeds) == 0:
      for key in best:
        for val in best[key]:
          seeds.append(val)

    if tuple(seed) not in traversed:
      traversed.add(tuple(seed))

      # Store any positions whose accuracy is greater than the smallest
      # in the dict best. Best is initialized containing 0 which is
      # removed after more than 20 items are inserted into the list.
      tries = explore(seq1, seq2, seed, .1)
      for key in tries:
        if key > min(best):
          if key not in best:
            best[key].append(list(np.random.permutation(tries[key])[0]))
          else:
            best[key].append(list(np.random.permutation(tries[key])[0]))

      # Remove all keys outside the top 20. Keep only top 20
      if len(best) > 20:
        keys = best.keys()
        keys.sort(reverse = True)
        tempdict = defaultdict(list)
        for i in range(20):
          tempdict[keys[i]] = best[keys[i]]
        best = tempdict
        # print best, numiterations
      # print best

    print np.average(best.keys())
    numiterations += 1
    if numiterations > 400:
      break


  return

# Explores seed in all directions in 3D Space (27 possibilities) with 
# _step size
# Output:
#   Returns a random seed of the highest accuracy 
def explore(seq1, seq2, seed, step):
  best = defaultdict(list)

  for i in [-step, 0, step]:
    for j in [-step, 0, step]:
      for k in [-step, 0, step]:
        mms_t = seed[0] + i
        go_t = seed[1] + j
        ge_t = seed[2] + k
        if mms_lowbound < mms_t < mms_highbound and go_t < go_highbound and ge_t < ge_highbound:
          # print mms_t, go_t, ge_t
          stats = locAL.external(seq1, seq2, ms, mms_t, go_t, ge_t)
          accuracy = float(stats[1]*100)/float(stats[0])
          # print stats, accuracy
          best[accuracy].append((mms_t, go_t, ge_t))

  return best

  # print max(best), best[max(best)]
  # print best
  if max(best) > 60:
    direction = np.random.permutation(best[max(best)])[0]
    bestdir = [seed[0] + direction[0], seed[1] + direction[1], seed[2] + direction[2]]
    return bestdir, max(best)
  else:
    return None, max(best)


# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
