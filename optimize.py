import sys
import string
import datetime
import random
import numpy as np
import locAL

from collections import defaultdict

def main():
  if len(sys.argv) != 4:
    print 'Usage: python optimize <input file> <input file 2> <cullingPct>'
    sys.exit(0)

  # Read in seqs
  with open(sys.argv[1]) as f:
    lines = f.readlines()
    if len(lines) == 0 or lines[0][0] != '>':
      print 'Error: Bad fasta file'
      sys.exit(0)
    seq1name = lines[0].strip()
    seq1 = ''.join(lines[1:])

  with open(sys.argv[2]) as f:
    lines = f.readlines()
    if len(lines) == 0 or lines[0][0] != '>':
      print 'Error: Bad fasta file'
      sys.exit(0)
    seq2name = lines[0].strip()
    seq2 = ''.join(lines[1:])

  if float(sys.argv[3]) <= 0 or float(sys.argv[3]) > 1:
    print 'Bad cullingPct'
    sys.exit(0) 

  print 'Training on:\n\t', seq1name, '\n\t', seq2name, '\n'

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
  global gridsize
  global cullingPct
  global cutMM
  global cutGaps
  global threshold
  global cutPercent
  ms = 1
  mms = -2
  go = -2
  ge = -1
  mms_lowbound = -4.0
  mms_highbound = -0.6
  go_highbound = 0
  go_lowbound = -4.0
  ge_highbound = 0
  ge_lowbound = -3
  # gridsizeAll = [0.5, 0.3, 0.2, 0.1]
  # gridsizeAll = [0.5, 0.25]
  gridsizeAll = [0.5]
  gridsize = 0.5
  cullingPct = float(sys.argv[3])
  cutMM = False
  cutGaps = False
  threshold = 1.50
  cutPercent = 0.50

  for i in range(len(gridsizeAll)):
    gridsize = gridsizeAll[i]
    print 'Training Iteration', i + 1
    optimize(seq1, seq2)

# Performs search space reduction using seq1 and seq2 as training
# Input:
#   seq1: A DNA string
#   seq2: A DNA string
# Output:
#   Modifies mms_lowbound, mms_highbound, go_lowbound, go_highbound,
#     ge_lowbound, and ge_highbound based on the most accurate
#     alignments found in the search space.
def optimize(seq1, seq2):
  global cutMM
  global cutGaps
  global threshold
  global cutPercent

  seeds = []
  traversed = set()
  numiterations = 0

  if cutMM is True:
    mms_width = mms_lowbound - mms_highbound
    mms_cutHigh = mms_highbound + (mms_width * (cutPercent / 2))
    mms_cutLow = mms_highbound + (mms_width * (1 - cutPercent / 2))
    print 'New mms:', mms_cutLow, mms_cutHigh
    for i in np.arange(mms_cutLow, mms_cutHigh, gridsize):
      for j in np.arange(go_lowbound, go_highbound, gridsize):
        for k in np.arange(ge_lowbound, ge_highbound, gridsize):    
          seeds.append([i, j, k])
  elif cutGaps is True:
    go_width = go_lowbound - go_highbound
    go_cutHigh = go_highbound + (go_width * (cutPercent / 2))
    go_cutLow = go_highbound + (go_width * (1 - cutPercent / 2)) 
    
    ge_width = ge_lowbound - ge_highbound
    ge_cutHigh = ge_highbound + (ge_width * (cutPercent / 2))
    ge_cutLow = ge_highbound + (ge_width * (1 - cutPercent / 2))
    print 'New go:', go_cutLow, go_cutHigh
    print 'New ge:', ge_cutLow, ge_cutHigh  
    for i in np.arange(mms_lowbound, mms_highbound, gridsize):
      for j in np.arange(go_cutLow, go_cutHigh, gridsize):
        for k in np.arange(ge_cutLow, ge_cutHigh, gridsize):   
          seeds.append([i, j, k])    
  else:
    for i in np.arange(mms_lowbound, mms_highbound, gridsize):
      for j in np.arange(go_lowbound, go_highbound, gridsize):
        for k in np.arange(ge_lowbound, ge_highbound, gridsize):
          seeds.append([i, j, k])
  print len(seeds), 'starting seeds at step size', gridsize
  print 'Current Search Space:\n\tmms:', mms_lowbound, mms_highbound, '\tgo:', go_lowbound, go_highbound, '\tge:', ge_lowbound, ge_highbound
  seeds = list(np.random.permutation(seeds))
  totalseeds = len(seeds)

  totalGaps = 0
  totalMM = 0

  best = defaultdict(list)
  smallest = 0
  overflow = False
  # best[0] = []
  while len(seeds) > 0:
    seed = seeds[0]
    seeds = seeds[1:]
    if len(seeds) == 0:
      break

    if tuple(seed) not in traversed:
      traversed.add(tuple(seed))
      if overflow:
        smallest = min(best.keys())

      # Store any positions whose accuracy is greater than the smallest
      # in the dict best. Best is initialized containing 0 which is
      # removed after more than 15 items are inserted into the list.
      # Do not explore
      stats = locAL.external(seq1, seq2, ms, seed[0], seed[1], seed[2])
      accuracy = float(stats[1]*100)/float(stats[0])
      totalMM += stats[2]
      totalGaps += stats[3]

      if accuracy > smallest:
        best[accuracy].append(seed)

      # Keep only the best keys
      if len(best) > 15:
        overflow = True
        keys = best.keys()
        keys.sort(reverse = True)
        tempdict = defaultdict(list)
        for i in range(15):
          tempdict[keys[i]] = best[keys[i]]
        best = tempdict
        # print best, numiterations
      # print best


    # print numiterations, len(seeds), np.average(best.keys())
    # print numiterations, np.average(best.keys())
    sys.stdout.write('\r' + str(numiterations) + ' / ' + str(totalseeds) + ' '*10)
    sys.stdout.flush()
    numiterations += 1

  print 'Done.\nAvg top accuracy:', np.average(best.keys())
  print 'totalGaps:', totalGaps, 'totalMM:', totalMM
  print best.keys()

  setRange(best)
  cutMM = False
  cutGaps = False
  if totalMM > threshold * totalGaps:
    cutGaps = True
    print 'Alignment mutations are primarily mismatches. Next iteration, gap search space will be reduced to', 100 * cutPercent, '%. Threshold:', threshold 
  elif totalGaps > threshold * totalMM:
    cutMM = True
    print 'Alignment mutations are primarily gaps. Next iteration, mismatch search space will be reduced to', 100 * cutPercent, '%. Threshold:', threshold

  return

# Input: best, a dictionary.
#   Keys are accuracy percents
#   Values are coordinates
# Output:
#   Lowest and highest values in best for mms, go, and ge
def setRange(best):
  global mms_lowbound
  global mms_highbound
  global go_lowbound
  global go_highbound
  global ge_lowbound
  global ge_highbound
  global cullingPct

  _mms = []
  _go = []
  _ge = []
  _mmsCulled = []
  _goCulled = []
  _geCulled = []

  for key in best:
    for value in best[key]:
      _mms.append(value[0])
      _go.append(value[1])
      _ge.append(value[2])

  cullingNum = int(round(cullingPct * len(_mms)))

  arr = np.arange(len(best.keys()))

  np.random.shuffle(arr)
  for i in range(cullingNum):
    _mmsCulled.append(_mms[i])

  np.random.shuffle(arr)
  for i in range(cullingNum):
    _goCulled.append(_go[i])

  np.random.shuffle(arr)
  for i in range(cullingNum):
    _geCulled.append(_ge[i])

  initVolume = (mms_highbound - mms_lowbound) * (go_highbound - go_lowbound) * (ge_highbound - ge_lowbound)

  # print 'mms:', min(_mms), max(_mms)
  # print 'go:', min(_go), max(_go)
  # print 'ge:', min(_ge), max(_ge)
  naiveVolume = (max(_mms) - min(_mms)) * (max(_go) - min(_go)) * (max(_ge) - min(_ge))

  # print 'mms:', min(_mmsCulled), max(_mmsCulled)
  # print 'go:', min(_goCulled), max(_goCulled)
  # print 'ge:', min(_geCulled), max(_geCulled)
  culledVolume = (max(_mmsCulled) - min(_mmsCulled)) * (max(_goCulled) - min(_goCulled)) * (max(_geCulled) - min(_geCulled))

  print 'No Culling : Search space reduced to', 100 * naiveVolume / initVolume, '%'
  print 'Culling at', cullingPct, ': Search space reduced to', 100 * culledVolume / initVolume, '%'
  print 'Culling Improvement:', (naiveVolume - culledVolume) / initVolume, '%' 

  mms_lowbound = min(_mmsCulled)
  mms_highbound = max(_mmsCulled)
  go_lowbound = min(_goCulled)
  go_highbound = max(_goCulled)
  ge_lowbound = min(_geCulled)
  ge_highbound = max(_geCulled)

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
  print 'Efficient Exploration of Rational Sequence Alignment (EERSA) v1.0'
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start
