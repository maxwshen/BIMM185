# locAL.py
# BIMM 185
# Max Shen A10082759
#
# python locAL.py <seq file1> <seq file2> -m <match> -s <mismatch>
# -go <gap-open> -ge gap-extend
#
# locAL gives the best local alignment of two seq files based on an affine
# cost function provided by the user, using the Smith-Waterman Algorithm.
# It outputs the best score, the resulting sequence, and its length.
#
# locAL is not space efficient.

import sys
import copy
import string
import datetime
import random

from collections import defaultdict

# Backtracking variable names
bt_zero = 0
bt_up = 1
bt_left = 2
bt_diag = 3

def main():
  if len(sys.argv) < 3 or len(sys.argv) % 2 == 0:
    print 'Usage: python locAL.py <seq file1> <seq file2> -m <match> -s <mismatch> -go <gap-open> -ge gap-extend'
    sys.exit(0)

  # Default values
  global ms
  global mms
  global go
  global ge
  global silenced
  ms = 1
  mms = -2
  go = -2
  ge = -1
  silenced = False

  # Detect flags and change parameters
  for param in sys.argv[3:]:
    if param == '-m':
      ms = float(sys.argv[sys.argv.index('-m') + 1])
    if param == '-s':
      mms = float(sys.argv[sys.argv.index('-s') + 1])
    if param == '-go':
      go = float(sys.argv[sys.argv.index('-go') + 1])
    if param == '-ge':
      ge = float(sys.argv[sys.argv.index('-ge') + 1])

  print 'Match = ', ms, '\tMismatch = ', mms
  print 'Gap Open = ', go, '\t\tGap Extend = ', ge

  # Ensure scoring values are valid
  if ms < 0 or mms > 0 or go > 0 or ge > 0:
    print 'ERROR: Bad scoring parameter values'
    sys.exit(0)
  if ms + 3 * mms > 0:
    print 'ERROR: Positive expected score for background'
    print ms - 3*mms
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

  return locAL(seq1, seq2)

# Used when calling locAL from another python script
def external(seq1, seq2, mscore, mmscore, goscore, gescore):
  global ms
  global mms
  global go
  global ge
  global silenced
  ms = float(mscore)
  mms = float(mmscore)
  go = float(goscore)
  ge = float(gescore)
  silenced = True

  return locAL(seq1, seq2)

# Returns the numeric alignment score between the input bases
def score(base1, base2):
  global ms, mms, go, ge
  if base1 == 'go' or base2 == 'go':
    return go
  if base1 == 'ge' or base2 == 'ge':
    return ge
  if base1 == base2:
    return ms
  else:
    return mms

# Given the backtracking matrix, prints out the best local alignment
def printAlignment(seq1, seq2, btMatrix, best_xy):
  if not silenced:
    print 'Printing alignment...'
  (i, j) = best_xy
  query = ''
  db = ''
  numgaps = 0
  mismatches = 0
  matches = 0

  while btMatrix[i][j] != bt_zero:
    if btMatrix[i][j] == bt_diag:
      query = string.join((query, seq1[i - 1]), '')
      db = string.join((db, seq2[j - 1]), '')
      if seq1[i - 1] == seq2[j - 1]:
        matches += 1
      else:
        mismatches += 1
      i, j = i - 1, j - 1

    if btMatrix[i][j] == bt_up:
      query = string.join((query, seq1[i - 1]), '')
      db = string.join((db, '-'), '')
      numgaps += 1
      i, j = i - 1, j

    if btMatrix[i][j] == bt_left:
      query = string.join((query, '-'), '')
      db = string.join((db, seq2[j - 1]), '')
      numgaps += 1
      i, j = i, j - 1      

  query = query[::-1]
  db = db[::-1]

  numGapExtends = printSeqs(query, db)
  if not silenced:
    print 'Alignment Length:', len(query), '\nMatches:', matches, '\tMismatches:', mismatches
    print 'Total Gaps:', numgaps, '\tGap Opens:', numgaps - numGapExtends, '\tGap Extends:', numGapExtends

  return (len(query), matches, mismatches, numgaps, numGapExtends)


def locAL(seq1, seq2):
  # s is 2d array of size (seq1 + 1 by seq2 + 1), init. to 0's
  s_mat = [[0]*(len(seq2) + 1) for i in range(len(seq1) + 1)]
  d_mat = [[0]*(len(seq2) + 1) for i in range(len(seq1) + 1)]
  bt_mat = [[0]*(len(seq2) + 1) for i in range(len(seq1) + 1)]

  # initialize d and i, which can have negative scores
  d_mat[1][0] = go + ge
  d_mat[0][1] = go + ge
  for i in xrange(2, len(seq1) + 1):
    d_mat[i][0] = d_mat[i - 1][0] + ge
  for j in xrange(2, len(seq2) + 1):
    d_mat[0][j] = d_mat[0][j - 1] + ge
  i_mat = copy.deepcopy(d_mat)

  best = (0, 0, 0)

  if not silenced:
    print 'Generating alignment... Sequence len:', len(seq1)
  # Generate local alignment
  for i in xrange(1, len(seq1) + 1):
    # print i # track progress
    for j in xrange(1, len(seq2) + 1):
      diag = s_mat[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1])

      d_up = d_mat[i - 1][j] + ge
      s_up = s_mat[i - 1][j] + go + ge
      d_mat[i][j] = max(d_up, s_up)

      i_left = i_mat[i][j - 1] + ge
      s_left = s_mat[i][j - 1] + go + ge
      i_mat[i][j] = max(i_left, s_left)

      s_mat[i][j] = max(0, d_mat[i][j], i_mat[i][j], diag)

      # update bt matrix
      if max(0, d_mat[i][j], i_mat[i][j], diag) == 0:
        bt_mat[i][j] = bt_zero
      if max(0, d_mat[i][j], i_mat[i][j], diag) == d_mat[i][j]:
        bt_mat[i][j] = bt_up
      if max(0, d_mat[i][j], i_mat[i][j], diag) == i_mat[i][j]:
        bt_mat[i][j] = bt_left
      if max(0, d_mat[i][j], i_mat[i][j], diag) == diag:
        bt_mat[i][j] = bt_diag

      if s_mat[i][j] >= best[0]:
        best = (s_mat[i][j], i, j)

  (alignLen, matches, mismatches, numgaps, numGapExtends) = printAlignment(seq1, seq2, bt_mat, best[1:])
  if not silenced: 
    print 'Score =', best[0]
  # print s_mat, '\n'
  # print bt_mat, '\n'
  # print d_mat, '\n'
  # print i_mat

  return (alignLen, matches, mismatches, numgaps, numGapExtends)

# Returns the number of gap extensions
def printSeqs(query, db):
  m = ''
  numGapExtends = 0
  for i in range(len(query)):
    if query[i] == db[i]:
      m += '|'
    elif query[i] == '-' or db[i] == '-':
      m += ' '
      if query[i-1] == '-' or db[i-1] == '-':
        numGapExtends += 1
    else:
      m += '*'
  # if not silenced:
  print query, '\n', m, '\n', db
  return numGapExtends

# main
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start
  main()
  end = datetime.datetime.now()
  print '\nEnd:', end, '\nTotal:', end - start
