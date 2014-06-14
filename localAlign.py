# BIMM 185, Spring 2014
# Max Shen
# localAlign.py
#
# localAlign.py:
#   python localAlign <input file> <scoringfile>
#     <input file>: A fasta filename. The 1st two sequences in the file are used.
#     <scoringfile>: A scoring matrix for protein alignment.
#   Summary: localAlign.py performs protein sequence alignment. It is not used
#     in the final project, but represents preliminary work on extending the
#     project to proteins.


import sys
import string
import datetime
import random

from collections import defaultdict

def main():
  if len(sys.argv) != 3:
    print 'Usage: python localAlign <input file> <scoringfile>'
    sys.exit(0)

  with open(sys.argv[1]) as f:
    lines = f.read().splitlines()

  q1 = lines[1]
  q2 = lines[3]

  localAlign(q1, q2, sys.argv[2])

def localAlign(q1, q2, scoringfile):
  aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
   'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'Z', 'X', '*']
  scoringmat = readMatrix(scoringfile)
  print 'Aligning with', scoringfile
  indel = float(scoringmat[-1][1])
  q1 = ' ' + q1
  q2 = ' ' + q2
  s = [[0]*(len(q1)) for i in range(len(q2))]
  bt = [['']*len(q1) for i in range(len(q2))]

  # print scoringmat

  # Base case
  for i in range(1, len(q2)):
    s[i][0] = 0
    bt[i][0] = 'none'
  for j in range(1, len(q1)):
    s[0][j] = 0
    bt[0][j] = 'none'
  
  best = {'i': 0, 'j': 0, 'score': 0}

  for i in range(1, len(q2)):
    # print i
    for j in range(1, len(q1)):
      nones = 0
      ups = s[i-1][j] + indel
      lefts = s[i][j-1] + indel
      q1i = aa.index(q1[j])
      q2i = aa.index(q2[i])
      diags = s[i-1][j-1] + int(scoringmat[q1i][q2i])
      s[i][j] = max(ups, lefts, diags, nones)
      if s[i][j] > best['score']:
        best['i'] = i
        best['j'] = j
        best['score'] = s[i][j]

      # Update bt
      if s[i][j] == 0:
        bt[i][j] = 'none'
      elif s[i][j] == ups:
        bt[i][j] = 'up'
      elif s[i][j] == lefts:
        bt[i][j] = 'left'
      elif s[i][j] == diags:
        bt[i][j] = 'diag'

  # print s
  # for i in bt:
  #   print i
  f1 = []
  f2 = []
  printBT(s, bt, q1, q2, best['i'], best['j'], scoringmat, aa)
  print 'Score =', s[best['i']][best['j']]

def readMatrix(scoringfile):
  mat = []
  save = False
  with open(scoringfile) as f:
    for i, line in enumerate(f):
      if save:
        mat.append(line.split()[1:])
      if line[0] != '#':
        save = True
  return mat

def printBT(s, bt, q1, q2, i, j, scoringmat, aa):
  f1 = []
  f2 = []
  while True:
    # print s[i][j]
    if bt[i][j] == 'none':
      break
    elif bt[i][j] == 'up':
      f1.insert(0, '-')
      f2.insert(0, q2[i])
      i -= 1
    elif bt[i][j] == 'left':
      f1.insert(0, q1[j])
      f2.insert(0, '-')
      j -= 1
    else:
      f1.insert(0, q1[j])
      f2.insert(0, q2[i])
      i -= 1
      j -= 1

  seq1 = ''
  mid = ''
  seq2 = ''
  matches = float(0)
  posmatches = float(0)
  mismatches = 0
  gapopen = 0
  gapextend = 0
  totalLen = float(len(f1))
  for i in range(len(f1)):
    seq1 += f1[i]
    seq2 += f2[i]
    if f1[i] == f2[i]:
      mid += '|'
      matches += 1
    elif f1[i] == '-' or f2[i] == '-':
      mid += ' '
      if i > 0:
        if f1[i-1] == '-' or f2[i-1] == '-':
          gapextend += 1
        else:
          gapopen += 1
    else:
      f1i = aa.index(f1[i])
      f2i = aa.index(f2[i])
      if int(scoringmat[f1i][f2i]) > 0:
        mid += '+'
        posmatches += 1
      else:
        mid += '*'
        mismatches += 1

  print seq1, '\n', mid, '\n', seq2
  print 'Accuracy =', 100*(matches+posmatches) / totalLen, '%'
  print 'Length =', len(f1)
  print 'M:', int(matches), '\tPM:', int(posmatches), '\tMM:', mismatches, '\tGO:', gapopen, '\tGE:', gapextend
      

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\nEnd:', end, '\nTotal:', end - start