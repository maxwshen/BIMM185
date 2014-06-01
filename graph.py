# Run on Windows, not debruijn
import numpy as np
import matplotlib.pyplot as plt
import sys
import datetime
import string
import random


def main():
  if len(sys.argv) != 2:
    print 'Usage: python graph.py <input_file>'
    sys.exit(0)
  graph(sys.argv[1])

def graph(inputfile):
  with open(inputfile) as f:
    lines = f.readlines()

  plt.plot(lines)

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start
  main()
  end = datetime.datetime.now()
  print '\nEnd:', end, '\nTotal:', end - start