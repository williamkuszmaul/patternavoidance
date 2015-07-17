# To run, type in terminal:
# sage
# attach("sagetest.sage")
# test()

import sys
import time
from sage.all import *

def test():
  t0 = time.time() # Note this is not very precise timing
  print len(Permutations(8, avoiding=[[3,4,1,2], [4,2,3,1]]));
  t1 = time.time()
  print "Time elapsed (s): ", (t1 - t0);