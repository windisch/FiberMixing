#!/usr/bin/python
import numpy as np
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import itertools
import argparse
from random import randint
from subprocess import call
from sage.all import *
from sage.interfaces.four_ti_2 import FourTi2
from multiprocessing import Process
from multiprocessing import Pool
from multiprocessing import freeze_support

C_THREADS=1
C_LATTEBIN='count'
C_LATTEDIR=''
C_LATTEPAR='--count-lattice-points'
C_LATTEOUT='numOfLatticePoints'

def countFiber(A,v):
   b=np.array([[i for i in A.dot(v)]])
   global C_LATTEBIN
   global C_LATTEDIR
   global C_LATTEPAR
   global C_LATTEOUT
   #latteBin="/opt/latte-integrale-1.7.1/bin/count"
   #lattePara="--count-lattice-points"
   #latteOut="numOfLatticePoints"
   tDir=tmp_dir()
   tFile=tDir+'fiberCount.latte'
   with open(tFile,'wb') as f:
      numRows=(A.shape)[0]
      numCols=(A.shape)[1]
      f.write(str(numRows)+' '+str(numCols+1)+'\n')
      np.savetxt(f, np.hstack((b.T,-A)),fmt='%d')
      #write equality constraints
      f.write('linearity '+str(numRows)+' ')
      for i in range(1,numRows+1):
         f.write(str(i)+' ')
      #write nonegative contrains
      f.write('\nnonnegative '+str(numCols)+' ')
      for i in range(1,numCols+1):
         f.write(str(i)+' ')
   os.chdir(tDir)
   #run LattE
   call([C_LATTEDIR+C_LATTEBIN,C_LATTEPAR,tFile])
   #parse output
   with open(tDir+C_LATTEOUT,'r') as f:
      return int(f.read())


def averageMixingTime(A,M,u,m,verbose=False):
   global C_THREADS
   n=countFiber(A,u)
   #x=np.array([uniformWalk(A,M,u,n,verbose=verbose) for i in xrange(m)])

   A_args=itertools.repeat(A,m)
   M_args=itertools.repeat(M,m)
   u_args=itertools.repeat(u,m)
   n_args=itertools.repeat(n,m)
   v_args=itertools.repeat(verbose,m)
   p = Pool(C_THREADS)
   x=np.array(p.map(uniformWalk_par,itertools.izip(A_args,M_args,u_args,n_args,v_args)))

   plt.hist(x,facecolor='c',bins=20)
   plt.xlabel('Number of steps')
   plt.ylabel('Number of occurences')
   plt.title(r'Mixing time')
   plt.axvline(x.mean(), color='b', linestyle='dashed', linewidth=2)
   #plt.axis([40, 160, 0, 0.03])
   plt.grid(True)
   plt.show()
   return x.mean() 

def hypergeoWalk(M,u,n):
   return 'TODO'

#thining; burn-in


def uniformWalk_par(args):
    return uniformWalk(*args)

def uniformWalk(A,M,u,n=0,verbose=False):
   if n==0:
       n=countFiber(A,u)
       print n
   numMoves=len(M)
   T={}
   i=0
   d=1
   while d>0.25:
      i+=1
      #update hash table
      if tuple(u) in T:
         T[tuple(u)]+=1
      else:
         T[tuple(u)]=1

      #check distance to uniform
      d=0.5*sum([abs((float(T[x])/i)-1/float(n)) for x in T])+ 0.5*(float(n-len(T)))/n
      if verbose:
         if i%1000==0:
            print d,i

      #sample move
      m=(2*randint(0,1)-1)*M[randint(0,numMoves-1)]
      if all(i>=0 for i in u+m):
           u=u+m
   return i


def main(argv):
   global C_THREADS
   global C_LATTEDIR

   #Parse arguments
   parser = argparse.ArgumentParser(description='Fiber walks mixing time estimation')
   parser.add_argument('--matrix',
                       dest='matrix',
                       metavar='mat.file',
                       type=str,
                       help='the constraint matrix')
   parser.add_argument('--markov',dest='markov',
                       metavar='mar.file',
                       type=str,
                       help='the Markov basis')
   parser.add_argument('--initial',dest='initial', metavar='ini.file', type=str,
                   help='the initial node')
   parser.add_argument('--runs',dest='runs',metavar='N',type=int,default=1,
                   help='the number of times the random walk is runned')
   parser.add_argument('--threads',dest='threads',metavar='N',type=int,default=1,
                   help='the number of threads used')
   parser.add_argument('--latte',dest='latte',metavar='dir',type=str,default='',
                   help='path to LattE binaries')
   parser.add_argument('--verbose', help="turn on verbose-mode",action="store_true")
   args = parser.parse_args()

   #read input
   f=FourTi2(os.curdir)
   A=f.read_matrix(args.matrix)
   M=f.read_matrix(args.markov)
   u=f.read_matrix(args.initial)
   C_THREADS=args.threads
   C_LATTEDIR=args.latte

   #convert input
   A=np.array(A)
   u=(np.array(u))[0]
   M=[m for m in M]
   print averageMixingTime(A,M,u,m=args.runs,verbose=args.verbose)

if __name__ == "__main__":
    freeze_support()
    main(sys.argv)
