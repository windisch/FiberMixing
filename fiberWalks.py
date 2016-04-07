#!/usr/bin/python
import numpy as np
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import itertools
import argparse
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
C_CURDIR=os.path.dirname(os.path.abspath(__file__))+'/'
C_OUTFILE_FIG='out.eps'
C_OUTFILE_DAT='out.txt'

#walk specific parameters
DEF_THINNING=1
DEF_VERBOSE=False
DEF_RUNS=1
DEF_BURNIN=0
DEF_FIBERSIZE=-1

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

def writeOutfile(A,M,u,fibersize,runs,thinning,burnin,mt):
   global C_CURDIR
   global C_OUTFILE_DAT
   #save array to outfile
   with open(C_CURDIR+C_OUTFILE_DAT,'wb') as f:
      f.write('################################\n')
      f.write('Average mixing time:\t\t\t'+str(mt.mean())+'\n')
      f.write('Number of walks:\t\t\t\t'+str(runs)+'\n')
      f.write('Number of lattice points:\t'+str(fibersize)+'\n')
      f.write('Thinning of walk:\t'+str(thinning)+'\n')
      f.write('Burn-in:\t'+str(burnin)+'\n')
      f.write('################################\n')
      f.write('##Constraint Matrix\n')
      np.savetxt(f,A,fmt='%d')
      f.write('##Markov moves\n')
      np.savetxt(f,M,fmt='%d')
      f.write('##Initial vertex\n')
      np.savetxt(f,u,fmt='%d')
      f.write('##Observed mixing times\n')
      np.savetxt(f,mt.T,fmt='%d')

def estimateMixing(A,M,u,fibersize,runs,verbose,thinning,burnin):
   global C_THREADS
   global C_CURDIR

   A_args=itertools.repeat(A,runs)
   M_args=itertools.repeat(M,runs)
   u_args=itertools.repeat(u,runs)
   fibersize_args=itertools.repeat(fibersize,runs)
   verbose_args=itertools.repeat(verbose,runs)
   thinning_args=itertools.repeat(thinning,runs)
   burnin_args=itertools.repeat(burnin,runs)
   p = Pool(C_THREADS)
   mt=np.array(p.map(walk_par,itertools.izip(A_args,M_args,u_args,fibersize_args,verbose_args,thinning_args,burnin_args)))
   print mt
   mean=mt.mean()

   plt.hist(mt,facecolor='c',bins=20)
   plt.xlabel('Number of steps')
   plt.ylabel('Number of occurences')
   plt.title(r'Mixing time')
   plt.axvline(mean, color='b', linestyle='dashed', linewidth=2)
   #plt.axis([40, 160, 0, 0.03])
   plt.grid(True)
   #save figure to outfile
   plt.savefig(C_CURDIR+C_OUTFILE_FIG,format='eps',dpi=1000)
   #write output file 
   writeOutfile(A,M,u,fibersize,runs,thinning,burnin,mt)
   #plt.show()
   return mt.mean()

def next(u,M):
   m=(2*np.random.randint(0,2)-1)*M[np.random.randint(0,len(M))]
   if all(j>=0 for j in u+m):
      u=u+m
   return u

def walk_par(args):
   #seed random number generator
   np.random.seed()
   return walk(*args)

def walk(A,M,u,n,verbose,thinning,burnin):
   seed()
   T={}
   i=0
   k=0
   d=1
   #burnin
   for i in xrange(burnin):
      u=next(u,M)
   while d>0.25:
      k+=1
      if k%thinning==0:
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

      #sample next fiber element
      u=next(u,M)
   return k

def main(argv):
   global C_CURDIR
   global C_THREADS
   global C_LATTEDIR
   global DEF_THINNING
   global DEF_BURNIN
   global DEF_VERBOSE
   global DEF_RUNS

   #Parse arguments
   parser = argparse.ArgumentParser(description='Fiber walks mixing time estimation')

   #required arguments
   requiredArgs = parser.add_argument_group('required arguments')
   requiredArgs.add_argument('-c','--matrix',
                       dest='matrix',
                       metavar='mat.file',
                       type=str,
                       help='path to constraint matrix')
   requiredArgs.add_argument('-m','--markov',dest='markov',
                       metavar='mar.file',
                       type=str,
                       help='path to Markov basis')
   requiredArgs.add_argument('-i','--initial',
                       dest='initial',
                       metavar='ini.file',
                       type=str,
                       help='path to initial node')
   #optional arguments
   parser.add_argument('-l','--latte',dest='latte',metavar='/latte/path/',type=str,default=C_LATTEDIR,
                   help='path to LattE binaries')
   parser.add_argument('-r','--runs',dest='runs',metavar='N',type=int,default=DEF_RUNS,
                   help='number of random walk runs, default is '+str(DEF_RUNS))
   parser.add_argument('-f','--fiber-size',dest='fibersize',metavar='N',type=int,default=DEF_FIBERSIZE,
                   help='the size of the fiber')
   parser.add_argument('-t','--threads',dest='threads',metavar='N',type=int,default=C_THREADS,
                   help='number of threads used, default is '+str(C_THREADS))
   parser.add_argument('-s','--thinning',dest='thinning',metavar='N',type=int,default=DEF_THINNING,
                   help='take only every Nth sample, default is '+str(DEF_THINNING))
   parser.add_argument('-b','--burn-in',dest='burnin',metavar='N',type=int,default=DEF_THINNING,
                   help='do a burn-in with N steps before sampling, default is '+str(DEF_BURNIN))
   parser.add_argument('-v','--verbose', help="turn on verbose-mode",action="store_true",default=DEF_VERBOSE)
   args = parser.parse_args()

   #read input
   f=FourTi2(os.curdir)
   A=f.read_matrix(args.matrix)
   M=f.read_matrix(args.markov)
   u=f.read_matrix(args.initial)
   C_THREADS=args.threads
   C_LATTEDIR=args.latte
   thinning=args.thinning
   burnin=args.burnin
   #convert input
   A=np.array(A)
   u=(np.array(u))[0]
   M=[m for m in M]
   if args.fibersize < 0:
      args.fibersize=countFiber(A,u)
   print estimateMixing(A,M,u,args.fibersize,args.runs,args.verbose,thinning,burnin)

if __name__ == "__main__":
    freeze_support()
    main(sys.argv)
