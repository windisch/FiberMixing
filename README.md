# FiberMixing

### Usage
The following performs the simple walk on the fiber of the
matrix specified in `matrix.mat` which starts at `initial.mat` and
which uses the Markov moves as specified in `markov.mat`

```bash
sage -python fiberWalks.py --matrix examples/matrix.mat 
                           --markov examples/markov.mat 
                           --initial examples/initial.mat 
```

*All files containing matrices or vectors have to be in the
[4ti2](www.4ti2.de)/[LattE](https://www.math.ucdavis.edu/~latte/)-format, for example*

```R
2 4
1 2 3 4
5 6 7 8
9 0 1 2
3 4 5 6

```


The random walk runs until the approximated distribution differs from
the uniform distribution by at most 0.25 (i.e. the mixing time). The
number of steps needed is then saved in `out.txt`. In case that the
random walk is runned mulitple times (with the parameter `-r,--runs`),
a histogram is saved in `out.eps`.

### Requirements
The software [LattE](https://www.math.ucdavis.edu/~latte/) is used to
count the number of integer points in the fiber. Unless the binaries
of Latte are not in a path from `$PATH`, the path to LattE can be
specified with `-l,--latte`.


### Multiple runs and parallel computing
With the parameter `-r,--runs`, the number of runs of the random walk
can be specified, which will be processed consecutively. These random
walks can be perfomed in parallel by specifying the number of threads
with `-t,--threads`.

```bash
sage -python fiberWalks.py --matrix examples/matrix.mat 
                           --markov examples/markov.mat  
                           --initial examples/initial.mat 
                           --runs 1000 
                           --threads 5
```

##Walk-tuning
The parameters `-s,--thinning` and `-b,--burn-in` allow to thin out
the walk or to do a burn-in first. The size of the fiber can be
specified directly with `-f,--fiber-size` and in this case the LattE
computation is ommited.

```bash
sage -python fiberWalks.py --matrix examples/matrix.mat 
                           --markov examples/markov.mat  
                           --initial examples/initial.mat 
                           --runs 1000 
                           --threads 5
                           --thinning 10
                           --burn-in 50
                           --fiber-size 55
```
