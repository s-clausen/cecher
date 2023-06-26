# Cecher

Copyright © 2023 Sönke Clausen


### Description

The Cecher is a C++ software for efficient computation of Čech persistence barcodes. It outperforms current software ([Dionysus], [GUDHI]) and is especially adapted to handle high ambient dimensions. Its main features are:

  - a core persistence algorithm adapted from the C++ code of the ultrafast [Ripser], featuring its support for coefficients in prime fields, its memory-efficient philosophy, and use of a combination of cohomology, clearing, and zero-apparent pairs shortcut
  - a version of the union-find algorithm in [Ripser] with compression
  - a symbolic perturbation of the Čech filtration, allowing to completely skip the construction of many columns of the coboundary matrix corresponding to zero-persistence pairs
  - a novel minimal enclosing sphere implementation, outperforming established software like [Miniball] in the context of a persistence algorithm and featuring efficient computation in high ambient dimensions
  - a lazy-exact type computation, where computations are performed primarily with an interval arithmetic type
    

### Assumption for Efficient Computation

The Cecher is efficient on noisy 'real world data' (more specifically on datasets in *d*-dimensional Euclidean space, where all simplices of dimension *d* or lower are affinely independent and have distinct minimal circumspheres with distinct radii). The implementation is robust to violation of this assumption.


### Input Format

Currently supported is an input in the form of a rational point cloud, where numerator and denominator are seperated by a forward-slash, and values are seperated by a whitespace ([examples]).


### Dependencies 

The implementation currently utilizes the (double-precision) interval arithmetic type and exact rational type of [CGAL].


### Options

Cecher supports the following preprocessor macros:

  - `USE_COEFFICIENTS`: enable support for coefficients in a prime field
  - `EXPLICIT_CENTERS`: minimal enclosing sphere algorithm explicitly computes the center (can increase efficiency for low ambient dimensions)



The following options can be passed in the file name or as program argument:

  - `--dim k`: compute persistent homology up to dimension *k* (default: *k*=1).
  - `--modulus p`: compute homology with coefficients in the finite prime field Z/*p*Z (when built with `USE_COEFFICIENTS`, default: *p*=2).
  - `--threshold t`: compute Čech filtration up to radius *t* (default: no threshold ).





[Dionysus]: <http://www.mrzv.org/software/dionysus/>
[GUDHI]: <http://gudhi.gforge.inria.fr>
[Ripser]: <https://github.com/Ripser/ripser>
[CGAL]: <https://github.com/CGAL/cgal>
[Miniball]: <https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html>
[examples]: <https://github.com/s-clausen/cecher/tree/main/examples>


