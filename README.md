# Cecher

Copyright © 2023 Sönke Clausen


### Description

The Cecher is a C++ software developed for efficient computation of Čech persistence barcodes. It heavily outperforms current software ([Dionysus], [GUDHI]) and is especially adapted to handle high ambient dimensions. Its main features are:

  - a core persistence algorithm adapted from the C++ code of the ultrafast [Ripser], featuring its memory-efficient philosophy, its support for coefficients in prime fields and use of a combination of cohomology, clearing, and zero-apparent pairs shortcut
  - a version of the union-find algorithm in [Ripser] featuring compression
  - a symbolic perturbation of the Čech filtration, allowing to skip the construction of many columns of the coboundary matrix (and with them the corresponding zero-persistence pairs)
  - a novel minimal enclosing sphere implementation, featuring efficient computation in high ambient dimensions and outperforming in the context of a persistence algorithm established software like [Miniball].
  - and a lazy-exact computation, primarily computing with the interval arithmetic, and recomputing with the rational exact type of [CGAL] when necessary 


### General Position Assumption for Efficient Computation

The Cecher is efficient on noisy 'real world data' (more specifically on datasets in *d*-dimensional Euclidean space, where all simplices of dimension *d* or lower are affinely independent and have distinct minimal circumspheres with distinct radii). The implementation is robust to violation of this assumption.


### Input Format

Currently supported is an input in the form of a rational point cloud, where numerator and denominator are seperated by a forward-slash, and values are seperated by a comma ([examples]).


### Dependencies 

The implementation currently utilizes the (double-precision) interval arithmetic type and exact rational type of [CGAL].


### Options

Cecher supports the following preprocessor macros:

  - `USE_COEFFICIENTS`: enable support for coefficients in a prime field
  - `EXPLICIT_CENTERS`: minimal enclosing sphere algorithm explicitly computes the center (can increase efficiency for low ambient dimensions)



The following options can be passed in the file name or as program argument:

  - `--dim k`: compute persistent homology up to dimension *k* (default *k*=1).
  - `--modulus p`: compute homology with coefficients in the finite prime field Z/*p*Z (when built with `USE_COEFFICIENTS`, default *p*=2).
  - `--threshold t`: compute Čech filtration up to radius *t* (default: no threshold ).





[Dionysus]: <http://www.mrzv.org/software/dionysus/>
[GUDHI]: <http://gudhi.gforge.inria.fr>
[Ripser]: <https://github.com/Ripser/ripser>
[CGAL]: <https://github.com/CGAL/cgal>
[Miniball]: <https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html>
[examples]: <https://github.com/s-clausen/cecher/tree/main/examples>


