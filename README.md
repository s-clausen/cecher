![Cecher](logo.png)

Copyright © 08.07.2023 Sönke Clausen

The Cecher is a C++ software for efficient computation of Čech persistence barcodes. It outperforms current software ([Dionysus], [GUDHI]) and is especially adapted to handle high ambient dimension. Its main features are:

  - a core persistence algorithm adapted from the C++ code of the ultrafast [Ripser], featuring its support for coefficients in prime fields, its memory-efficient philosophy, and use of a combination of cohomology, clearing, and zero-apparent pairs shortcut.
  - a version of the union-find algorithm in [Ripser] with compression.
  - a symbolic perturbation of the Čech filtration, allowing to completely skip the construction of many columns of the coboundary matrix corresponding to zero-persistence pairs.
  - a novel minimal enclosing sphere implementation, outperforming in the context of a persistence algorithm established software like [Miniball] with a construction using the Euclidean distance matrix and featuring efficient computation in high ambient dimension.
  - a lazy-exact type implementation, where computation are performed primarily with an interval arithmetic type.
    

### Assumption for Efficient Computation

The Cecher is robust to degenerate input, but requires exact computations when flat simplices exist or when minimal circumspheres of simplices have coinciding radii. However on most 'real world data' all operations utilize the interval arithmetic type.


### Input Format

Currently supported is an input in the form of a rational point cloud, where numerator and denominator are seperated by a forward-slash, and values are seperated by a whitespace ([examples]).


### Dependencies 

The implementation currently utilizes the (double-precision) interval arithmetic type and exact rational type of [CGAL].


### Options

The Cecher supports the following preprocessor macros:

  - `USE_COEFFICIENTS`: enable support for coefficients in a prime field Z/*p*Z (default: *p*=2).
  - `EXPLICIT_CENTERS`: enable variant of the minimal enclosing sphere algorithm for low ambient dimension


The Cecher supports the following program arguments (can be passed in the file name):

  - `--dim k`: compute persistent homology up to dimension *k* (default: *k*=1).
  - `--threshold t`: compute Čech filtration up to radius *t* (default: no threshold).
  - `--modulus p`: compute homology with coefficients in the prime field Z/*p*Z (requires `USE_COEFFICIENTS`). 






[Dionysus]: <http://www.mrzv.org/software/dionysus/>
[GUDHI]: <https://gudhi.inria.fr/>
[Ripser]: <https://github.com/Ripser/ripser>
[CGAL]: <https://github.com/CGAL/cgal>
[Miniball]: <https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html>
[examples]: <https://github.com/s-clausen/cecher/tree/main/examples>


