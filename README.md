![Cecher](logo.png)

Copyright © 08.07.2023 Sönke Clausen

The Cecher is a C++ software for efficient computation of Čech persistence barcodes outperforming current software. Its main features are:

  - a core persistence algorithm adapted from the C++ code of the ultrafast [Ripser], featuring its memory-efficient philosophy, and use of a combination of cohomology, clearing, and zero-apparent pairs shortcut.
    
  - a symbolic perturbation of the Čech filtration, allowing to completely skip the construction of many columns of the boundary matrix corresponding to zero-persistence pairs.
    
  - a minimal enclosing sphere algorithm especially efficient in high ambient dimension, outperforming established software like [Miniball] in the context of a persistence algorithm with a construction using the distance matrix.
    
  - a lazy-exact implementation, where computations are performed primarily with an interval arithmetic (double-precision) type.
    

### Input

The Cecher currently supports an input in the form of a rational point cloud, where numerator and denominator are seperated by a forward-slash, and values are seperated by a whitespace ([examples]).

### Lazy-Exact Implementation

The Cecher falls back on exact computations when flat simplices exist or when minimal circumspheres of simplices have coinciding radii. On most 'real world data' all operations utilize the interval arithmetic (double-precision) type.

### Dependencies 

The Cecher depends on the library [CGAL] for its lazy-exact implementation. 

### Options

The Cecher supports the following preprocessor macros:

  - `USE_COEFFICIENTS`: enable support for coefficients in a prime field Z/*p*Z (default: *p*=2).
  - `EXPLICIT_CENTERS`: enable variant of the minimal enclosing sphere algorithm for low ambient dimension.


The Cecher supports the following program arguments (can be passed in the file name):

  - `--dim k`: compute persistent homology up to dimension *k* (default: *k*=1).
  - `--threshold t`: compute Čech filtration up to radius *t* (default: no threshold).
  - `--modulus p`: compute homology with coefficients in the prime field Z/*p*Z (requires `USE_COEFFICIENTS`). 





[Ripser]: <https://github.com/Ripser/ripser>
[CGAL]: <https://github.com/CGAL/cgal>
[Miniball]: <https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html>
[examples]: <https://github.com/s-clausen/cecher/tree/main/examples>


