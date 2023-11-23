![Cecher](logo.png)

Copyright © 08.07.2023 Sönke Clausen

The Cecher is a C++ software for efficient computation of Čech persistence barcodes outperforming current software ([Dionysus], [GUDHI]). Its main features are:

  - a core persistence algorithm adapted from the C++ code of the ultrafast [Ripser], featuring its memory-efficient philosophy, and use of a combination of cohomology, clearing, and zero-apparent pairs shortcut.
    
  - a symbolic perturbation of the Čech filtration, allowing to completely skip the construction of many columns of the boundary matrix corresponding to zero-persistence pairs.
    
  - a minimal enclosing sphere algorithm especially efficient in high ambient dimension, outperforming established software like [Miniball] in the context of a persistence algorithm with a construction using the distance matrix.
    
  - a lazy-exact implementation, where computations are performed primarily with an interval arithmetic type.

  - a version of the union-find algorithm in [Ripser] with compression.


    
### Input & Options

The Cecher supports an input in the form of a decimal point cloud, where a whitespace seperates values ([examples]). The input is given in a file whose name is passed as a command line argument. 

The following preprocessor macros are supported:

  - `USE_COEFFICIENTS`: enable support for coefficients in a prime field Z/*p*Z (default: *p*=2).
  - `EXPLICIT_CENTERS`: enable variant of the minimal enclosing sphere algorithm for low ambient dimension.

The following command line arguments are supported:

  - `--dim k`: compute persistent homology up to dimension *k* (default: *k*=1).
  - `--threshold t`: compute Čech filtration up to radius *t* (default: no threshold).
  - `--modulus p`: compute homology with coefficients in the prime field Z/*p*Z (requires `USE_COEFFICIENTS`).

### Lazy-Exact Implementation

The Cecher falls back on exact computations when flat simplices exist or when minimal circumspheres of simplices have coinciding radii. On most 'real world data' all operations utilize the interval arithmetic (double-precision) type.

### Dependencies 

The Cecher depends on the library [CGAL] for its lazy-exact implementation. 

### How to Build & Run

On a preconfigured system (Windows, macOS or Linux, where [CGAL] and its dependencies are installed) building and running the Cecher should be as easy as:

```sh
cd path\to\cecher
cmake -DUSE_COEFFICIENTS=off -DEXPLICIT_CENTERS=off .
make
.\cecher examples\dragon.txt --dim 2

```
Adapt preprocessor macros and program arguments to your needs.


### Credits 

Thank you to Ulrich Bauer and Fabian Roll for their contributions.


[Ripser]: <https://github.com/Ripser/ripser>
[CGAL]: <https://www.cgal.org/download.html>
[Dionysus]: <http://www.mrzv.org/software/dionysus/>
[GUDHI]: <https://gudhi.inria.fr/>
[Miniball]: <https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html>
[examples]: <https://github.com/s-clausen/cecher/tree/main/examples>


