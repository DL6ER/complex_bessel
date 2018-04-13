# Library: Complex Bessel and Hankel functions
A library to efficiently evaluate Bessel and Hankel functions of all kinds in C. This library includes full support for complex input and output parameters.

## Introduction

This project provides a standard C interface for the library available [here](https://github.com/valandil/complex_bessel).

With the provided library interface you can compute

1. Bessel functions of first and second kind of any order (also negative orders)
2. Hankel functions of first and second kind of any order (also negative orders)
3. Derivatives of any order of the functions described in 1. and 2.

## Installation

1. Clone this project to your machine
  ```bash
  git clone https://github.com/DL6ER/complex_bessel.git
  cd complex_bessel
  ```

2. Install and use `cmake`
  ```bash
  sudo apt-get install cmake
  bash build.sh
  ```

3. Install library (the library will be installed to `/usr` by default)
  ```bash
  cd build
  sudo make install
  cd ..
  ```

4. Run test program to verify the library is working correctly
  ```bash
  gcc test.c -lcomplex_bessel -lm
  ./a.out
  ```

## Implemented routines

### Bessel functions of first and second kind
![Bessel function of first kind](https://github.com/DL6ER/complex_bessel/blob/master/img/bessel-1.svg)

![Bessel function of second kind](https://github.com/DL6ER/complex_bessel/blob/master/img/bessel-2.svg)
```c
double complex bessel(double order, int kind, double complex z);
double complex besselJ(double order, double complex z);
double complex besselY(double order, double complex z);
```
These routines compute the Bessel functions for any given order. `besselJ` (first kind) and `besselY` (second kind) are provided as short handles for the more generic `bessel` function.

### Hankel functions of first and second kind
![Hankel functions of first and second kind](https://github.com/DL6ER/complex_bessel/blob/master/img/hankel.svg)
```c
double complex hankel(double order, int kind, double complex z);
double complex hankelH1(double order, double complex z);
double complex hankelH2(double order, double complex z);
```
These routines compute the Hankel functions for any given order. `hankelH1` (first kind) and `hankelH2` (second kind) are provided as short handles for the more generic `hankel` function.

### Derivation functions
```c
double complex besselDiff(double order, int kind, int n, double complex z);
double complex besselJdiff(double order, int n, double complex z);
double complex besselYdiff(double order, int n, double complex z);
double complex hankelDiff(double order, int kind, int n, double complex z);
double complex hankelH1diff(double order, int n, double complex z);
double complex hankelH2diff(double order, int n, double complex z);
```
These routines compute the n-th derivative of the Bessel/Hankel functions for any given order. `besselJdiff` (Bessel, first kind), `besselYdiff` (Bessel, second kind), `hankelH1diff` (Hankel, first kind), and `hankelH2diff` (Hankel, second kind) are provided as short handles for the more generic `besselDiff` and `hankelDiff` functions.

## License
This project is copyrighted under the latest version of the European Union Public License (EUPL).
Please see `LICENSE` file for your rights under this license.
