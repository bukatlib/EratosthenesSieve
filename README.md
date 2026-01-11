# EratosthenesSieve

Parallel segmented Eratosthenes sieve in Modern C++.

## Goals of the implementation

1. Educational purposes (sieve has just about 1k lines of code).
2. Efficient, only about 19% to 45% slower than [[prime sieve](https://github.com/kimwalisch/primesieve)].
3. Possible to store/load sieve binary image to file (all in memory).

## Technology

- C++20 bit manipulation (rotl, popcount)
- Parallelization with OpenMP (dynamic schedule).

TODO
