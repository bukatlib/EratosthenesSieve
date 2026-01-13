# EratosthenesSieve

Parallel segmented Eratosthenes sieve in Modern C++ that is cache-friendly
and uses wheel 2, 2x3, or 2x3x5 based on the build configuration.

## Goals of the implementation

1. Educational purposes (sieve has just about 1k lines of code).
2. Efficient, only about 13% to 45% slower than primesieve [[1](https://github.com/kimwalisch/primesieve)].
3. Possible to store/load sieve binary image to/from file.
4. A good CPU/memory stress test (I detected a faulty hardware based on the wrong prime count).

## Technology

- C++20 bit manipulation (rotl, popcount, and countr\_zero)
- Parallelization with OpenMP (dynamic schedule)
- CMake build system

## How to build

```
$ mkdir build
$ cd build
$ cmake ..
$ make -j4
```

As an example of a more complex setup consider a usage of a different compiler,
Profile Guided Optimization (PGO), and wheel 2x3 instead of 2x3x5 default.

```
$ mkdir build_pgo_gcc
$ cd build_pgo_gcc
$ cmake -DCMAKE_C_COMPILER=/usr/bin/gcc-12 -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 -DWHEEL_TYPE="2x3" -DPGO_STAGE="Collect" ../
$ make -j5 profile
$ cmake -DPGO_STAGE="Optimize" ../
$ make -j5
```

PGO is only supported for GCC compilers. Note that it might not improve the
performance, for example, GCC 14 produces slower code with PGO on Raspberry Pi 5
but GCC 12 achieves a slight speedup, therefore, the best is to experiment.

## Usage examples

```
$ ./prime_sieve_bench 100000000000 -
sieve size: 26666666666 bits (3333333336 bytes, 3178.91 MBs)
Sieve initialized in 0.61024 secs
Primes sieved in 11.0368 secs.
Found 4118054813 primes in 0.270697 secs (bandwidth 12313.9 MB/s).

$ ./prime_sieve_bench 100000000000 sieve_img.bin
sieve size: 26666666666 bits (3333333336 bytes, 3178.91 MBs)
Sieve initialized in 0.596821 secs
Primes sieved in 10.9931 secs.
Found 4118054813 primes in 0.267827 secs (bandwidth 12445.9 MB/s).
Primes written to 'sieve_img.bin' file in 2.12945 secs.

$ ./prime_sieve_bench sieve_img.bin
sieve size: 26666666666 bits (3333333336 bytes, 3178.91 MBs)
Sieve up to 100000000000 loaded in 0.884763 secs
Found 4118054813 primes in 0.266247 secs (bandwidth 12519.7 MB/s).

$ ./prime_sieve_bench 1000 primes_1k.txt
sieve size: 266 bits (40 bytes, 3.8147e-05 MBs)
Sieve initialized in 0.00300503 secs
Primes sieved in 8.58307e-05 secs.
Found 168 primes in 2.19345e-05 secs (bandwidth 1.82361 MB/s).
Primes written to 'primes_1k.txt' file in 0.000384569 secs.
```

The above results are on Raspberry Pi5 16 GB (GCC 14 without PGO, wheel 2x3x5).
There are also `show_cpu_caches` and `factorization_wheel_demo` binaries,
the first one extract CPU caches info from `/sys` directory (Linux only)
and the second is a small demonstration of wheel 2x3x5 (removal of 2, 3,
and 5 multiples).

## Performance Notes

The most important optimization is to process parts of sieve in parallel
where each part is a number of segments (`WORKCHUNK_SEGMENTS`) where the
start indices are calculated at the beginning (expensive divides) and then
each segment is processed separately. Segment should fit to the cache, which
level of cache and fraction can be adjusted in `Eratosthenes::Eratosthenes`
or default `DEFAULT_SEGMENT_SIZE` if cache info is not available.

The best on Raspberry Pi5 is to use L1d cache with 64 kB size, compared
to Ryzen 7950X, where the segment size is 0.75 of L2d cache size.
Note that AMD SMT results in 2 HW threads per core, therefore,
each physical core handles two segments of 1.5 L2d cache size in total.
In other words, AMD 7950X manages to feed the data partially from L3 cache
without significant deterioriation in Instructions Per Cycle (IPC).

There are three core methods to sieve, i.e., reset composites:

1. `Eratosthenes::sieve_bit_range_small`
2. `Eratosthenes::sieve_bit_range_medium`
3. `Eratosthenes::sieve_bit_range_large`

The suffix indicates whether the sieving step is small, medium, or large with
respect to the segment size. A small step usually means that multiple composites
can be reset in one uint64\_t element, therefore, we prepare reset masks and
use SIMD to efficiently sieve with small primes. The method for sieving with
medium primes expect multiple full loops per segment, therefore, the loop
unrolling is applied and small work before the main loop is acceptable.
Sieving with large primes is an optimized method where only a few resets
are done per segment.

Note that the wheel 2x3x5 employs steps compression to better fit `sqrt(N)`
sieve records in L2/L3 caches (beneficial effects on Rasperry Pi 5).

## Experimental results

You can find the experimental results for Raspberry Pi 5 16 GB and
AMD Ryzen 7950X in `results` directory. There are performance metrics,
interesting fragments of assembly (hotspots), and prime counts.

## Limitations

There are currently the following limitations:

1. Detection of cache size is Linux only (contribution welcomed).
2. Binary image of sieve is not Endianess aware (portability).
3. Sieving of large ranges could be better optimized.

## Contributions

Contributions are welcomed :-).

You might contribute experimental results, improve performance/documentation,
add Windows/MacOS support...
