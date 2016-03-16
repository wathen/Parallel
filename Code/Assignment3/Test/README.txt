README.txt: source files for cs418 hw3.

The files:
  Makefile      a makefile
  examples.cu   some simple CUDA examples
  perc.cu       percolation simulation
  time_it.cu    utilities for measuring the time to execute code
  time_it.h     header file corresponding to time_it.cu
  timing.cu     simple demo of time_it

Some examples:
  % make timing
  % ./timing empty 100
  empty kernel: mean(T) =  6.640e-06, stddev(T) =  7.980e-07
    The empty kernel does nothing.  It launches one block with one thread that
    immediately exits.  This measures the time to coordinate between the CPU
    and GPU for a kernal launch.  Apparently, that's 6 or 7 microseconds.

  % ./timing sync 100
  synchronize threads: mean(T) =  1.549e-04, stddev(T) =  5.551e-07
    Each thread in this kernel calls __syncthreads() 1000 times.  There
    are 1024 threads in the block.  Apparently, it takes about 155ns
    to perform a __syncthreads().  This means that __syncthreads() is
    about 40 times faster than launching a new kernel.  In other words,
    it is much faster to synchronize within a block than across blocks.

  % make perc
  % ./perc 1024 65536 0.71
  371 vertices reachable after 65536 generations
  perc(1024, 65536, 0.710): mean(T) =  2.747e-01, stddev(T) =  3.038e-03
    Compute 65536 generations of percolation on a graph of width 1024.
    This means that the GPU is computing about 3.8 million vertex updates
    per second.
  % ./perc 1024 65536 0.1
  0 vertices reachable after 65536 generations
  perc(1024, 65536, 0.100): mean(T) =  5.508e-02, stddev(T) =  2.243e-05
    With the lower marking probability, all vertices are unreachable after
    just a few generation.  This eliminates most of the calls to the
    random number generator, and the code runs much faster even if the
    resulting graph is boring.  This suggest that we should focus on speeding
    up the random number generator.

  % make examples
  % ./examples 1000000 4
  The GPU is a GeForce GTX 550 Ti
  Cuda capability 2.1.
  x = 0.123, 0.410, 0.919, 0.282, 0.770, 0.673, 0.836, 0.520, 0.948, 0.186
  y = 0.548, 0.966, 0.128, 0.434, 0.958, 0.156, 0.514, 0.974, 0.098, 0.345
  z = 380142.000
  The results match!
    Run test-case 4 (dot-product) on vectors of length 1000000.
    The code prints the first 10 elements of the x and y vectors and their
    dot-product.  The result from the GPU is compared with the result from a
    sequential implementation in C to see if they are "near".  Floating-point
    arithmetic has terrible round-off error.  So, the results only match to
    about a part in 10^4.  See the function "near" in examples.cu for more
    details.
