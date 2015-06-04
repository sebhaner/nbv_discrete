Discrete Next Best View Planning

This MATLAB code implements the method presented in
  S. Haner and A. Heyden, "Discrete Optimal View Path Planning", VISAPP 2015.
If you use this code in a publication, please cite the above paper.

To run the code, you need YALMIP (http://users.isy.liu.se/johanl/yalmip) and it is highly recommended
to also install the MOSEK convex optimizer, free for academic use.
The code also depends on a number of mex-files, some of which depend on the header library Eigen.
Binaries for Win64 are provided, on Linux you will need to compile them with

    mex *.c* -I"/path/to/eigen/"

To solve the example problems, run the file run_demo.m, or just the individual code cells within.

