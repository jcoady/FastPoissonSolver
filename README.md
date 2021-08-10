# FastPoissonSolver
This repository contains examples of using the Fastest Fourier Transform in the West (FFTW) library. 

1) There is an example of a 2D Fast Poisson Solver in the file.

   Poisson2D_FFTW_type1_DST_Strang.cpp
   
Compile the file on the command line using

   g++ -O2 Poisson2D_FFTW_type1_DST_Strang.cpp -o Poisson2D_FFTW_type1_DST_Strang -lfftw3 -Wall

and then run it with the command

 ./Poisson2D_FFTW_type1_DST_Strang
 
 Unpon running the program you shoould see a result of an error between the exact solution and the FFT solution to around 15 decimal places.
 
2) There is an example of a 3D Fast Poisson Solver in the file.

   Poisson3D_FFTW_type1_DST_Strang.cpp
   
Compile the file on the command line using

   g++ -O2 Poisson3D_FFTW_type1_DST_Strang.cpp -o Poisson3D_FFTW_type1_DST_Strang -lfftw3 -Wall

and then run it with the command

 ./Poisson3D_FFTW_type1_DST_Strang
 
Unpon running the program you shoould see a result of an error between the exact solution and the FFT solution to around 15 decimal places.
   
3) There is a 3D example of a program that executes a Gauss-Seidel solver in 3D.

    Gauss_Siedel_3D.cpp

Compile the file on the command line using

   g++ -O2 Gauss_Siedel_3D.cpp -o Gauss_Siedel_3D
  
and then run it with the command

 ./Gauss_Siedel_3D
 
 This will display a number of results comparing the exact solution to the Gauss-Siedel solver's solution.
 
 Compare the results of the output of the Gauss_Siedel_3D program with the Poisson3D_FFTW_type1_DST_Strang program. They should give the same results to several decimal places.
 

