MPI Parallel CG Solver for PDEs

Overview
This project implements a parallel Conjugate Gradient (CG) solver for solving Partial Differential Equations (PDEs) using MPI for inter-process communication. It employs a Cartesian communicator to manage a grid of processes and exchange boundary data, optimizing the computation of PDE solutions over a 2D grid.

Requirements
•	C++ compiler – mpic++ compiler
•	MakeFile
•	Grid.hpp header file
•	MPI Library

Functionality
The main program cgsolve.cpp consists of 3 basic functions.
1)	Getdims : This function computes the optimal dimensions for the Cartesian communicator based on the world size and domain dimensions.
2)	Updateneighbours : This function handles the communication of boundary data between neighboring processes using non-blocking send and receive operations.
3)	main : This function initializes MPI, sets up the grid and Cartesian
communicator, manages the grid partitioning, and runs the CG solver iteration


Dependencies
In order to use the program properly, the use Grid.hpp is mandatory.
1) Grid.hpp :The Grid class is a utility to manage and operate on 2D arrays in a structured manner. It provides functionalities for:

•	Memory allocation and initialization.
•	Assigning values to the grid using quick set methods.
•	Boundary condition setting.
•	Data access through overloaded operator() for both getting and setting values.

This header file consists 8 functionalities:

1)	Pointwisesum:Adds elements of alpha * B to corresponding elements in A to compute the pointwise sum.
2)	NodeResidualnorm & Residualnorm: Computes the norm of the residual error between the computed solution and the true solution.
3)	WriteSolution: Writes the solution to a text file, including boundary conditions.
4)	Residual: Computes the residual error for the linear system.
5)	Errornorm: Calculates the error norm for a given solution.
6)	Stencil_Grid_Mult:Performs stencil grid multiplication for calculating the residual error.
7)	Grid_dot_product: Computes the dot product of two grids, with an option for transpose.
8)	L2norm: Calculates the L2 norm of the grid data.



Usage
The program accepts command-line arguments for setting parameters:

mpirun -np N ./cgsolve nx ny iter print_status	


•	nx= number of nodes in x direction
•	ny =number of nodes in y direction
•	iter = number of iteration (0 to run until convergence)
•	print_status= integer value for writing solution.txt. (1 to print solution file)
•	N= number of processes - MPI

Command Line Arguments	: The program can be executed with command line
arguments directly.

Example: mpirun -np 2 ./cgsolve 18000 18000 4 0
•	nx= 18000
•	ny =18000
•	iter = 4
•	print_status= does not create the “solution.txt” file
•	N= 2 - MPI


Interactive mode:  If no command line arguments are provided, the program runs with the users input.

Example: mpirun -np 1 ./cgsolve
>>You can choose the number of nodes in x and y directions and print status repectively, by running the software with command parameters
>> Nx(0 to exit)
<<
>>Ny
<<
>>Number of iterations
<<

Output
•	solution.txt: Contains the final grid values of the solution after computation (if print status is 1.)
•	report.txt: It contains wall-clock time, number of nodes domain size. 

References
•	MPI Documentation: https://www.mpich.org/documentation
•	OpenMP Documentation: https://www.openmp.org/specifications
•	Conjugate Gradient Method: https://en.wikipedia.org/wiki/Conjugate_gradient_method

License
This project is open-source. It is free to use, modify, and distribute this software.

Authors
This program was authored by 
Anish Adgaonkar @ FAU
Arash Moradian @ FAU
Mehmet Arif Bagci @ FAU



 
 
