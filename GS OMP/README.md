**README**

**Gauss-Seidel Solver** 

**Introduction**

This program solves the Poisson equation using the Gauss-Seidel iterative method. Computations for an approximate solution to a discretization done by an elliptic partial differential equation: 

` `−Δu = f in Ω and u = g on ∂Ω, 

Where f is f(x, y) = 2π2 cos(πx) cos(πy) and, g is g(x, y) = cos(πx) cos(πy).

**Functionality**

- **Grid Class**: The **Grid** class creates 2D arrays and provides methods for setting and accessing grids’ values. This class has also 2 functionalities; the *quick assign method* and the *even quicker assign method*. Thanks to these 2 inner methods, it is possible to create A, x, and b vectors and matrices.
- **Gauss-Seidel Solver**: The **GSS** function solves the Poisson equation 𝐴𝑥=b using the Gauss-Seidel method. In detail, this function assign the formula of the stencil grid as;

  N= North

  E= East

  S=South

  W= West

  NE= Northeast

  NW= Northwest

  SE=Southeast

  SW=Southwest

  C=Central

  Even though the problem uses a stencil that has values on N,E,W,S, and C, it is also possible to use different stencils thanks to this structure. As a result, this function solves the Gauss-Seidel problem with red and black properties. This calculation is performed by different loops that are parallelized.

- **Residual Calculation**: The **Residual** function calculates the  L<sub>2</sub>  norm of the residual r=AU-b after the final iteration.  This calculation is also parallized by OpenMP using a reduction statement.
- **Main Function**: The **main** function initializes the problem, sets boundary conditions, calls the solver, and saves results.

  This function begins with observing the arguments which have to be defined by users. In the first part of the function, there is a communication with users in case of failure of defining these arguments. 

  Secondly, the initialization of the A, x, and b is done by the using Grid class and the definition of the boundary condition and inner points calculation is done. 

  Lastly, the Gauss-Seidel Solver and Residual calculation functions are called in order to solve the problem. It should be noted that, this part uses the *chrono* library in order to calculate the time that is spent by Gauss-Seidel and Residual functions.  Additionally, the results are written in the text file for possible analysis. 

**Usage**

- **Command Line Arguments**: The program can be executed with command line arguments specifying grid size, number of iterations, and optionally, the number of threads.
  - Example: **./gssolve 100 100 100**

    The first argument is the number of points in the x-direction.

    The second argument is the number of points in the y-direction.

    The third argument is the number of iterations.

- **Interactive Mode**: If no command line arguments are provided, the program prompts the user to input parameters interactively. This part is coded in the main function. In this mode, users can also provide the number of threads.

Example: **./gssolve 100 100 100 16** 

The first argument is the number of points in the x-direction.

The second argument is the number of points in the y-direction.

The third argument is the number of iterations.

The fourth argument is the number of threads or processors.

**Compilation**

- **Dependencies**: The program requires a C++ compiler with support for OpenMP. It is recommended to use the g++ compiler in this project. In some cases, the CodeBlock IDE makes it difficult to use OpenMP on computers with Windows operating systems. In this case, it is recommended to run it on a Linux node.
- **Compilation Command**: In order to compile the code, the following command could be used.
  - Example: -O3 -fopenmp -Wall -Winline -Wshadow -std=c++17

**Output**

- **Solution File**: The solution is saved in the **solution.txt** file, containing the calculated values of 𝑥 at each grid point.
- **Error File**: Error norms, runtime, grid size, and thread count are logged in the **error.txt** file. This file allows users to analyze the accuracy of the solver. The structure of this file is described as follows.
  - **Run Time**: This column indicates the time taken for a particular computation or process to execute. It's measured in seconds.
  - **Error Norm**: This column represents the level of error in the computation or simulation. It's a measure of the deviation from the expected or true value.
  - **Grid Size**: This column specifies the size of the computational grid used in the code. 
  - **Number of Cores**: This column indicates the number of processor cores or computing units utilized during the computation. It's a measure of parallelization.

**Parallelization**

- **OpenMP Parallelization**: Parallelization is implemented using OpenMP directives to distribute computation across multiple threads. 

  Calculation of Gauss-Seidel Solver and Residual functions are parallelized with OpenMP.

**Performance**

In order to analyze the performance of the code, some experiments are done. This code is run in 34\*34 and 1026\*1026 matrix sizes using 100 iterations. In addition, different processor numbers are used in these experiments and, the results such as parallel efficiency, speed-up, and wall-clock time values are stated in the report file.

Additionally, it is also possible to use bash.txt. file in order to analyze the performance data. 


**Authors**

This program was authored by 

Anish Adgaonkar @ FAU

Arash Moradian @ FAU

Mehmet Arif Bagci @FAU  

**Licence** 

This project is open-source. It is free to use, modify, and distribute this software.

