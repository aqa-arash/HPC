#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iomanip>
using namespace std;

// Grid class to control 2D arrays
class Grid
{
public:
  // allocate memory in constructor
  Grid(int __nx__, int __ny__)
      : nx(__nx__)
      , ny(__ny__)
  {
      data = std::vector<double>(nx * ny); // Initialize vector with size
  }

  // free memory in destructor (not required for std::vector)
  ~Grid() {} // Empty destructor

  // Quick assign method
    void quickset(std::vector<vector<double>> input)
    {   int k=0;

        for (int j=0;j<ny;j++) {
            for (int i=0;i<nx;i++) {
                data[k]=input[i][j];
                k++;
            }
        }
    }

    //even quicker - unity
     void evenquickerset(double input)
    {
        #pragma omp parallel for
        for (int k=0;k<ny*nx;k++) {
           data[k]=input;
        }
    }
  // getter
  double operator()(int i, int j) const { return data[i + j * nx]; }

  // setter
  double& operator()(int i, int j) { return data[i + j * nx]; }

  // access dimensions- Domain size is considered 1x1
  int nyp() const { return ny; }//number of points in y direction
  int nxp() const { return nx; }// number of points in x direction
  double hx() const {return 1.0/(nx-1);}
  double hy() const {return 1.0/(ny-1);}


private:
  // internal data
  int nx, ny;
  std::vector<double> data;
};

//Gauss seidel function, changes the value of initial estimate x, based on B,A and number of iterations. where B is a stencil
void GSS(Grid& x,const Grid& b,int iter,const Grid& A){
//solve Ax=b with Gauss Seidel method

double C,W,E,N,S,NW,NE,SW,SE;
// calculating the stencil
SW=A(0,0)*(1.0/(x.hx()*x.hy()));
S= A(1,0)*(1.0/(x.hy()*x.hy()));
SE=A(2,0)*(1.0/(x.hx()*x.hy()));
W=A(0,1)*(1.0/(x.hx()*x.hx()));
E=A(2,1)*(1.0/(x.hx()*x.hy()));
NW=A(0,2)*(1.0/(x.hx()*x.hy()));
N=A(1,2)*(1.0/(x.hy()*x.hy()));
NE=A(2,2)*(1.0/(x.hx()*x.hy()));
C=(A(1,1)*((1.0/(x.hx()*x.hx()))+(1.0/(x.hy()*x.hy()))));

//main iteration of GS
for (int k=0;k<iter;k++){ //main iteration
    //------red
    int i,j;
    #pragma omp parallel for private(i)
    for ( j=1; j<x.nyp()-1;j++)// y direction iteration
        {
            if(j%2==0){
                   i=1; //jumping over
            }
            else {
                i=2; // switching the jump
            }

        for (;i<x.nxp()-1;i+=2){//iterate x direction
            double newparts=0; // repository for the sum values
            newparts+=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1); //calculate based on old estimates
            x(i,j)=((b(i,j)- newparts) /C); //calculate new estimate of xi
            }
        }

        //-----black
        #pragma omp parallel for private(i)
        for (j=1; j<x.nyp()-1;j++)// y direction iteration
        {
            if(j%2!=0){
                   i=1; //jumping 1 point
            }
            else {
                i=2; // switching on the next line
            }

        for (;i<x.nxp()-1;i+=2){//iterate x direction
            double newparts=0; // repository for sum
            newparts+=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1); //calculate based on old estimates
            x(i,j)=((b(i,j)- newparts) /C); //calculate new estimate of xi
            }
        }}
    }

// calculates the norm of the residual error, based on x,b and A
double Residualnorm (const Grid& x,const Grid& b,const Grid& A){
    double C,W,E,N,S,NW,NE,SW,SE;
    double er=0;
    int n=x.nxp()*x.nyp();//number of grid nodes
    // calculating the stencil
    SW=A(0,0)*(1.0/(x.hx()*x.hy()));
    S= A(1,0)*(1.0/(x.hy()*x.hy()));
    SE=A(2,0)*(1.0/(x.hx()*x.hy()));
    W=A(0,1)*(1.0/(x.hx()*x.hx()));
    E=A(2,1)*(1.0/(x.hx()*x.hy()));
    NW=A(0,2)*(1.0/(x.hx()*x.hy()));
    N=A(1,2)*(1.0/(x.hy()*x.hy()));
    NE=A(2,2)*(1.0/(x.hx()*x.hy()));
    C=(A(1,1)*((1.0/(x.hx()*x.hx()))+(1.0/(x.hy()*x.hy()))));

    #pragma omp parallel for reduction(+:er)
    for (int j=1; j<x.nyp()-1;j++)// y direction iteration
        {
        for (int i =1 ;i<x.nxp()-1;i++){//iterate x direction
            double newparts=0;// repository for the Ax
            newparts=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1)+C*x(i,j); //calculate based on old estimates
            newparts=(b(i,j)-newparts); //adds the second degree of error term for each point
            er+=newparts*newparts;//r square sum
            }            }
        //returns square root of the sum above, divided by the number of grid points.
        return sqrt(er/n);}


int main(int argc,char** argv) //can be used with command line variables, optionally, can also define threads
{//input pharsing
    int nx,ny,iter,threadcount;
    if (argc==4) {
    nx=stoi(argv[1]);
    ny=stoi(argv[2]);
    iter=stoi(argv[3]);
    }
    else if (argc==5){
    nx=stoi(argv[1]);
    ny=stoi(argv[2]);
    iter=stoi(argv[3]);
    omp_set_dynamic(0);
    threadcount=stoi(argv[4]);
    omp_set_num_threads(threadcount);
    }
    else { // if no command line arguments are used, the code does not crash or stop. it asks for them
        cout<< "You can choose the number of nodes, iteration, and thread count by running the software with command parameters"<<endl;
        cout<< "Nx, Ny, Iter, Thread count.  The last input is optional."<<endl;
        cout<< "Would you like to input the variable manually now? Y/N"<<endl;
        char A;
        cin >> A;
        if (A=='Y' || A=='y') {
            cout << "\n Nx = ";
            cin >>nx;
            if (cin.fail()) {
                    cout << "That was not an integer." << endl;
                    return -1;
                }
        cout << "\n Ny = ";
        cin >>ny;
        if (cin.fail()) {
                    cout << "That was not an integer." << endl;
                    return -1;
                }
        cout << "\n Number of iterations Iterations = ";
        cin >>iter;
        if (cin.fail()) {
                    cout << "That was not an integer." << endl;
                    return -1;
                }
        cout << "\n Number of threads (choose 0 to assign dynamically) = ";
        cin >>threadcount;
                if (cin.fail()) {
                    cout << "That was not an integer." << endl;
                    return -1;
                }
                else if (threadcount>0) {
                    omp_set_num_threads(threadcount);
                }
    }
    else {return -1;}
}

Grid A(3,3); //define variables as Grid A for the stencil
Grid x(nx,ny);
Grid b(nx,ny);
double er; // for error
A.quickset({{0,-1,0},{-1,2,-1},{0,-1,0}});// Coefficient matrix
x.evenquickerset(0);//x vector (initial estimate)
b.evenquickerset(0);//b vector (right hand side)
double pi=3.141592653589793238462643383279502884;

//define the boundry condition - basically g function
//in x direction
#pragma omp parallel for
for(int i=0;i<x.nxp();i++ ){
    x(i,x.nyp()-1)=-cos(pi*i*x.hx());
    x(i,0)=cos(pi*i*x.hx());
}
#pragma omp parallel for
// in y direction
for(int i=0;i<x.nyp();i++){
    x(0,i)=cos(pi*i*x.hy());
    x(x.nxp()-1,i)=-cos(pi*i*x.hy());
}
//define the b for the internal point -- basically f function
#pragma omp parallel for
for(int j=0;j<x.nyp();j++ ){
        for(int i=0; i<x.nxp();i++){
            b(i,j)= 2*(pi*pi)*cos(pi*i*x.hx())*cos(pi*j*x.hy());
        }
}
    // start time of the actual GS function
    auto start= std::chrono::system_clock::now();
    GSS(x,b,iter,A); //calling the GS function
    auto fin= std::chrono::system_clock::now(); // End time of actual GS function and error calculations
    er=Residualnorm(x,b,A); //calculate residual norm
    auto ellapse= (fin-start);//wall clock time of GS function
    cout<< "calculation complete for "<< iter<< " iterations, with grid size of "<<x.nxp()<< "x" << x.nyp()<< endl;
    cout << setw(20)<< "run time (s)" << setw(20)<< "Error norm"<< endl;
    cout <<setw(20)<<ellapse.count()/pow(10,9)<<setw(20)<<er<< endl;

/* To facilitate testing, outputs test results to file.
    ofstream errorfile;
    errorfile.open ("error.txt",std::ios_base::app);
    //errorfile << setw(20)<< "run time" << setw(20)<< "Error norm"<< setw(20)<<"gridsize"<<setw(20)<<"number of cores"<<endl;
    errorfile <<setw(20)<<ellapse.count()/pow(10,9)<<setw(20)<<er<< setw(20)<<nx<<"*"<<ny<<setw(20)<< threadcount <<endl;
    errorfile.close();
*/
    //writes the solution to the text file
    ofstream solutionfile;
    solutionfile.open ("solution.txt");
    for (int j=x.nyp()-1; j>=0 ;j--) { //printing results of x
        for (int i=0;i<x.nxp();i++){
               solutionfile << i*x.hx() <<" " <<j*x.hy() << " "<<x(i,j)<<"\n";
        }
    }
    solutionfile.close();
    cout<<"The result is available in the file solution.txt \n";


    return 0;
}
