#include <iostream>
#include <vector>
#include <cmath>
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
    {
        //#pragma omp parallel for
        for (int j=0;j<ny;j++) {
            for (int i=0;i<nx;i++) {
                data[i + j * nx]=input[i][j];
                }
        }
    }

    //even quicker - unity
     void evenquickerset(double input)
    {
        //#pragma omp parallel for
        for (int k=0;k<ny*nx;k++) {
           data[k]=input;
        }
    }

    //Boundary set
    void boundaryset (double x){
    for (int i=0;i<nx;i++){
        data[i + ny * nx]=x;
        data[i]=x;}
    for (int j=0;j<ny;j++){
        data[j*nx]=x;
        data[(j+1)*nx]=x;
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



//adds the elements of alpha * B, to corresponding elements in x
Grid Pointwisesum (Grid&A,const Grid& B,double aplha =1){
    Grid sum(A.nxp(),A.nyp());
    sum=A;
    //#pragma omp parallel for
    for (int j=1;j<A.nyp()-1;j++){
        for (int i=1;i<A.nxp()-1;i++){
            sum(i,j)=A(i,j)+ (aplha*B(i,j));
        }
    }
    return sum;
}

// calculates the norm of the residual error, based on x,b and A
double NodeResidualnorm (const Grid& x,const Grid& b,const Grid& A,double hx,double hy){
double C,W,E,N,S,NW,NE,SW,SE;
double er=0.0;
// calculating the stencil
SW=A(0,0)*(1.0/(hx*hy));
S= A(1,0)*(1.0/(hy*hy));
SE=A(2,0)*(1.0/(hx*hy));
W=A(0,1)*(1.0/(hx*hx));
E=A(2,1)*(1.0/(hx*hx));
NW=A(0,2)*(1.0/(hx*hy));
N=A(1,2)*(1.0/(hy*hy));
NE=A(2,2)*(1.0/(hx*hy));
C=(A(1,1)*((1.0/(hx*hx))+(1.0/(hy*hy))));


//#pragma omp parallel for reduction(+:er)
for (int j=1; j<x.nyp()-1;j++)// y direction iteration
    {
        for (int i =1 ;i<x.nxp()-1;i++){//iterate x direction
            double newparts=0.0;// repository for the Ax
            newparts=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1)+C*x(i,j); //calculate based on old estimates
            newparts=(b(i,j)-newparts); //adds the second degree of error term for each point
            er+=newparts*newparts;//r square sum
            }
    }
//returns square root of the sum above, divided by the number of grid points.
return er;
}

// calculates the norm of the residual error, based on x,b and A
double Residualnorm (const Grid& x,const Grid& b,const Grid& A,double hx,double hy){
double C,W,E,N,S,NW,NE,SW,SE;
double er=0.0;
// calculating the stencil
SW=A(0,0)*(1.0/(hx*hy));
S= A(1,0)*(1.0/(hy*hy));
SE=A(2,0)*(1.0/(hx*hy));
W=A(0,1)*(1.0/(hx*hx));
E=A(2,1)*(1.0/(hx*hx));
NW=A(0,2)*(1.0/(hx*hy));
N=A(1,2)*(1.0/(hy*hy));
NE=A(2,2)*(1.0/(hx*hy));
C=(A(1,1)*((1.0/(hx*hx))+(1.0/(hy*hy))));


//#pragma omp parallel for reduction(+:er)
for (int j=1; j<x.nyp()-1;j++)// y direction iteration
    {
        for (int i =1 ;i<x.nxp()-1;i++){//iterate x direction
            double newparts=0.0;// repository for the Ax
            newparts=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1)+C*x(i,j); //calculate based on old estimates
            newparts=(b(i,j)-newparts); //adds the second degree of error term for each point
            er+=newparts*newparts;//r square sum
            }
    }
//returns square root of the sum above, divided by the number of grid points.
return sqrt(er/(x.nxp()*x.nyp()));
}

//writes the solution to the text file
void Writesolution(const Grid &x,double hx,double hy, int globalIstart, int globalJstart,std::vector<int> dims,std::vector<int>my_coords){
    ofstream solutionfile;
    solutionfile.open ("solution.txt",std::ios_base::app);
double pi=3.141592653589793238462643383279502884;
                           for (int j=0; j<x.nyp() ;j++) { //printing results of x
                                if (j%2==0){
                                solutionfile << 0 <<" " <<(j+globalJstart)*hy << " "<< cos(pi*(j+globalJstart)*hy)<<"\n";
                                for (int i=1;i<x.nxp()-1;i++){
                                        solutionfile << i*hx <<" " <<(j+globalJstart)*hy << " "<<x(i,j)<<"\n";
                                        }
                                solutionfile << 1 <<" " <<(j+globalJstart)*hy << " "<< -cos(pi*(j+globalJstart)*hy)<<"\n";
                                }
                                else{
                                    solutionfile << 1 <<" " <<(j+globalJstart)*hy << " "<< -cos(pi*(j+globalJstart)*hy)<<"\n";
                                    for (int i=x.nxp()-1;i>0;i--){
                                        solutionfile << i*hx <<" " <<(j+globalJstart)*hy << " "<<x(i,j)<<"\n";
                                        }

                                    solutionfile << 0 <<" " <<(j+globalJstart)*hy << " "<< cos(pi*(j+globalJstart)*hy)<<"\n";

                                }
                                }

    solutionfile.close();
cout<<" Node " << my_coords[0]<<" " <<my_coords[1]<< " wrote the result in the file solution.txt \n";
}



// calculates the residual error, based on x,b and A
Grid Residual (const Grid& x,const Grid& b,const Grid& A,double hx,double hy){
double C,W,E,N,S,NW,NE,SW,SE;
Grid er(x.nxp(),x.nyp());
// calculating the stencil
SW=A(0,0)*(1.0/(hx*hy));
S= A(1,0)*(1.0/(hy*hy));
SE=A(2,0)*(1.0/(hx*hy));
W=A(0,1)*(1.0/(hx*hx));
E=A(2,1)*(1.0/(hx*hx));
NW=A(0,2)*(1.0/(hx*hy));
N=A(1,2)*(1.0/(hy*hy));
NE=A(2,2)*(1.0/(hx*hy));
C=(A(1,1)*((1.0/(hx*hx))+(1.0/(hy*hy))));


er.evenquickerset(0); //sets the er for all nodes to zero
//#pragma omp parallel for
        for (int j=1; j<x.nyp()-1;j++)// y direction iteration
        {
			for (int i =1 ;i<x.nxp()-1;i++){//iterate x direction
				double newparts=0.0;// repository for the Ax
				newparts=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1)+C*x(i,j); //calculate based on old estimates
				er(i,j)=(b(i,j)-(newparts)); //adds the second degree of error term for each point
            }
        }
return er;}



// calculates error norm for this particular problem, when given x as estimate
double Errornorm (const Grid& x){
double errornorm=0,pi=3.141592653589793238462643383279502884;

//#pragma omp parallel for reduction (+:errornorm)
for (int j=0;j<x.nyp();j++){
    for (int i=0;i<x.nxp();i++){
        double error=0.0;
        error = cos(pi*i*x.hx())*cos(pi*j*x.hy()) - x(i,j);
        errornorm+=error*error;
    }
}
return sqrt(errornorm/(x.nyp()*x.nxp()));
}


// calculates the residual error, based on x,b and A
Grid Stencil_Grid_Mult (const Grid& x,const Grid& A,double hx,double hy){
double C,W,E,N,S,NW,NE,SW,SE;
Grid result(x.nxp(),x.nyp());
// calculating the stencil
SW=A(0,0)*(1.0/(hx*hy));
S= A(1,0)*(1.0/(hy*hy));
SE=A(2,0)*(1.0/(hx*hy));
W=A(0,1)*(1.0/(hx*hx));
E=A(2,1)*(1.0/(hx*hx));
NW=A(0,2)*(1.0/(hx*hy));
N=A(1,2)*(1.0/(hy*hy));
NE=A(2,2)*(1.0/(hx*hy));
C=(A(1,1)*((1.0/(hx*hx))+(1.0/(hy*hy))));



result.evenquickerset(0); //sets the er for all nodes to zero
//#pragma omp parallel for
        for (int j=1; j<x.nyp()-1;j++)// y direction iteration
        {
			for (int i =1 ;i<x.nxp()-1;i++){//iterate x direction
				double newparts=0.0;// repository for the Ax
				newparts=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1)+C*x(i,j); //calculate based on old estimates
				result(i,j)=newparts;
            }
        }
return result;}


double Grid_dot_product (const Grid& x,const Grid& b,char A='y'){

double result=0;
if (A =='T' || A =='t'){
        //cout<<"using the transpose"<<endl;
//#pragma omp parallel for reduction (+:result)
        for (int j=0; j<x.nyp();j++)// y direction iteration
        {
			for (int i =0 ;i<x.nxp();i++){//iterate x direction
				result+=x(j,i)*b(i,j);
            }
        }}
else {
    //#pragma omp parallel for reduction (+:result)
        for (int j=1; j<x.nyp()-1;j++)// y direction iteration
        {
			for (int i =1 ;i<x.nxp()-1;i++){//iterate x direction
				result+=x(i,j)*b(i,j);
            }
        }

}
return result;}


//Grid L2 norm calculator
double L2norm (const Grid &x){
double norm=0;
//#pragma omp parallel for reduction(+:norm)
for (int j=0; j<x.nyp();j++)// y direction iteration
    {
        for (int i =0 ;i<x.nxp();i++){//iterate x direction
            norm+=x(i,j)*x(i,j);//r square sum
            }
    }
//returns square root of the sum above, divided by the number of grid points.
return sqrt(norm/(x.nxp()*x.nyp()));
}
















