#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <mpi.h>

#include "Timer.h"
#include "Grid.hpp"

using namespace std;

// calculates optimal dimensions of the cart communicator, based on world and domain size
vector <int> Getdims(int worldsize,int nx=1, int ny=1){
int dim1 = sqrt(worldsize);
int dim2 = worldsize/dim1;
double ration=nx/ny;
//checking if dim1*dim2 equals to world size
if (worldsize== dim1*dim2){
    if (ration>= 1){
        return {dim2,dim1};
    }
    else{
        return {dim1,dim2};
    }
}
else {
    while (worldsize!=dim1*dim2 ){
        if (dim1*dim2<worldsize) dim2++;
        else dim1--;
        }
    }
    if (ration >= 1){
        return {dim2,dim1};
    }
    else{
        return {dim1,dim2};
    }
}

//MPI to update neighbours of a grid
int Updateneighbours(Grid &x,MPI_Comm Cart_comm, MPI_Datatype &columns,std::vector<int> neighbours) {
int dump,tag=1;

MPI_Request request[8];
MPI_Status status[8];
//update ghost zones
//up
dump= MPI_Issend( &x(0,x.nyp()-2), x.nxp(),
MPI_DOUBLE, neighbours[1],
tag,Cart_comm,
&request[0] );

dump= MPI_Irecv( &x(0,x.nyp()-1), x.nxp(),
MPI_DOUBLE, neighbours[1],
tag, Cart_comm,
&request[1]);

//down
dump= MPI_Issend( &x(0,1), x.nxp(),
MPI_DOUBLE, neighbours[0],
tag,Cart_comm,
&request[2] );

dump= MPI_Irecv( &x(0,0), x.nxp(),
MPI_DOUBLE, neighbours[0],
tag, Cart_comm,
&request[3]);

//left
dump= MPI_Issend( &x(1,0), 1,
columns, neighbours[2],
tag,Cart_comm,
&request[4] );

dump= MPI_Irecv( &x(0,0), 1,
columns, neighbours[2],
tag, Cart_comm,
&request[5]);

//right

dump= MPI_Issend( &x(x.nxp()-2,0), 1,
columns, neighbours[3],
tag,Cart_comm,
&request[6] );

dump= MPI_Irecv( &x(x.nxp()-1,0), 1,
columns, neighbours[3],
tag, Cart_comm,
&request[7]);


MPI_Waitall(8, request, status);
return dump;

}



//main program, CG MPI parallel
int main(int argc,char** argv) //can be used with command line variables, optionally, can also define threads
{

    // Definition of the variables
    int size;  //The total number of processes
    int rank;  //The rank/number of this process

    int rank2;  //The rank/number of the whole process WORLD_COMM


    // MPI initialization
	MPI_Init(&argc, &argv);

    // Determining the number of CPUs and the rank for each CPU
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank2);


//input pharsing
int nx,ny,iter,print=1;
if (argc==4) {
nx=stoi(argv[1]);
ny=stoi(argv[2]);
iter=stoi(argv[3]);
if (iter==0) iter=1000000000;
}
else if (argc==5){
nx=stoi(argv[1]);
ny=stoi(argv[2]);
iter=stoi(argv[3]);
if (iter==0) iter=1000000000;
print=stoi(argv[4]);
}
else { // if no command line arguments are used, the code does not crash or stop. it asks for them
    cout<< "You can choose the number of nodes in x and y direction, iteration, and print status respectively, by running the software with command parameters"<<endl;
    cout<< "Nx (0 to exit) = ";
    cin >>nx;
    if (cin.fail() || nx<3) {
          return -1;
          }
    cout << "\n Ny = ";
    cin >>ny;
    if (cin.fail() || ny<3) {
          return -1;
          }
    cout << "\n Number of iterations = ";
    cin >>iter;
    if (cin.fail()||iter<1) {
          cout << "That was not acceptable." << endl;
          return -1;
          }
    }

// calculate optimal cartesian dimensions
vector <int> dims=Getdims(size,nx,ny);
vector<int> periods={0,0};
//generate cartasian communicator
MPI_Comm Cart_comm;
MPI_Cart_create( MPI_COMM_WORLD, 2, &dims[0],&periods[0], 1, &Cart_comm);

if (rank2==0) cout<< "Cart comm created with size ("<<dims[0]<<","<<dims[1]<<")" <<endl;

// calculate node rank in cartcomm
   MPI_Comm_size(Cart_comm, &size);
    MPI_Comm_rank(Cart_comm, &rank);



//get neighbours
int D,R,L,U;
MPI_Cart_shift (Cart_comm,0,1,&L,&R);
MPI_Cart_shift (Cart_comm,1,1,&D,&U);
std::vector<int> neighbours={D,U,L,R};



//get node coordinates
vector<int> my_coords(2);
MPI_Cart_coords( Cart_comm, rank, 2, &my_coords[0] );

//calculating the mapping between local and global domain
int localnx,localny,globalIstart,globalJstart,globalIend,globalJend;

globalIstart=my_coords[0]*(nx-2)/dims[0];
globalIend=2+(my_coords[0]+1)*(nx-2)/dims[0];

globalJstart=my_coords[1]*(ny-2)/dims[1];
globalJend=2+(my_coords[1]+1)*(ny-2)/dims[1];

//calculating  the local grid size
localnx=globalIend-globalIstart;
localny=globalJend-globalJstart;

//calculating global hx and hy (cell size)
double hx=1.0/(nx-1);
double hy=1.0/(ny-1);

//defining the holders of stencil, results, right hand side and residual and stencil * residual (Ares)
Grid A(3,3); //define variables as Grid A for the stencil
Grid x(localnx,localny);
Grid b(localnx,localny);
Grid resold(localnx,localny),resnew(localnx,localny),Ares(localnx,localny);

//create a datatype for columns
MPI_Datatype columns;
int elements_in_column=x.nyp();
int stride=x.nxp();

MPI_Type_vector(elements_in_column,1,stride,MPI_DOUBLE,&columns);
MPI_Type_commit(&columns);

double resnormold,resnormnew,epsilon=0.0000000001,alpha=0,beta=0; // for error

A.quickset({{0,-1,0},{-1,2,-1},{0,-1,0}});// Coefficient matrix
x.evenquickerset(0);//x vector (initial estimate)
b.evenquickerset(0);//b vector (right hand side)


double pi=3.141592653589793238462643383279502884;



//define the boundary condition - basically g function
//in x direction with y=0
if (my_coords[1]==0) {
//pragma omp parallel for
for(int i=0;i<x.nxp();i++ ){
    x(i,0)=cos(pi*(i+globalIstart)*hx);
}
}
//in x direction with y=ymax
if (my_coords[1]==dims[1]-1){
//pragma omp parallel for
for(int i=0;i<x.nxp();i++ ){
    x(i,x.nyp()-1)=-cos(pi*(i+globalIstart)*hx);
}
}

//in y direction with x=0
if (my_coords[0]==0){
//pragma omp parallel for
for(int i=0;i<x.nyp();i++ ){
    x(0,i)=cos(pi*(i+globalJstart)*hy);
}
}

//in y direction with x=xmax
if (my_coords[0]==dims[0]-1){
//pragma omp parallel for
for(int i=0;i<x.nyp();i++ ){
    x(x.nxp()-1,i)=-cos(pi*(i+globalJstart)*hy);
}
}


//define the b for the internal point -- basically f function
//pragma omp parallel for
for(int j=0;j<x.nyp();j++ ){
        for(int i=0; i<x.nxp();i++){
            b(i,j)= 2*(pi*pi)*cos(pi*(i+globalIstart)*hx)*cos(pi*(j+globalJstart)*hy);
        }
}

//calculating residual
resold = Residual(x,b,A,hx,hy);
//updating ghost layers of res
Updateneighbours(resold,Cart_comm,columns,neighbours);
resnew=resold;
resnormold=Grid_dot_product(resnew,resold);


MPI_Allreduce(&resnormold, &resnormnew, 1, MPI_DOUBLE, MPI_SUM, Cart_comm);

if (rank==0)  cout<<"Problem initialized, initial residual norm is "<<sqrt(resnormnew/(nx*ny))<<endl;

resnormold= resnormnew;


// starting iterations here *****************************
if (rank==0) cout<< "Starting iteratoion " <<endl;
siwir::Timer counter;
int k=0;
while (resnormnew*hx*hy>epsilon*epsilon && k<iter){


if(rank==0) cout<< "*" ;

// A * d
Ares=Stencil_Grid_Mult(resnew,A,hx,hy);

// res.d
double temp =Grid_dot_product(resnew,Ares);
double temp3;
MPI_Allreduce(&temp,&temp3,1,MPI_DOUBLE,MPI_SUM,Cart_comm);

// alpha = resnorm/res.d
alpha = resnormnew/temp3;
//updating x
x=Pointwisesum(x,resnew,alpha);
//calculating new res norm
resold = Pointwisesum(resold,Ares,-1*alpha);
resnormold=resnormnew;
//calculating new resnorm
resnormnew=Grid_dot_product(resold,resold);
double temp2;
MPI_Allreduce(&resnormnew, &temp2, 1, MPI_DOUBLE, MPI_SUM, Cart_comm);
resnormnew=temp2;
//calculating beta
beta = resnormnew/resnormold;
//updating res
resnew=Pointwisesum(resold,resnew,beta);
Updateneighbours(resnew,Cart_comm,columns,neighbours);

if(k%100==0 && rank==0) cout<< " \nResnorm is " << sqrt(resnormnew/(nx*ny))<<" Iteration continues\n"  ;
k++;
}
//making sure all nodes are done
MPI_Barrier(Cart_comm);
auto wallclocktime=counter.elapsed();
//calculation done here ********

if (print==0) {
        if (rank ==0){
                if (rank==0) cout<< "\nUpdating report.txt\n" << endl;
        // To facilitate testing, outputs test results to file.
    ofstream report;
    report.open ("report.txt",std::ios_base::app);
        //report << setw(20)<< "run time" << setw(20)<< "Error norm"<< setw(20)<<"gridsize"<<setw(20)<<"number of cores"<<endl;
    report <<setw(20)<<setprecision(15)<<"nx: "<<nx<<setw(20)<<setprecision(15)<<"ny: "<<ny<<setw(20)<<setprecision(15)<<"Np: "<< size<<setw(20)<<setprecision(15)<< wallclocktime<< " s"<<setw(20)<<setprecision(15)<< "Iterations: "<< k  <<endl;
    report.close();}
}
else{
if (rank==0) cout<< "Time to complete " << k << " iterations was " <<wallclocktime<< " seconds" << endl;
if (rank==0) cout<<"Final residual norm is : "<<sqrt(resnormnew/(nx*ny))<<endl;

//initialize solution.txt
if (rank==0){
    printf("writing solution file\n");
    ofstream solutionfile;
    solutionfile.open ("solution.txt");
    solutionfile.close();
}
Updateneighbours(x,Cart_comm,columns,neighbours);

if (size==1) {
    Writesolution(x, hx,hy, globalIstart, globalJstart,dims,my_coords);
}
else {
//writing solution.txt
MPI_Barrier(Cart_comm);
for (int rankj=0;rankj<dims[1];rankj++){ //iterating through nodes in j direction
    if (rankj==my_coords[1]){ //getting all the nodes in x direction with the same coordinates in Y direction

            if (my_coords[0]==0){//the first node in x direction in the row (j=n) is responsible for gathering information and writing to file- called master
                    //printf("Node %d writing for itself\n",rank);
                    Grid solution(nx,x.nyp()); // a placeholder for all the variables in that j=n part of the domain from x=0 to x=nx
                    for (int j = 0;j<x.nyp();j++){
                        for(int i =0;i<x.nxp();i++){
                            solution(i,j)=x(i,j); //the master writes its own values to the placeholder for solution
                        }
                    } if (dims[0]>1){
                    //printf("Node %d writing for rankj %d\n",rank,rankj);
                    for (int ii=1;ii<dims[0];ii++){ //master iterates through neighbouring nodes in x direction and fetches the x values
                        //calculating the domain of the neighbour node
                        int neighbourglobalIstart=ii*(nx-2)/dims[0];
                        int neighbourglobalIend=2+(ii+1)*(nx-2)/dims[0];
                        MPI_Request req;
                        MPI_Status stat;
                        int neighbourrank;
                        std::vector<int> neighbourcoords;
                        neighbourcoords={ii,rankj};
                        MPI_Cart_rank(Cart_comm, &neighbourcoords[0], &neighbourrank); //getting the rank of the neighbour node
                        //calculating  the local grid size of neighbour node
                        int neighbournx=neighbourglobalIend-neighbourglobalIstart;
                        int number_of_neighbour_nodes=neighbournx*x.nyp();
                        Grid temp(neighbournx,x.nyp()); // temporary holder for the values in neighbour node
                        MPI_Irecv(&temp(0,0),number_of_neighbour_nodes , MPI_DOUBLE, neighbourrank, 0, Cart_comm, &req);
                        MPI_Wait(&req,&stat);
                        //printf("Uptading solution holder with values from node %d\n",neighbourrank);
                        int imax = neighbournx-1;
                        if (ii==dims[0]-1) imax=neighbournx;// setting boundary condition
                        //printf("Node %d recieved data from rankj %d\n",rank,neighbourrank);
                        for (int j = 0;j<x.nyp();j++){ //iterating through the values of the neighbour node and filing the solution placeholder
                             for(int i =1;i<imax;i++){
                                 solution(i+neighbourglobalIstart,j)=temp(i,j);
                                }
                            }
                    }

                        //after the solution placeholder is filled with all the values from neighbouring nodes
                       ofstream solutionfile;// writing the solution placeholder to solution.txt file
                        solutionfile.open ("solution.txt",std::ios_base::app);
                        int jstart=1,jmax=x.nyp()-1;
                        if (rankj==0) jstart=0;
                        if (rankj==dims[1]-1) jmax=x.nyp();


                        for (int j=jstart; j<jmax ;j++) { //printing results of x
                                if (j%2==0){
                                solutionfile << 0 <<" " <<(j+globalJstart)*hy << " "<< cos(pi*(j+globalJstart)*hy)<<"\n";
                                for (int i=1;i<nx-1;i++){
                                        solutionfile << i*hx <<" " <<(j+globalJstart)*hy << " "<<solution(i,j)<<"\n";
                                        }
                                solutionfile << 1 <<" " <<(j+globalJstart)*hy << " "<< -cos(pi*(j+globalJstart)*hy)<<"\n";
                                }
                                else{
                                    solutionfile << 1 <<" " <<(j+globalJstart)*hy << " "<< -cos(pi*(j+globalJstart)*hy)<<"\n";
                                    for (int i=nx-1;i>0;i--){
                                        solutionfile << i*hx <<" " <<(j+globalJstart)*hy << " "<<solution(i,j)<<"\n";
                                        }

                                    solutionfile << 0 <<" " <<(j+globalJstart)*hy << " "<< cos(pi*(j+globalJstart)*hy)<<"\n";

                                }
                                }
                        solutionfile.close();}
                        //printf("node %d,%d Done writing solution file\n",my_coords[0],my_coords[1]);

            }
            else {
                        //all the node except for the master (left most node) pass their x values to the master
                        MPI_Request req;
                        MPI_Status stat;
                        int neighbourrank;
                        std::vector<int> neighbourcoords;
                        neighbourcoords={0,rankj};
                        MPI_Cart_rank(Cart_comm, &neighbourcoords[0], &neighbourrank); //getting the rank of the master
                        //calculating  the local grid size
                        int number_of_nodes=x.nxp()*x.nyp();
                        MPI_Isend(&x(0,0),number_of_nodes , MPI_DOUBLE, neighbourrank, 0, Cart_comm, &req); // sending to master
                        //printf("Node %d sending informaiton to node %d\n",rank,neighbourrank);
                        MPI_Wait(&req,&stat);

            }
}
MPI_Barrier(Cart_comm); //all nodes wait for the j=n nodes to finish writing to the solution.txt file
}
}}


MPI_Type_free( &columns );
MPI_Comm_free(&Cart_comm);
MPI_Finalize();
return 0;
}


