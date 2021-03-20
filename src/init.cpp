#include <cmath>
#include <iostream>
#include "init.h"
//#include "vars.h"

void yexpress(){
    double yexpres = 0.;
    double X = 0.5*dx;
    double Y = X;
    double sigma = 0.05;
    yexpres = -(pow((Y-0.25),2.)/(2.* pow(sigma,2.)));
    yexpres +=  -(pow((Y-0.75),2.)/(2.*sigma*sigma));
    //std::cout << "Y express: " << yexpres << std::endl;

}
void init_mesh(meshcell*& mesh){
    int cellno = 0;
    double currentx = 0.;
    double currenty = 0.;
    // SUPER UGLY!!!
    while (cellno < N*N){
        for (int y=0; y < N; y++)
        {
            for (int x=0; x < N; x++)
            {
                //std::cout << "Mesh pointer inside: " << &mesh << std::endl;
                // set x and y pos
                mesh[cellno].setPos((x+1)*dx-0.5*dx,(y+1)*dx-0.5*dx);
                //mesh[cellno].printPos();
                //std::cout<< "X Pos: " << x*dx*0.5 << "Y Pos: "<< y*dx*0.5 << std::endl;
                cellno++;
            }
        }
       
    }
}

void init_fluid(meshcell *& mesh){
    int cellno = 0;
    double vx; 
    double vy;
    double P;
    double rho;
    double sigma = 0.05/sqrt(2.);
    double  w0 = 0.1;
    double X;
    double Y;
    double yexpres1;
    double yexpres2;

    // SUPER UGLY!!!
       while (cellno < N*N){
           for (int y=0; y < N; y++)
           {
               for (int x=0; x < N; x++)
               {
                   //std::cout << "x + 1: " << x+1;
                   //std::cout << " y + 1: " << y+1 << std::endl;
                   // Setup fluid values
                   X = (dx*(x+1)-0.5*dx);
                   Y = (dx*(y+1)-0.5*dx);
                   // Assumes boxsize of 1
                   // Should be more general in future revision
                   // Hacky but simple
                   rho = 1. + (double)(abs(Y-0.5) < 0.25);
                   P = 2.5;
                   vx = -0.5 + (double)(abs(Y-0.5) < 0.25);
                   yexpres1 = -((Y-0.25)*(Y-0.25)/(2.*sigma*sigma));
                   yexpres2 = - ((Y-0.75)*(Y-0.75)/(2.*sigma*sigma));
                   //std::cout << yexpres << std::endl;
                   //yexpress();
                   //exit(0);
                   vy = w0*sin(4*M_PI*X)*(exp(yexpres1) + exp(yexpres2));
                   std::cout << vy << std::endl;
                   mesh[cellno].setPrim(vx,vy,rho,P);
                   mesh[cellno].xcoord = x;
                   mesh[cellno].ycoord = y;
                   //if (rho == 1.){
                   //std::cout<<"rho: "<< rho << std::endl; 
                   //}
                   if (cellno == 5598){
                       std::cout << "x: " << X << std::endl;
                       std::cout << "y: " << Y << std::endl;
                       std::cout << "xcoord: " << x << std::endl;
                       std::cout << "ycoord: " << y << std::endl;
                       std::cout << "rho: " << rho << std::endl;
                       std::cout << "vx:  " << vx << std::endl;
                       std::cout << "vy: " << vy << std::endl;
                       std::cout << "P: " << P << std::endl;
                       //exit(0);
                  }
                  
                   cellno++;
               }
           }
      }
}

