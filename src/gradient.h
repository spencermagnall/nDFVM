#ifndef GRADIENT_H
#define GRADIENT_H
#include "meshdata.h"

void getCellcoord(int cellno,int &xcoord, int &ycoord)
{
    // NOT ENOUGH INFORMATION
    // x and y coords
    int i;
    int j;

    i = (cellno - j)/N;
    j = cellno - i*N;

    xcoord = i;
    ycoord = j;


}

// c++ mod is implemntation specific, roll our own
int modulus(int coord){
    
    int mod;
    if (coord < 0){
        mod = N + (coord);
    } else {
        mod =  coord % N;
    }

    return mod;
} 
    

void getGradients(meshcell* mesh)
{
    double rho_x;
    double rho_y;
    double rho_l;
    double rho_r;
    double vx_x;
    double vx_y;
    double vx_l;
    double vx_r;
    double vy_x;
    double vy_y;
    double vy_l;
    double vy_r;
    double P_x;
    double P_y;
    double P_l;
    double P_r;
    int leftcell;
    int rightcell;
    int leftcellx;
    int leftcelly;
    int rightcellx;
    int rightcelly;

    for (int cellno = 0; cellno < N*N; cellno++)
    {
        // get position of variables in array
        int xcoord = mesh[cellno].xcoord;
        int ycoord = mesh[cellno].ycoord;
        // have cellcoord as a class method for meshdata?
        //getCellcoord(cellno,xcoord,ycoord);
    
        // dx calc
        // Get left and right cellno's
        //  NEED to account for periodic boundaries
        
        if (xcoord == 0) {
            leftcellx = N-1;
        } else{
           leftcellx =  xcoord-1;
        }

        leftcell = ycoord*N + leftcellx;
        
        if (xcoord == N-1) {
            rightcellx = 0;
        } else{
           rightcellx = xcoord + 1; 
        } 
        
        rightcell = ycoord*N + rightcellx;
        // i-1
        mesh[leftcell].getPrim(rho_l,vx_l,vy_l,P_l);
        mesh[rightcell].getPrim(rho_r,vx_r,vy_r,P_r);
        
        rho_x = (rho_r - rho_l)/(2.*dx);
        vx_x  = (vx_r - vx_l)/(2.*dx);
        vy_x = (vy_r - vy_l)/(2.*dx);
        P_x = (P_r - P_l)/(2.*dx);

        // dy calc 
        if (ycoord == 0) {
            leftcelly = N-1;
        } else{
           leftcelly =  ycoord-1;
        }

        leftcell = leftcelly*N + xcoord;
        
        if (ycoord == N-1) {
            rightcelly = 0;
        } else{
           rightcelly = ycoord + 1; 
        }
        rightcell = rightcelly*N + xcoord;

        mesh[leftcell].getPrim(rho_l,vx_l,vy_l,P_l);
        mesh[rightcell].getPrim(rho_r,vx_r,vy_r,P_r);
        
        rho_y = (rho_r - rho_l)/(2.*dx);
        vx_y = (vx_r - vx_l)/(2.*dx);
        vy_y = (vy_r - vy_l)/(2.*dx);
        P_y = (P_r - P_l)/(2.*dx);

        mesh[cellno].setGradients(rho_x,rho_y,vx_x,vx_y,vy_x,vy_y,P_x,P_y);


    }
} 
#endif 
