#ifndef MESHDATA_H
#define MESHDATA_H
#include  "vars.h"
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
// TODO MAKE CHANGE TO MESH CLASS 
// PLACE INIT AND MEMORY ALLOCATION METHODS INSIDE CLASS 
// TODO Check prims are all correct
// 
// IMPORTANT!!!
// MESH ONLY STORES RIGHT FACES VALUES, SINCE LEFT VALUES ARE CALCULATED 
// BY PRECEEDING CELL
class  meshcell
{
    private:
        // Unfortunately these have to be fixed for arrays
        
        double prim[4]; // Primative variables 
        double cons[4]; // Conserved variables 
        double pos[2]; // Position of the center of the mesh cell
        double J[8]; // Jacobian Matrix/Gradients
        
    public:
        double F[16]; // Flux function
        double faces[16]; // Left and Right face extraps of prims
        // The coordinates of the cell in the 2d array
        int xcoord;
        int ycoord;
        meshcell()
        {

        }

        void setPos(double x , double y=0., double z=0)
        {
            // TODO assertions for setting array pos
            //std::cout<<"X and Y in: "<< x << y<<std::endl;
            if (y == 0){
                // 1D 
                pos[0] = x;
            }
            else if (z==0){ 
                // 2D
                 pos[0] = x;
                 pos[1] = y;
            }
            else{
                // 3D
                pos[0] = x;
                pos[1] = y;
                pos[2] = z;
            }                
           
            // std::cout<< "x pos: "<< pos[0] << " y pos: " << pos[1] << std::endl;
        }

        std::string printCell()
        {
            /*std::string celldata;
            celldata = std::to_string(pos[0]) + " "  + std::to_string(pos[1]) + " ";
            celldata += std::to_string(prim[1]) + " ";
            celldata += std::to_string(prim[2]) + " ";
            celldata += std::to_string(prim[0]) + " ";
            celldata += std::to_string(prim[3]);
            return celldata;
            */

            // New version that allows full precision
            std::stringstream celldata;
            celldata << std::setprecision(17) << pos[0] << " ";
            celldata << std::setprecision(17) << pos[1] << " ";
            celldata << std::setprecision(17) << prim[1] << " ";
            celldata << std::setprecision(17) << prim[2] << " ";
            celldata << std::setprecision(17) << prim[3] << " ";
            celldata << std::setprecision(17) << prim[0];
            return celldata.str();
        }

        void printCons()
        {
            for (int i=0; i<4; i++){
                std::cout << cons[i] << " ";
            }
            std::cout << std::endl;
        }
        
        void printPrims()
        {
            for (int i=0; i<4; i++){
                std::cout << prim[i] << " ";
            }
            std::cout << std::endl;
        }
        void printGrad()
        {
            for (int i=0; i<8; i++){
                std::cout << J[i] << " ";
            }
            std::cout << std::endl;
        }

        void printFaces()
        {
            for (int i=0; i<16; i++){
                std::cout << faces[i] << " ";
            }
            std::cout << std::endl;
        }
       
        void printFlux()
        {
            for (int i=0; i<16; i++){
                std::cout << F[i] << " ";
            }
            std::cout << std::endl;
        }
        void setPrim(double vx, double vy, double rho, double P)
        {
            prim[0] = rho;
            prim[1] = vx;
            prim[2] = vy;
            prim[3]= P;

        }

        void printPos()
        {
            std::cout<< "x pos: "<< pos[0] << " y pos: " << pos[1] << std::endl; 
        }

        double getSoundspeed()
        {
            //std::cout << "Gamma: " << gam << " P: " << prim[3] << " rho: " << prim[0] << std::endl;
            double cs = sqrt(gam*prim[3]/prim[0]);
            //std::cout << cs << std::endl;
            return cs;
        }

        double getVelocity()
        {
            double v = sqrt(prim[1]*prim[1] + prim[2]*prim[2]);
            //std::cout << v << std::endl;
            return v;
        }

        void cons2prim()
        {
            double momx = cons[1];
            double momy = cons[2];
            double energy = cons[3];
            
            double rho = cons[0]/vol;
            double vx;
            double vy;
            // vx
            prim[1] = momx / rho / vol;
            // vy
            prim[2] = momy / rho / vol;
            // rho
            prim[0] = rho;
            vx = prim[1];
            vy = prim[2];
            // P 
            prim[3] = (energy/vol - 0.5*rho * (vx*vx + vy*vy))*(gam-1);
            //std::cout << "P: " << prim[3] << std::endl;    
        }

        void prim2cons()
        {
            // Maybe better to define variables first i.e rho
            // less code changes that way
            double rho;
            double vx;
            double vy;
            double P;
            rho = prim[0];
            vx = prim[1];
            vy = prim[2];
            P = prim[3];
            
            // Mass 
            cons[0] = rho * vol;
            // Mom x
            cons[1] = rho*vx*vol;
            // Mom y
            cons[2] = rho*vy*vol;
            // Energy
            cons[3] = (P/(gam-1.) + 0.5*rho*(vx*vx + vy*vy))*vol;
        }

        void set_cons()
        {

        }

        // Should this be in a new class that inherits from meshcell?
        // Stick get gradients, FL and FR into class 
        // That way it is unique to solving method i.e KT/Rusanov
        void setGradients(double rho_x, double rho_y,double vx_x,double vx_y,double vy_x,double vy_y,double P_x, double P_y)
        {
            J[0] = rho_x;
            J[1] = rho_y;
            J[2] = vx_x;
            J[3] = vx_y;
            J[4] = vy_x;
            J[5] = vy_y;
            J[6] = P_x;
            J[7] = P_y;

        }

        void getPrim(double &rho, double &vx, double &vy, double &P)
        {
            rho = prim[0];
            vx = prim[1];
            vy = prim[2];
            P = prim[3];

        }

        void extraptime(double dt)
        {
            double rhoold;
            double vxold;
            double vyold;
            double Pold;

            rhoold = prim[0];
            vxold = prim[1];
            vyold = prim[2];
            Pold = prim[3];

            prim[0] = rhoold - 0.5*dt *(vxold*J[0] + rhoold*J[2] + vyold*J[1] + rhoold*J[5]);
            prim[1] = vxold - 0.5*dt *(vxold*J[2] + vyold*J[3] + (1./rhoold)*J[6]);
            prim[2] = vyold - 0.5*dt *(vxold*J[4] + vyold*J[5] + (1./rhoold)*J[7]);
            prim[3] = Pold - 0.5*dt *(gam*Pold *(J[2] + J[5]) + vxold*J[6] + vyold*J[7]); 
        }

        void extrapface()
        {
            /*
            double Jl[8];
            Jl[0] = rho_xl;
            Jl[1] = rho_yl;
            Jl[2] = vx_xl;
            Jl[3] = vx_yl;
            Jl[4] = vy_xl;
            Jl[5] = vy_yl;
            Jl[6] = P_xl;
            Jl[7] = P_yl;
            */
            // TODO COMPLETE RE-WRITE 
            // Must take values from other mesh cells 
            // Either a public method which grabs prims and 
            // sets faces manually or calculate L values and shift array
            //for (int i=0; i<4; i++){
                // Should clean this up later
                // These are backward!!!
                // R.H.S face (L.H.S of interface)
           //     faces[i*4+1] = prim[i] + J[2*i]*dx/2.;
                // L.H.S face (R.H.S of interface)
           //     faces[i*4] = prim[i] - J[2*i]*dx/2.;
                // DOWN face (UP face of interface)
           //     faces[i*4+3] = prim[i] + J[2*i+1]*dx/2.;
                // UP face (Down face of interfaces)
           //     faces[i*4+2] = prim[i] - J[2*i+1] *dx/2;
           // }
           // Rho values at faces 
           faces[0] = prim[0] - J[0]*dx/2.;
           faces[1] = prim[0] + J[0]*dx/2.;
           faces[2] = prim[0] - J[1]*dx/2.;
           faces[3] = prim[0] + J[1]*dx/2.;
           faces[4] = prim[1] - J[2]*dx/2.;
           faces[5] = prim[1] + J[2]*dx/2.;
           faces[6] = prim[1] - J[3]*dx/2.;
           faces[7] = prim[1] + J[3]*dx/2.;
           faces[8] = prim[2] - J[4]*dx/2.;
           faces[9] = prim[2] + J[4]*dx/2.;
           faces[10] = prim[2] - J[5]*dx/2.;
           faces[11] = prim[2] + J[5]*dx/2.;
           faces[12] = prim[3] - J[6]*dx/2.;
           faces[13] = prim[3] + J[6]*dx/2.;
           faces[14] = prim[3] - J[7]*dx/2.;
           faces[15] = prim[3] + J[7]*dx/2.;
        
        }

        void shiftLeft(){

        }
        void calcflux()
        {
            // Faces guide
            // 0-3 rho_lx, rho_rx, rho_ly,rho_ry
            // 4-7 vx_lx, vx_rx, vx_ly, vx_ry
            // 8-11 vy_lx, vy_rx, vy_ly,vy_ry
            // 12-15 P_lx, P_rx, P_ly, P_ry
            
            double eng_xl = faces[12]/(gam-1) + 0.5*faces[0] * (faces[4]*faces[4] + faces[8]*faces[8]);
            double eng_xr = faces[13]/(gam-1) + 0.5*faces[1] * (faces[5]*faces[5] + faces[9]*faces[9]);
            double eng_yl = faces[14]/(gam-1) + 0.5*faces[2] * (faces[6]*faces[6] + faces[10]*faces[10]);
            double eng_yr = faces[15]/(gam-1) + 0.5*faces[3] * (faces[7]*faces[7] + faces[11]*faces[11]);


            double rho_starx = 0.5*(faces[0] + faces[1]);
            double rho_stary = 0.5*(faces[2] + faces[3]);

            double momx_starx = 0.5*(faces[0]*faces[4] + faces[1]*faces[5]);
            double momx_stary = 0.5*(faces[2]*faces[6] + faces[3]*faces[7]);
            
            double momy_starx = 0.5*(faces[0]*faces[8] + faces[1]*faces[9]);
            double momy_stary = 0.5*(faces[2]*faces[10] + faces[3]*faces[11]);
            

            double eng_starx = 0.5*(eng_xl + eng_xr);
            double eng_stary = 0.5*(eng_yl + eng_yr);

            double P_starx = (gam-1)*(eng_starx -0.5*(momx_starx*momx_starx +momy_starx + momy_starx)/rho_starx);
            double P_stary = (gam-1)*(eng_stary -0.5*(momx_stary*momx_stary +momy_stary + momy_stary)/rho_stary);
            
            F[0] = momx_starx;
            F[1] = momx_stary;
            F[2] = (momx_starx*momx_starx)/rho_starx + P_starx;
            F[3] = (momx_stary*momx_stary)/rho_stary + P_stary;
            F[4] = momx_starx * momy_starx/rho_starx;
            F[5] = momx_stary * momy_stary/rho_stary;
            F[6] = (eng_starx + P_starx) * momx_starx/rho_starx;
            F[7] = (eng_stary + P_stary) * momx_stary/rho_stary;

            // Find max signal speed
            double C_lx = sqrt(gam*faces[12]/faces[0]) + abs(faces[4]);
            double C_rx = sqrt(gam*faces[13]/faces[1]) + abs(faces[5]);
            double C_ly = sqrt(gam*faces[14]/faces[2]) + abs(faces[6]);
            double C_ry = sqrt(gam*faces[15]/faces[3]) + abs(faces[7]);

            double C_x = fmax(C_lx,C_rx);
            double C_y = fmax(C_ly,C_ry);

            // add diffusive term
            F[0] -= C_x * 0.5*(faces[0] - faces[1]);
            F[1] -= C_y * 0.5*(faces[2] - faces[3]);
            F[2] -= C_x * 0.5*(faces[0]*faces[4] - faces[1]*faces[5]);
            F[3] -= C_y * 0.5*(faces[2]*faces[6] -  faces[3]*faces[7]);
            F[4] -= C_x * 0.5*(faces[0]*faces[8] - faces[1]*faces[9]);
            F[5] -= C_y * 0.5*(faces[2]*faces[10] - faces[3]*faces[11]);
            F[6] -= C_x * 0.5*(eng_xl - eng_xr);
            F[7] -= C_y * 0.5*(eng_yl - eng_yr);

        }

        void setCons(double mass, double momx, double momy, double energy)
        {
            cons[0] = mass;
            cons[1] = momx;
            cons[2] = momy;
            cons[3] = energy;
        }

        void addFlux(double dt)
        {
            std::cout << "Momy fluxes: " << F[8] << " " << F[9] << " " << F[10]<< " " << F[11] << std::endl;
            // take sum of neighbours 
            double massold = cons[0];
            double momxold = cons[1];
            double momyold = cons[2];
            double energyold = cons[3];

            cons[0]  = massold -dt*dx*F[1]+dt*dx*F[0] - dt*dx*F[3] + dt*dx*F[2];
            cons[1]  = momxold -dt*dx*F[5]+dt*dx*F[4] - dt*dx*F[7] + dt*dx*F[6];
            cons[2]  = momyold -dt*dx*F[9] + dt*dx*F[8] - dt*dx*F[11] + dt*dx*F[10];
            cons[3]  = energyold -dt*dx*F[13]+ dt*dx*F[12] - dt*dx*F[15] + dt*dx*F[14];

            std::cout << "Momy: " << cons[2] << std::endl;
        }
};
#endif 
