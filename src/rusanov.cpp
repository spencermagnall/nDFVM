#include "flux.h"

void getLeftFace(meshcell* mesh, int currentcell){

    // Get x and y coords of current cell 
    int xcoord = mesh[currentcell].xcoord;
    int ycoord = mesh[currentcell].ycoord;

    // Get the cell bellow and to the right of the current cell 
    int rightcell;
    int bellowcell;
    
    // right cell
    if (xcoord + 1 > N-1) {
        xcoord = 0;
    } else {
        xcoord = xcoord + 1;
    }
    rightcell = (ycoord*N) + xcoord;
    
    mesh[currentcell].faces[0] = mesh[rightcell].faces[0];
    mesh[currentcell].faces[4] = mesh[rightcell].faces[4];
    mesh[currentcell].faces[8] = mesh[rightcell].faces[8];
    mesh[currentcell].faces[12] = mesh[rightcell].faces[12];

    // cell bellow 
    xcoord = mesh[currentcell].xcoord;
    if (ycoord + 1 > N-1){
        ycoord = 0;
    } else {
        ycoord += 1;
    }
    
    bellowcell = (ycoord*N) + xcoord;

    mesh[currentcell].faces[2] = mesh[rightcell].faces[2];
    mesh[currentcell].faces[6] = mesh[rightcell].faces[6];
    mesh[currentcell].faces[10] = mesh[rightcell].faces[10];
    mesh[currentcell].faces[14] = mesh[rightcell].faces[14];



}
void getLeftFaces(meshcell* mesh){
    double tempFaces[8];
    // Store the first values 
    int k = 0;
    for (int j=0; j < 8; j++){
        tempFaces[j] = mesh[0].faces[k];
        k+=2;
    }
    // Get left face value from cell on right
    for (int i=0; i < N*N - 1; i++){
        mesh[i].faces[0] = mesh[i+1].faces[0];
        mesh[i].faces[2] = mesh[i+1].faces[2];
        mesh[i].faces[4] = mesh[i+1].faces[4];
        mesh[i].faces[6] = mesh[i+1].faces[6];
        mesh[i].faces[8] = mesh[i+1].faces[8];
        mesh[i].faces[10] = mesh[i+1].faces[10];
        mesh[i].faces[12] = mesh[i+1].faces[12];
        mesh[i].faces[14] = mesh[i+1].faces[14];

    }

    // Set last value 
    /*mesh[N*N-1].faces[0] = tempFaces[0];
    mesh[N*N-1].faces[2] = tempFaces[1];
    mesh[N*N-1].faces[4] = tempFaces[2];
    mesh[N*N-1].faces[6] = tempFaces[3];
    mesh[N*N-1].faces[8] = tempFaces[4];
    mesh[N*N-1].faces[10] = tempFaces[5];
    mesh[N*N-1].faces[12] = tempFaces[6];
    mesh[N*N-1].faces[14] = tempFaces[7];
    */


}
void rusanov(double (&prims)[8], double  (&flux)[4] ){
   // prims guide
   // rho_l, rho_r
   // vx_l, vx_r
   // vy_l, vy_r
   // P_l, P_r
   // 8-15 the same but for the rhs boundadry 
    double engl = prims[6]/(gam-1) + 0.5*prims[0]*(prims[2]*prims[2] + prims[4]*prims[4]);
    double engr = prims[7]/(gam-1) + 0.5*prims[1]*(prims[3]*prims[3] + prims[5]*prims[5]);
    
    /*std::cout << "P_l " << prims[6] << std::endl;
    std::cout << "P_r " << prims[7] << std::endl;
    std::cout << "rho l " << prims[0] << std::endl;
    std::cout << "rho r " << prims[1] << std::endl;
    std::cout << "vx l " << prims[2] << std::endl;
    std::cout << "vx r " << prims[3] << std::endl;
    std::cout << "vy l " << prims[4] << std::endl;
    std::cout << "vy r " << prims[5] << std::endl;
    std::cout << "eng l " << engl << std::endl;
    std::cout << "eng r " << engr << std::endl;*/
    double rho_star = 0.5*(prims[0] + prims[1]);
    double momx_star = 0.5*(prims[0]*prims[2] + prims[1]*prims[3]);
    double momy_star = 0.5*(prims[0]*prims[4] + prims[1]*prims[5]);
    double eng_star = 0.5*(engl + engr);
    
    double P_star = (gam-1)*(eng_star-0.5*(momx_star*momx_star + momy_star*momy_star)/rho_star);
    //std::cout << "Pstar: " << P_star > 0 << std::endl;
    // Fluxes
    flux[0] = momx_star;
    flux[1] = (momx_star*momx_star)/rho_star + P_star;
    flux[2] = momx_star * momy_star/rho_star;
    flux[3] = (eng_star+P_star)*momx_star/rho_star;
    /*
    std::cout << "Mass: " << flux[0] << std::endl;
    std::cout << "momx: " << flux[1] << std::endl;
    std::cout << "momy: " << flux[2] << std::endl;
    std::cout << "Energy: " << flux[3] << std::endl;
    */
    //double englr = prims[14]/(gam-1) + 0.5*prims[8]*(prims[10]*prims[10] + prims[12]*prims[12]);
    //double engrr = prims[15]/(gam-1) + 0.5*prims[9]*(prims[11]*prims[11] + prims[13]*prims[13]);
    
    //rho_star = 0.5*(prims[8] + prims[9]);
    //momx_star = 0.5*(prims[8]*prims[10] + prims[9]*prims[11]);
    //momy_star = 0.5*(prims[8]*prims[12] + prims[9]*prims[13]);
    //eng_star = 0.5*(englr + engrr);
    
    //P_star = (gam-1)*(eng_star-0.5*(momx_star*momx_star + momy_star*momy_star)/rho_star);

    // Fluxes
    //flux[4] = momx_star;
    //flux[5] = (momx_star*momx_star)/rho_star + P_star;
    //flux[6] = momx_star * momy_star/rho_star;
    //flux[7] = (eng_star+P_star)*momx_star/rho_star;

    double C_L = sqrt(gam*prims[6]/prims[0]) + abs(prims[2]);
    double C_R = sqrt(gam*prims[7]/prims[1]) + abs(prims[3]);
    double C = fmax(C_L,C_R);
    
    //std::cout << "C: " << C << std::endl;
    double expres = -(C * 0.5 * (prims[0]*prims[4] - prims[1]*prims[5]));
    //std::cout << "momy exp: " << expres << std::endl; 
    // Add dissipation 
    flux[0] -= C*0.5*(prims[0]-prims[1]);
    flux[1] -= C*0.5*(prims[0]*prims[2]-prims[1]*prims[3]);
    flux[2] -= C*0.5*(prims[0]*prims[4]-prims[1]*prims[5]);
    flux[3] -= C*0.5*(engl-engr);
    
    //std::cout << "Mass: " << flux[0] << std::endl;
    //std::cout << "momx: " << flux[1] << std::endl;
    //std::cout << "momy: " << flux[2] << std::endl;
    //std::cout << "Energy: " << flux[3] << std::endl;
}
void extrap_to_time(meshcell* mesh, double dt)
{
    for (int cellno = 0; cellno < N*N; cellno++)
    {
        mesh[cellno].extraptime(dt);
    }
}

void extrap_to_face(meshcell* mesh)
{
    for (int cellno = 0; cellno < N*N; cellno++)
    {
        mesh[cellno].extrapface();
    }

    // Shift faces to the left after extrap
    // Since ul_j+1/2 = u_j+1/2 - dx/2 grad(u_j+1)
    //shiftFaceLeft(mesh);
    //
    //getLeftFaces(mesh);
    // Rather than shift we just do this in the getFlux section 
    // Shifting is a pain in the ass 
}


void getFlux(meshcell* mesh, double dt)
{
    // coords for the cells adjacent (periodic)
    int xcoord;
    int ycoord;
    int xcoordadj;
    int ycoordadj;
    int adjcellno;
    double rhocell;
    double rhoadj;
    double vxcell;
    double vxadj;
    double vycell;
    double vyadj;
    double Pcell;
    double Padj;
    double prims[8];
    double flux[4];
    double fluxes[16];

    for (int cellno = 0; cellno < N*N; cellno++)
    {
        xcoord = mesh[cellno].xcoord;
        ycoord = mesh[cellno].ycoord;
        // start with left boundary
        // get flux at right face of left adjacent cell 
        xcoordadj = (xcoord == 0) ? N-1 : xcoord - 1; 
        ycoordadj = ycoord;
        adjcellno = (ycoordadj*N) + xcoordadj;
        /*
        std::cout << "Xcoord: " << xcoord << std::endl;
        std::cout << "Ycoord: " << ycoord << std::endl;
        std::cout << "Xcoord adj: "<< xcoordadj << std::endl;
        std::cout << "Ycoord adj: " << ycoordadj << std::endl;
        std::cout << "Adj cellno: " << adjcellno << std::endl;
        */
        // flux calculation
        // Do x together 
        // Do y together
        //  LEFT IS LEFT, RIGHT IS RIGHT!
        // cell vars (left face of cell) left in calc
        // this sequential layout is wrong 
        prims[0] = mesh[cellno].faces[0];
        prims[2] = mesh[cellno].faces[4];
        prims[4] = mesh[cellno].faces[8];
        prims[6] = mesh[cellno].faces[12];
        // adj cell vars (right face) right in calc
        prims[1] = mesh[adjcellno].faces[1];
        prims[3] = mesh[adjcellno].faces[5];
        prims[5] = mesh[adjcellno].faces[9];
        prims[7] = mesh[adjcellno].faces[13];
        
        // calc flux for bounadry
        //if (cellno == 0){
            //std::cout << "True" << std::endl;
           rusanov(prims,flux);
        //}
        //exit(0);      
        fluxes[0] = flux[0];
        fluxes[4] = flux[1];
        fluxes[8] = flux[2];
        fluxes[12] = flux[3];
        
        
        // get coords of right adjcell
        xcoordadj = (xcoord == N-1) ?  0 : xcoord + 1;
        ycoordadj = ycoord;
        adjcellno = (ycoordadj*N) + xcoordadj;
        // cell vars (right face of cell) right in calc
        prims[1] = mesh[cellno].faces[1];
        prims[3] = mesh[cellno].faces[5];
        prims[5] = mesh[cellno].faces[9];
        prims[7] = mesh[cellno].faces[13];
        
        // adj cell vars (left face) left in calc 
        prims[0] = mesh[adjcellno].faces[0];
        prims[2] = mesh[adjcellno].faces[4];
        prims[4] = mesh[adjcellno].faces[8];
        prims[6] = mesh[adjcellno].faces[12];
       
        //if (cellno == 0)
            //std::cout << "True" << std::endl;
            //std::cout << "Xcoord is: " << xcoordadj << std::endl;
            rusanov(prims,flux);
        
        fluxes[1] = flux[0];
        fluxes[5] = flux[1];
        fluxes[9] = flux[2];
        fluxes[13] = flux[3];
        
        // SWAP X AND Y FOR VERTICAL VALUES   
        // Get coords for up adjacent face
        xcoordadj = xcoord;
        ycoordadj = (ycoord == 0 ) ? N-1 : ycoord - 1;
        adjcellno = (ycoordadj*N) + xcoordadj;
        // cell vars (up face (left)) left in calc
        prims[0] = mesh[cellno].faces[2];
        prims[4] = mesh[cellno].faces[6];
        prims[2] = mesh[cellno].faces[10];
        prims[6] = mesh[cellno].faces[14];
        
        // adj vars (down face (right)) right in calc
        prims[1] = mesh[adjcellno].faces[3];
        prims[5] = mesh[adjcellno].faces[7];
        prims[3] = mesh[adjcellno].faces[11];
        prims[7] = mesh[adjcellno].faces[15];
        
        //if (cellno == 0)
            rusanov(prims,flux);
        
        // AGAIN X AND Y ARE SWAPED
        fluxes[2] = flux[0];
        fluxes[6] = flux[2];
        fluxes[10] = flux[1];
        fluxes[14] = flux[3];


        // Get coords for down adjacent face
        xcoordadj = xcoord;
        ycoordadj = (ycoord == N-1 ) ? 0 : ycoord + 1;
        adjcellno = (ycoordadj*N) + xcoordadj;
        // cell vars (down face (right)) left in calc
        prims[1] = mesh[cellno].faces[3];
        prims[5] = mesh[cellno].faces[7];
        prims[3] = mesh[cellno].faces[11];
        prims[7] = mesh[cellno].faces[15];
        
        // adj vars (up face (left)) right in calc
        prims[0] = mesh[adjcellno].faces[2];
        prims[4] = mesh[adjcellno].faces[6];
        prims[2] = mesh[adjcellno].faces[10];
        prims[6] = mesh[adjcellno].faces[14];
        
        //if (cellno == 0)
            rusanov(prims,flux);
        fluxes[3] = flux[0];
        fluxes[7] = flux[2];
        fluxes[11] = flux[1];
        fluxes[15] = flux[3];
        
        for (int i=0; i<16; i++){
            mesh[cellno].F[i] = fluxes[i];
        }
        // apply to conserved variables 
        mesh[cellno].addFlux(dt);
        //exit(0);
    }

}

//void addFlux(meshcell* mesh)
//{
//    for (int cellno = 0; cellno < N*N; cellno++)
//    {
//        int xcoord = mesh[cellno].xcoord;
//        int ycoord = mesh[cellno].xcoord;
        //double mass = mesh[cellno].;
//        double momx;
//        double momy;
//        double energy;

        
  //  }
//}

//void rusanov(meshcell* mesh, double dt)
//{
// Extrapolate half-step in time (predictor)
//         extrap_to_time(mesh,dt);
         // Extrapolate face (corrector)
//         extrap_to_face(mesh);
         // Compute flux
//         getFlux(mesh);
         // Add flux to cons
//         addFlux(mesh);
//}
