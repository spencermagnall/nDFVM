#include <iostream>
#include <stdlib.h>
#include "vars.h"
#include "application.h"
#include "init.h"
#include "output.h"
#include "meshdata.h"
#include "cfl.h"
#include "gradient.h"
#include "flux.h"
void run(void)
{
   // probably move this later 
    double dt;
    double t = 0.;
    int snapshot = 1;
    double tprint = (double)snapshot*tout;
    
    // Create Mesh,  N*N cells 
    meshcell* Mesh = (meshcell*)malloc(sizeof(meshcell)* N*N);
    init_mesh(Mesh);
    //std::cout<< "Mesh pointer outside: " << &Mesh << std::endl;
    //Mesh[1].printPos();
    //exit(0);
    // SETUP INITIAL CONDITIONS
    init_fluid(Mesh);
    //Mesh[0].printPos();
    //output_mesh(Mesh,t,snapshot);
    //exit(0);
    for (int i=0; i < N*N; i++){
        Mesh[i].prim2cons();
    }
    //Mesh[0].printCons();
    //exit(0);
    // MAIN LOOP
    while (t < tmax){
      
        Mesh[5998].printPrims(); 
        // Get primatives 
        for (int i=0; i < N*N; i++)
        {
            Mesh[i].cons2prim();
        }
        Mesh[5998].printPrims(); 
        if (t > tout){
            snapshot += 1;
            // output snapshot
            tprint = (double) snapshot*tout;
            std::cout << "snapshot num: " << snapshot << std::endl;
            output_mesh(Mesh,t,snapshot);
            snapshot += 1;
            // output snapshot
            tprint = (double) snapshot*tout;

        }

        // GET DT
        dt = get_dt(Mesh);
        std::cout<< "dt (CFL): "<< dt << std::endl;
        std::cout<< "t = "<< t << std::endl;
        //dt = 0.0001;
        //exit(0);
        //Mesh[0].printPrims(); 
                
        // GET Gradients
        getGradients(Mesh);
        //std::cout << "gradients "; 
        //Mesh[127].printGrad();
        //Mesh[0].printGrad();
        //exit(0);
        // Extrapolate half-step in time (predictor) 
        extrap_to_time(Mesh,dt);
        //Mesh[0].printPrims(); 
        //exit(0);
        // Extrapolate face (corrector)
        extrap_to_face(Mesh);
        //std::cout << "faces: " << std::endl;
        //Mesh[0].printGrad();
        //Mesh[0].printFaces();
        //Mesh[127].printFaces();
        //Mesh[1].printFaces();
        //exit(0);
        // Compute flux 
        getFlux(Mesh,dt);
        //std::cout << "Flux: " << std::endl;
        //exit(0);
        //Mesh[1].printFlux();
        //Mesh[127].printFlux();
        // Add flux to cons 
        Mesh[0].printCons();
        // Get primatives 
        for (int i=0; i < N*N; i++)
        {
            Mesh[i].cons2prim();
        }
        Mesh[0].printPrims();
        //exit(0); 
        t += dt;
    }
    
}
