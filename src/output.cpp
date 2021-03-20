#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
//#include <limits>
#include "vars.h"
#include "output.h"

void output_mesh(meshcell Mesh[], float time, int snapshot){
    // Get file name 
     std::string snapshotnum = std::to_string(snapshot);
     std::string snapshotnumnew = std::string(5-snapshotnum.length(), '0') +     snapshotnum;
     std::string filename  = "snapshot_" + snapshotnumnew;
    // Open/Create file
     std::cout << "Writing file: " << filename << std::endl;
     std::ofstream file (filename);
    
    // Write header
    if (file.is_open()){
        
        file << "x, "<< "y, " << "vx, " << "vy, " << "P, " << "rho" << std::endl; 
        file << time << std::endl;
        for (int cellno=0; cellno<N*N; cellno++){
            file << std::setprecision(17) << Mesh[cellno].printCell() <<std::endl;
        }       

    }
    file.close();
}
