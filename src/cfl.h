#ifndef CFL_H
#define CFL_H
#include <algorithm> 
#include <cmath>
#include "vars.h"
#include "meshdata.h"

double get_dt(meshcell* mesh)
{
    double courant_fac = 0.4;
    double cs;
    double v;
    double dt = HUGE_VAL;
    double csmax = -HUGE_VAL;
    for (int i=0; i < N*N; i++){
        cs = mesh[i].getSoundspeed();
        v = mesh[i].getVelocity();
        //std::cout << "soundspeed: " << cs << std::endl;
        //std::cout << "velocity: " << v << std::endl;
        if (dx/(cs + v) < dt){
            dt = dx/(cs+ v);
        }

        if (cs+v > csmax){
            csmax = cs+v;
        }
    }
    std::cout<< "Cs max is: " << csmax << std::endl;
    std::cout<< "dt: " << dt << std::endl;
    return dt*courant_fac;
}

#endif // CFL_H
