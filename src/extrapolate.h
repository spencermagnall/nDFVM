#ifndef EXTRAP_H
#define EXTRAP_H
#include "meshdata.h"
#include "vars.h"

void extrapolate_time(meshcell* mesh, double dt)
{
    for (int cellno=0; cellno<N*N; cellno++)
    {
        mesh[cellno].extraptime();
    }

}
void extrapolate_face()
#endif 
