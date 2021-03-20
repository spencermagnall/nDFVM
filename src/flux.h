#ifndef FLUX_H
#define FLUX_H
#include "vars.h"
#include "meshdata.h"

//void rusanov(meshcell* mesh, double dt);
void getFlux(meshcell* mesh, double dt);

void extrap_to_time(meshcell* mesh, double dt);

void extrap_to_face(meshcell* mesh);
#endif // FLUX_H
