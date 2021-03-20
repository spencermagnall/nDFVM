#ifndef INIT_H
#define INIT_H
#include "meshdata.h"
#include "vars.h"
// Setup the mesh spacing 
void init_mesh(meshcell* &mesh);
// Setup the fluid
void init_fluid(meshcell* &mesh);

#endif 
