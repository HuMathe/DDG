#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Types.h"
#include "Vector.h"
#include "Mesh.h"

namespace DDG
{
    /* calculate color according to potential phi */
    void calcPotentialColor(Mesh &mesh);

    /* randomly assign value to rho */
    void randomAssignRho(Mesh &mesh);

}

#endif
