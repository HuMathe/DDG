#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Vertex.h"

namespace DDG
{
    void setNorm(int state)
    {
        Vertex::normState = state;
        printf("normState has been changed to %d \n", Vertex::normState);
    }

}

#endif
