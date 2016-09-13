#include <stdlib.h>

#include "Types.h"
#include "Vector.h"
#include "Mesh.h"

namespace DDG
{
    void calcPotentialColor(Mesh &mesh)
    {
        double min = 1, max = -1;

        double abs_max = 0;
        for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
        {
            if (std::abs(v->phi) > abs_max) abs_max = std::abs(v->phi);
        }

        for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
        {
            if (v->phi > max) max = v->phi;
            if (v->phi < min) min = v->phi;

            if (v->phi > 0) {v->color = Vector(v->phi / abs_max, 0, 0);}
            else {v->color = Vector(0, 0, -v->phi / abs_max);}

        }

        printf("max: %f, min: %f\n", max, min);
    }

    void randomAssignRho(Mesh &mesh)
    {
      // initialize all the value
      for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
      {
          v->rho = 0;
          v->phi = 0;
      }

      // randomly add some potential on surface
      int N = 3;
      int VN = mesh.vertices.size();
      for (int i=0; i<N; i++)
      {
          int index = rand() % VN;
          mesh.vertices[index].rho = 1.0;
          mesh.vertices[index].phi = 1.0;
          index = rand() % VN;
          mesh.vertices[index].rho = -1.0;
          mesh.vertices[index].phi = -1.0;
      }
    }

}
