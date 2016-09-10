#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

#include "Types.h"
#include "Mesh.h"

typedef Eigen::VectorXd Vec;
typedef Eigen::Triplet<double> T;

namespace DDG
{
    void Mesh :: buildLaplacian( void )
    // build laplacian matrix of the mesh
    {
        /* assign index to each element */
        this->indexElements();
 
        Laplacian.resize(vertices.size(), vertices.size());

        /* prepare input quene to initialize Laplacian */
        std::vector<T> tripletlist;
        tripletlist.reserve( 4*halfedges.size() );

        for (HalfEdgeCIter he = halfedges.begin(); he != halfedges.end(); he++)
        {
            if ( he->onBoundary ) continue;

            int origin = he->vertex->index;
            int target = he->next->vertex->index;
            double cot_value = he->cotan() / 2.0 ;

            tripletlist.push_back( T(origin, target, cot_value) );
            tripletlist.push_back( T(target, origin, cot_value) );
            tripletlist.push_back( T(origin, origin, -cot_value));
            tripletlist.push_back( T(target, target, -cot_value));
        }

        /* assigin value to Laplacian */
        Laplacian.setFromTriplets(tripletlist.begin(), tripletlist.end());

        /* construct solver */
        this->LapSolver.compute(Laplacian);
        if (LapSolver.info() != Eigen::Success)
        {
            printf("fail to construct solver.\n");
        }
    }

    void Mesh::solveScalarPoissonProblem( void )
    /* solve the scalar poisson problem on mesh */
    {
       Vec rho(vertices.size());
       Vec phi(vertices.size());

       double A = 0;
       for (VertexCIter v = vertices.begin(); v != vertices.end(); v++)
       {
           rho( v->index ) = v->rho * v->area();
           A += v->area();
       }
       
       /* normalize rho */
       double rho_m = rho.sum() / A;
       for (VertexCIter v = vertices.begin(); v != vertices.end(); v++)
       {
           rho( v->index ) -= rho_m * v->area();
       }

       phi = LapSolver.solve(rho);

       if (LapSolver.info() != Eigen::Success)
       {
           printf("fail to solve equation.\n");
           return;
       }

       for (VertexIter v = vertices.begin(); v != vertices.end(); v++)
       {
           v->phi = phi(v->index);
       }
    }
}
