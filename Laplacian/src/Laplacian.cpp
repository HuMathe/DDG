#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

#include "Types.h"
#include "Mesh.h"

typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;
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
            if ( he->onBoundary ) {continue;}

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
        /* NOTE: here what we get is the operator d*d rather than *d*d, namely a 2-form.
         * The reason comes from the symmetry of this matrix, and the computation of many 
         * quantities ask for *(*d*d)= d*d 
         * */

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
       Vec A(vertices.size());  /* area of mesh */
       
       for (VertexCIter v = vertices.begin(); v != vertices.end(); v++)
       {
           rho( v->index ) = v->rho;
           A  ( v->index ) = v->area();
       }
       
       /* make integration of rho on mesh to be 0 */
       rho = rho - ( rho.dot(A) / A.sum() ) * Vec::Ones(rho.size());

       // debug use
       printf("normed rho: max: %f, min: %f\n", rho.maxCoeff(), rho.minCoeff());

       /* solve for 2-form of d*d rho */
       phi = LapSolver.solve( (rho.array() * A.array()).matrix());

       if (LapSolver.info() != Eigen::Success)
       {
           printf("fail to solve equation.\n");
           return;
       }

       for (VertexIter v = vertices.begin(); v != vertices.end(); v++)
       {
           v->phi = phi( v->index ); 
       }
    }


    /* calculate a step forward to soomth mesh according to curvature flow */
    void Mesh::soomthMesh(double h)
    {
        Mat Pos(vertices.size(), 3);
        Mat deltaPos(vertices.size(), 3);

        /* assign value to Pos */
        for (VertexCIter v = vertices.begin(); v != vertices.end(); v++)
        {
            Pos(v->index, 0) = v->position[0];
            Pos(v->index, 1) = v->position[1];
            Pos(v->index, 2) = v->position[2];
        }

        /* at each time the location of the vertices are different and thus 
         * Laplacian should be updated 
         * */
        buildLaplacian();

        Vec A(vertices.size());  /* area of mesh */
        for (VertexCIter v = vertices.begin(); v != vertices.end(); v++)
        {
           A(v->index) = v->area();
        }

        /* curvature flow */
        deltaPos = h * Laplacian * Pos;

        /* assign value to Pos */
        for (VertexIter v = vertices.begin(); v != vertices.end(); v++)
        {
            v->position[0] += deltaPos(v->index, 0) / A( v->index );
            v->position[1] += deltaPos(v->index, 1) / A( v->index );
            v->position[2] += deltaPos(v->index, 2) / A( v->index );
        }

        printf("Assigned new position to vertices.\n");
    }
}
