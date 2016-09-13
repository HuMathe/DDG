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
        this->negLapSolver.compute( -Laplacian );
        if (negLapSolver.info() != Eigen::Success)
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
       Vec hodgeRho = rho.array() * A.array();
       phi = - negLapSolver.solve( hodgeRho );

       if (negLapSolver.info() != Eigen::Success)
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
        Mat newPos(vertices.size(), 3);

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

        spMat A;  /* area of mesh */
        std::vector<T> area;
        area.reserve(vertices.size());
        for (VertexCIter v = vertices.begin(); v != vertices.end(); v++)
        {
           area.push_back( T(v->index, v->index, v->area()) );
        }
        A.resize(vertices.size(), vertices.size());
        A.setFromTriplets(area.begin(), area.end());

        /* curvature flow */
        cholSolver solver;
        solver.compute(A - h*Laplacian);
        newPos = solver.solve(A*Pos);

        /* assign value to Pos */
        for (VertexIter v = vertices.begin(); v != vertices.end(); v++)
        {
            v->position[0] = newPos(v->index, 0);
            v->position[1] = newPos(v->index, 1);
            v->position[2] = newPos(v->index, 2);
        }

        printf("Assigned new position to vertices.\n");
    }
}
