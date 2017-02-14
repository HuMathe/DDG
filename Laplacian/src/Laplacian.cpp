#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

#include "Types.h"
#include "Mesh.h"
#include "Laplacian.h"

typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;
typedef Eigen::Triplet<double> T;

namespace DDG
{
    void Laplacian :: buildLaplacian( Mesh& mesh )
    // build laplacian matrix of the mesh
    {
        Laplacian.resize(mesh.vertices.size(), mesh.vertices.size());

        /* prepare input quene to initialize Laplacian */
        std::vector<T> tripletlist;
        tripletlist.reserve( 4 * mesh.halfedges.size() );

        for (HalfEdgeCIter he = mesh.halfedges.begin(); he != mesh.halfedges.end(); he++)
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
        negLapSolver.compute( -Laplacian );
        if (negLapSolver.info() != Eigen::Success)
        {
            printf("fail to construct solver.\n");
        }
    }

    void Laplacian::solveScalarPoissonProblem( Mesh& mesh )
    /* solve the scalar poisson problem on mesh */
    {
       Vec rho(mesh.vertices.size());
       Vec phi(mesh.vertices.size());
       Vec A(mesh.vertices.size());  /* area of mesh */
       
       for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
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

       for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
       {
           v->phi = phi( v->index ); 
       }
    }


    /* calculate a step forward to soomth mesh according to curvature flow */
    void Laplacian::soomthMesh(Mesh& mesh, double h)
    {
        Mat Pos(mesh.vertices.size(), 3);
        Mat newPos(mesh.vertices.size(), 3);

        /* assign value to Pos */
        for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
        {
            Pos(v->index, 0) = v->position[0];
            Pos(v->index, 1) = v->position[1];
            Pos(v->index, 2) = v->position[2];
        }

        /* at each time the location of the vertices are different and thus 
         * Laplacian should be updated 
         * */
        buildLaplacian(mesh);

        spMat A;  /* area of mesh */
        std::vector<T> area;
        area.reserve(mesh.vertices.size());
        for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
        {
           area.push_back( T(v->index, v->index, v->area()) );
        }
        A.resize(mesh.vertices.size(), mesh.vertices.size());
        A.setFromTriplets(area.begin(), area.end());

        /* curvature flow */
        cholSolver solver;
        solver.compute(A - h*Laplacian);
        newPos = solver.solve(A*Pos);

        /* assign value to Pos */
        for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
        {
            v->position[0] = newPos(v->index, 0);
            v->position[1] = newPos(v->index, 1);
            v->position[2] = newPos(v->index, 2);
        }

        printf("Assigned new position to vertices.\n");
    }
}
