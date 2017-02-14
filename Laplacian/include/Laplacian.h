#pragma once

#include "Types.h"
#include "Vector.h"
#include "Mesh.h"

#include <Eigen/Sparse>

#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

namespace DDG
{
    
    /* calculate color according to potential phi */
    void calcPotentialColor(Mesh &mesh);

    /* randomly assign value to rho */
    void randomAssignRho(Mesh &mesh);

    class Laplacian
    {
        typedef Eigen::SparseMatrix<double> spMat;

        #define _CHOL_TYPE_ 3

        #if _CHOL_TYPE_ == 1
            typedef Eigen::SparseLU<spMat> cholSolver;
        #elif _CHOL_TYPE_ == 2
            typedef Eigen::SimplicialLLT<spMat> cholSolver;
        #elif _CHOL_TYPE_ == 3
            typedef Eigen::CholmodSimplicialLLT<spMat> cholSolver;
        #elif _CHOL_TYPE_ == 4
            typedef Eigen::CholmodSupernodalLLT<spMat> cholSolver;
        #endif

        /* build laplacian matrix of the mesh */
        void buildLaplacian( Mesh& mesh );

        /* solve the scalar poisson problem on mesh */
        void solveScalarPoissonProblem( Mesh& mesh );

        /* calculate a step forward to soomth mesh according to curvature flow */
        void soomthMesh(Mesh& mesh, double h);
      
        protected:

        spMat Laplacian;
        // Laplacian sprase matrix of the mesh

        cholSolver negLapSolver;
        // cholmode solver based on SuiteSparse
    };

}


