#include <fstream>
#include <iostream>
#include "../include/nonParallel.h"

void calcul(VectorXd u, VectorXd u_nouveau, Eigen::SparseMatrix<double> p_Kelim, Eigen::DiagonalMatrix<double, Eigen::Dynamic>
             diagonale, VectorXd felim)
{
    u_nouveau = diagonale.inverse() * (p_Kelim * u + felim);
}


