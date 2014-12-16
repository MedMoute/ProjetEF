#include <fstream>
#include <iostream>
#include "../include/nonParallel.h"

void calcul(VectorXd u, VectorXd u_nouveau, Eigen::SparseMatrix<double> p_Kelim, VectorXd felim)
{
    int taille = u.rows();
    Eigen::DiagonalMatrix<double, taille> diag;
    diag.diagonal() = (p_Kelim).diagonal();
    u_nouveau = diag.inverse() * (((p_Kelim)-diag) * u + felim);
}


