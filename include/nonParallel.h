#ifndef NONPARALLEL_H
#define NONPARALLEL_H

#include <Eigen/Sparse>
#include "probleme.h"


void calcul(VectorXd, VectorXd, Eigen::SparseMatrix<double>, Eigen::DiagonalMatrix<double, Eigen::Dynamic>, VectorXd);
void affich(Eigen::SparseMatrix<double>);
void affichVector(VectorXd);
void affiche_vector(vector<vector<int> >);

#endif // NONPARALLEL_H
