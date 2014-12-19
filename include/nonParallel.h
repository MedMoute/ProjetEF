#ifndef NONPARALLEL_H
#define NONPARALLEL_H

#include <Eigen/Sparse>
#include "probleme.h"

void calcul(VectorXd, VectorXd, Eigen::SparseMatrix<double>, Eigen::DiagonalMatrix<double, Eigen::Dynamic>, VectorXd);

#endif // NONPARALLEL_H
