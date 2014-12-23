#ifndef PARALLEL_H
#define PARALLEL_H

#include <Eigen/Sparse>
#include "probleme.h"

extern int rang;
extern int nb_procs;

void initialisation_mpi( int, char ** );
void communication(VectorXd , vector<vector<int> >, vector<vector<int> >);
double erreur_entre_etapes (VectorXd, VectorXd);
void finalisation_mpi();


#endif // PARALLEL_H
