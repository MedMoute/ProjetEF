#include "../include/assemblage.h"

void assemblage(Eigen::SparseMatrix<double> & mat, double* mat_elem, int* tab,int n){

  for(int sommet_1=0;sommet_1<3;sommet_1++)
  {
      int ind_global_1=tab[3*n+sommet_1];

      for(int sommet_2;sommet_2<3;sommet_2++)
      {
          int ind_global_2=tab[3*n+sommet_2];
          mat.coeffRef(ind_global_1,ind_global_2)+=mat_elem[3*sommet_1+sommet_2];
      }
  }
 
}
