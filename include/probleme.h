#include "maillage.h"
#include <Eigen/Sparse>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd ;

class Probleme {
  
  // membres
 public :
  Maillage* maillage;
  VectorXd* uexa ;
  VectorXd* g ;
  VectorXd* u ;
  VectorXd* felim ;

  //int *tab_local_global;
  Eigen::SparseMatrix<double> *p_K ; // matrice K
  Eigen::SparseMatrix<double> *p_Kelim ; // matrice apres pseudo elimination
  
  // methodes
 public :
  // constructeur a partir d'un maillage
  Probleme( Maillage const& ) ;
  void affich(Eigen::SparseMatrix<double>);
  
  ~Probleme() {
    delete p_K ;
  };




  void pelim(int) ;

  void assemblage(Eigen::SparseMatrix<double> &, double *, int*, int) ;
  void mat_K_elem(double *, double *, int *, int);

} ;
