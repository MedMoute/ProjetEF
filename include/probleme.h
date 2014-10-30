#include "maillage.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

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
  Probleme(Maillage & ) ;
  void affich(Eigen::SparseMatrix<double>);
  void affichUnElem(Eigen::SparseMatrix<double>, int, int);
  
  ~Probleme();




  void pelim(int) ;

  double base_loc(int, double, double);
  double calcul_f(double, double);
  double calcul_uexa(double, double);
  double calcul_g(double, double);

  void assemblage_par(Eigen::SparseMatrix<double> &, double *, int*, int) ;
  void mat_K_elem_par(double *, double *, int *, int);

} ;
