#ifndef PROBLEME_H
#define PROBLEME_H

#include "maillage.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>
#include <algorithm>
#define PARALLELE 0
#define it_max 10000
#define eps 2.e-5
#define PI 3.14159

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd ;

class Probleme {

    // membres
private :
    Maillage* maillage;
    VectorXd* uexa ;
    VectorXd* g ;
    VectorXd* u ;
    VectorXd* felim ;

    /*partition_noeud est un vecteur qui a en position i la partition du noeud i+1, valant 0 si
   * le noeud est sur l'interface */
    int* partition_noeud;

    /*voisins_interface est un vecteur de vecteurs ; le ième vecteur contient l'ensemble
   *  des points de la partition i voisins de l'interface */
    vector<vector<int> > voisins_interface;

    /*voisins_partition est un vecteur de vecteurs ; le ieme vecteur contient l'ensemble
    * des points de l'interface voisins de la partition i */
    vector<vector<int> > voisins_partition;

    //int *tab_local_global;
    Eigen::SparseMatrix<double> *p_K ; // matrice K
    Eigen::SparseMatrix<double> *p_M; // matrice de masse
    Eigen::SparseMatrix<double> *p_Kelim ; // matrice apres pseudo elimination
    Eigen::SparseMatrix<double> *diag; // diagonale de la matrice de rigidité

    // methodes
public :
    // constructeur a partir d'un maillage
    Probleme(Maillage , int);
    void assemblage(int);
    
    void mat_K_elem(double* , double* , double );
    void assemblage_pKelem(int , double* , int );

    /*Calcule, par parcours des noeuds de l'interface, voisins_interface et voisins_partition */

    void calcul_voisins();

    //Méthodes de Set pour l'encapsulation
    void Set_maillage (Maillage* _maillage);
    void Set_uexa (VectorXd* _uexa);
    void Set_g (VectorXd* _g);
    void Set_u (VectorXd* _u);
    void Set_felim (VectorXd* _felim);

    void Set_partition_noeud (int* _partition_noeud);
    void Set_voisins_interface (vector<vector<int> > _voisins_interface);
    void Set_voisins_partition (vector<vector<int> > _voisins_partition);

    void Set_p_K (Eigen::SparseMatrix<double>* _p_K);
    void Set_p_M (Eigen::SparseMatrix<double>* _p_M);

    void Set_p_Kelim (Eigen::SparseMatrix<double>* _p_Kelim);
    void Set_diag (Eigen::SparseMatrix<double>* _diag);

    //Méthodes "Get" pour l'encapsulation

    Maillage* Get_maillage ();
    VectorXd* Get_uexa ();
    VectorXd* Get_g ();
    VectorXd* Get_u ();
    VectorXd* Get_felim ();

    int* Get_partition_noeud ();
    vector<vector<int> > Get_voisins_partition ();
    vector<vector<int> > Get_voisins_interface ();

    Eigen::SparseMatrix<double>* Get_p_K ();
    Eigen::SparseMatrix<double>* Get_p_M ();

    Eigen::SparseMatrix<double>* Get_p_Kelim ();
    Eigen::SparseMatrix<double>* Get_diag ();

    //Destructeur
    ~Probleme();

    void pelim(int) ;

    double base_loc(int, double, double);
    double calcul_f(double, double);
    double calcul_uexa(double, double);
    double calcul_g(double, double);

    void assemblage_par(Eigen::SparseMatrix<double> &, double *, int*, int) ;
    void mat_K_elem_par(double *, double *, int *, int);

} ;
#endif //PROBLEME_H
