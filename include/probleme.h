#ifndef PROBLEME_H
#define PROBLEME_H

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

    /*partition_noeud est un vecteur qui a en position i la partition du noeud i, valant 0 si
   * le noeud est sur l'interface */
    int* partition_noeud;

    /*voisins_interfaces est un vecteur de vecteurs ; le ième vecteur contient l'ensemble
   *  des points de la partition i voisins de l'interface */
    vector<vector<int> > voisins_interface;

    /*voisins_partition est un vecteur de vecteurs ; le ieme vecteur contient l'ensemble
    * des points de l'interface voisins de la partition i */
    vector<vector<int> > voisins_partition;


    //int *tab_local_global;
    Eigen::SparseMatrix<double> *p_K ; // matrice K
    Eigen::SparseMatrix<double> *p_Kelim ; // matrice apres pseudo elimination

    // methodes
public :
    // constructeur a partir d'un maillage
    Probleme(Maillage & , int);
    void affich(Eigen::SparseMatrix<double>);
    void affichVector(VectorXd);
    void assemblage(int);
    void assemblage_felim(double* , double , double , double , double , double , double , double
                                    , int , int );
    void mat_K_elem(double* , double* , double );
    void assemblage_pKelem(int , double* , int );

    /*Calcule, par parcours des noeuds de l'interface, voisins_interface et voisins_partitions */
    void calcul_voisins();

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
