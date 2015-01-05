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


/*! \class Probleme
 * \brief Classe représentant le problème posé dans son intégralité, vu par une partition.
 */
class Probleme {

    // membres
private :
    Maillage* maillage; /*!< Pointeur vers un objet de la classe maillage*/
    VectorXd* uexa ; /*!< Solution exacte du problème étudié (Utilisé pour la visualisation de l'erreur)*/
    VectorXd* g ; /*!< Nombre de noeuds du maillage*/
    VectorXd* u ; /*!< Solution du problème obtenue itérativement par Jacobi*/
    VectorXd* felim ; /*!<Second membre après pseudo-élimination*/

   
    int* partition_noeud; /*!<Vecteur qui a en position i la partition du noeud i+1, valant 0 si
   * le noeud est sur l'interface */

    
    vector<vector<int> > voisins_interface; /*!<Vecteur de vecteurs ; le ième vecteur contient l'ensemble
   *  des points de la partition i voisins de l'interface */

   
    vector<vector<int> > voisins_partition; /*!<Vecteur de vecteurs ; le ieme vecteur contient l'ensemble
    * des points de l'interface voisins de la partition i */

    Eigen::SparseMatrix<double> *p_K ; /*!< Matrice de rigidité K*/
    Eigen::SparseMatrix<double> *p_M; /*!< matrice de Masse M*/
    Eigen::SparseMatrix<double> *p_Kelim ; /*!<Matrice de rigidité après pseudo elimination */
    Eigen::SparseMatrix<double> *diag; /*!< Diagonale de la matrice de rigidité */

    // methodes
public :
    /*!
    * \brief Constructeur à partir d'un maillage
    *
    * Constructeur du problème pour la géométrié donnée par le maillage, dans la partition donnée 
    * par le second argument rang
    *
    * \param monMaillage Maillage correspondant à la géométrie du problème
    * \param rang Rang de la tache en cours d'éxecution, correspond à une partition du maillage 
    */
    Probleme(Maillage , int);

    /*!
    *\brief Assemblage de la matrice de rigidité
    *
    *Calcule et assemble la matrice de rigidité du problème pour les sommets du maillage 
    * concernés par la partition correspondant à la tache en cours d'execution.
    *
    *\param rang Rang de la tache en cours d'execution = Partition considérée
    */
    void assemblage(int);
    
    /*!
    *\brief Calcul des coefficients du bloc élémentaire de la matrice de rigidité pour un élément du maillage
    *
    * Méthode des coefficients élémentaires de la matrice affectés par l'élément en cours d'analyse
    *
    *\param tab_interm Tableau des distances point à point entre les points le l'élément analysé
    *\param p_K_elem Tableau de 9 double, pouvant être vu comme une matrice 3x3 corrempondant à la matrice élémentaire de l'élément
    *\param aire_triangle Aire de l'élément analysé.
    */
    void mat_K_elem(double* , double* , double );

    /*!
    * \brief Assemblage des blocs élémentaires dans la matrice de rigidité
    *
    * Méthode d'assemblage des rigidités élémentaires pour former la matrice globale de la partition
    *
    *\param ind_triangle Index de l'élément triangulaire dont la matrice de rigidité élémentaire va être ajoutée à la matrice globale
    *\param p_K_elem Matrice de rigidité élémentaire de l'élément
    *\param rang Rang de la tache en cours d'execution = Partition considérée
    */
    void assemblage_pKelem(int , double* , int );

    

 /*!
 * \brief Fonction de calcul des voisins de chaque partition
 * 
 * Calcule les attributs de classe voisins_interface et voisins_partition en parcourant les noeuds de l'interface.
 */
    void calcul_voisins();

  /*!
   * \brief Méthode de Set pour l'attribut Maillage
   *  
   * Méthode de Set générique pour assurer l'encapsulation de la classe.
   *
   *Chaque attribut de la classe dispose de sa methode Set propre (cf source)
   */ 

    void Set_maillage (Maillage* _maillage);
     #ifndef DOXYGEN_SHOULD_SKIP_THIS
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
 #endif /* DOXYGEN_SHOULD_SKIP_THIS */
    //Méthodes "Get" pour l'encapsulation
    //** Méthodes GET
 /*!
  * \brief Méthode de Get pour l'attribut Maillage
  *  
  *Méthode de Get générique pour assurer l'encapsulation de la classe.
  *
  *Chaque attribut de la classe dispose de sa methode Set propre (cf source)
  */
    Maillage* Get_maillage ();
     #ifndef DOXYGEN_SHOULD_SKIP_THIS
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
 #endif /* DOXYGEN_SHOULD_SKIP_THIS */
    //Destructeur

 /*! \brief Destructeur de la classe maillage
  *
  */
    ~Probleme();

     /*! \brief Calcule le second membre du problème pour un point donné
     *
  *
  *Calcule la solution exacte pour un couple x,y de coordonées passé en argument
  *
  *\param coor_1 Abscisse du point dont on calcule la valeur exacte de la solution
  *\param coor_2 Ordonée du point dont on calcule la valeur exacte de la solution
  */
    double calcul_f(double, double);

     /*! \brief Calcule la solution exacte du problème pour un point donné
     *
  *
  *Calcule la solution exacte pour un couple x,y de coordonées passé en argument
  *
  *\param coor_1 Abscisse du point dont on calcule la valeur du second membre
  *\param coor_2 Ordonée du point dont on calcule la valeur du second membre
  */
    double calcul_uexa(double, double);

         /*! \brief Calcule la valeur du relèvenement pour un point donné
     *
  *
  *Calcule la solution exacte pour un couple x,y de coordonées passé en argument
  *
  *\param coor_1 Abscisse du point dont on calcule la valeur du relevement
  *\param coor_2 Ordonée du point dont on calcule la valeur du relevement
  */
    double calcul_g(double, double);


} ;
#endif //PROBLEME_H
