/*!
 * \file maillage.h
 * \brief Header du lecteur de fichier.msh issu de GMSH et parseur pour la résolution
 * \version 1.0
 */

#ifndef MAILLAGE_H
#define MAILLAGE_H

#include <iostream>
#include <fstream>
#include <cstdlib>

/*! \namespace std
*
* \brief Espace de nommage standard
*
 * Utilise le namespace de la bibliothèque standard 
 */
using namespace std;

 /*!
 * \fn int type2nnodes(int type)
 * \brief Renvoie le nombre de sommets d'un type d'éléments
 *
 * Fonction renvoyant la correspondance entre le type d'element passé en entée 
 * et le nombre de sommets de cet element.
 *
 * \param type Entier caractérisant un type d'élement dans la syntaxe GMSH
 */

int type2nnodes(int type);

/*! \class Maillage
 * \brief Classe représentant le maillage spatial décomposé
 */
class Maillage {
  
  // membres
 private :
  int n_nodes ;     /*!< Nombre de noeuds du maillage*/
  int n_elems ;     /*!< Nombre d'elements du maillage*/
  int n_triangles ; /*!< Nombre de triangles parmi les elements (exclut les segments) */
  int nb_partitions; /*<Nombre de partitions dans le maillage */
  double* nodes_coords ; /*!< Coordonées des noeuds */
  int* nodes_ref ;   /*!< Index reférencant les noeuds */
  int* elems_type ;  /*!< Liste du type de chaque element (donne son nombre de sommets via n2nodes() )*/
  int* elems_ref ;   /*!< Index référencant les elements*/
  int* elems_sommets ; /*!< Réference des sommets pour chaque element*/
  int* triangles_ref ; /*!< Index référencant les elements triangulaires (exclut les segments)*/
  int* triangles_sommets ; /*!< Donne accès à l'index global des trois sommets de chaque triangle*/
  int* n_partition; /*!< 5eme tag de l'element -> Nombre de partition de l'élément */
  int* partition_ref; /*!< 6eme tag de l'element ->Partition "principale" de l'élément*/


  // methodes
 public :
  /*! \brief Constructeur a  partir d'un fichier
  *
  * \param file : fichier d'entrée de type .msh
  */
  Maillage(ifstream& file);

    /*! \brief Destructeur de la classe maillage
  *
  */
  ~Maillage() ;
  /*!
  * \brief Méthode de Set pour l'attribut n_nodes
  *  
  *Méthode de Set générique pour assurer l'encapsulation de la classe.
  *
  *Chaque attribut de la classe dispose de sa methode Set propre (cf source)
  */


  void Set_n_nodes (int _n_nodes);
 #ifndef DOXYGEN_SHOULD_SKIP_THIS
  void Set_n_elems (int _n_elems);

  void Set_n_triangles (int _n_triangles);
 
  void Set_nb_partitions (int _n_partitions);

  void Set_nodes_coords (double* _nodes_coords);

  void Set_nodes_ref (int* _nodes_ref);

  void Set_elems_type (int* _elems_type);

  void Set_elems_ref (int* _elems_ref);
  
  void Set_elems_sommets (int* _elems_sommets);
  
  void Set_triangles_ref (int* _triangles_ref);
    
  void Set_triangles_sommets (int* _triangles_sommets);
 
  void Set_n_partition (int* _n_partition);
  
  void Set_partition_ref (int* _partition_ref);

  #endif /* DOXYGEN_SHOULD_SKIP_THIS */
  /*!
  * \brief Méthode de Get pour l'attribut n_nodes
  *
  *Méthode de Get générique pour assurer l'encapsulation de la classe
  *
  *Chaque attribut de la classe dispose de sa methode Get propre (cf source)
  */
  int Get_n_nodes ();
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  int Get_n_elems ();

  int Get_n_triangles();
  
  int Get_nb_partitions();
  double* Get_nodes_coords();
  int* Get_nodes_ref();
  int* Get_elems_type();
  int* Get_elems_ref();
  
  int* Get_elems_sommets();
  
  int* Get_triangles_ref();
  int* Get_triangles_sommets();

  int* Get_n_partition();
  
  int* Get_partition_ref();
  #endif /* DOXYGEN_SHOULD_SKIP_THIS */

};
#endif // MAILLAGE_H

