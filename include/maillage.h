#ifndef MAILLAGE_H
#define MAILLAGE_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <map>

using namespace std;

class Maillage {
  
  // membres
 private :
  int n_nodes ;     // nb de noeuds
  int n_elems ;     // nb d'elements 
  int n_triangles ; // nb de triangles parmi les elements
  int nb_partitions;
  double* nodes_coords ; // coord. des noeuds
  int* nodes_ref ;   // ref. des noeuds
  int* elems_type ;  // type de chaque element et nb de sommets
  int* elems_ref ;   // ref. des elements
  int* elems_sommets ; // ref. des sommets de chaque element
  int* triangles_ref ; // C EST QUOI
  int* triangles_sommets ; // donne acces au num global d un des trois sommets de chaque triangle
  int* n_partition; // tmp[5] = 5eme tag de l'element -> nbre de partition de l'élément X
  int* partition_ref; //tmp[6] = 6eme tag de l'element ->partition "principale" de l'élément

  /*
  //TODO
  std::map<int,int>** partition_map_voisins; // Pointeur vers des tableaux de Map (indice du point,partition du point)
                                     // qui à chaque point associe la map des voisins et leur partition associée avec 0 pour
                                     // l'interface

                                     */
  // methodes
 public :
  // constructeur a partir d'un fichier
  Maillage(ifstream&);
  ~Maillage() ;
  // Méthodes set pour l'encapsulation de la classe
  void Set_n_nodes (int _n_nodes);
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

  //Méthodes get pour l'encapsulation 
  int Get_n_nodes ();
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

  //Destructeur

};
#endif // MAILLAGE_H

