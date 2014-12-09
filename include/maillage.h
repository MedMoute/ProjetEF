#include <iostream>
#include <fstream>
#include <cstdlib>
#include <map>

using namespace std;

class Maillage {
  
  // membres
 public :
  int n_nodes ;     // nb de noeuds
  int n_elems ;     // nb d'elements 
  int n_triangles ; // nb de triangles parmi les elements
  int nb_partitions;
  double* nodes_coords ; // coord. des noeuds
  int* nodes_ref ;   // ref. des noeuds
  int* elems_type ;  // type de chaque element et nb de sommets
  int* elems_ref ;   // ref. des elements
  int* triangles_ref ; // C EST QUOI
  int* triangles_sommets ; // donne acces au num global d un des trois sommets de chaque triangle
  int* elems_sommets ; // ref. des sommets de chaque element
  int* n_partition; // tmp[5] = 5eme tag de l'element -> nbre de partition de l'élément X
  int* partition_ref; //tmp[6] = 6eme tag de l'element ->partition "principale" de l'élément

  //TODO
  std::map<int,int>** partition_map_voisins; // Pointeur vers des tableaux de Map (indice du point,partition du point)
                                     // qui à chaque point associe la map des voisins et leur partition associée avec 0 pour
                                     // l'interface

  // methodes
 public :
  // constructeur a partir d'un fichier
  Maillage( ifstream& );

  ~Maillage() ;
};
