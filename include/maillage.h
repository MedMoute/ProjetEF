#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

class Maillage {
  
  // membres
 public :
  int n_nodes ;     // nb de noeuds
  int n_elems ;     // nb d'elements 
  int n_triangles ; // nb de triangles parmi les elements
  double* nodes_coords ; // coord. des noeuds
  int* nodes_ref ;   // ref. des noeuds
  int* elems_type ;  // type de chaque element et nb de sommets
  int* elems_ref ;   // ref. des elements
  int* triangles_ref ; // C EST QUOI
  int* triangles_sommets ; // donne acces au num global d un des trois sommets de chaque triangle
  int* elems_sommets ; // ref. des sommets de chaque element
  int* n_partition;
  int* partition_ref;
  // methodes
 public :
  // constructeur a partir d'un fichier
  Maillage( ifstream& );

};
