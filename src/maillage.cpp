
/*!
 * \file maillage.cpp
 * \brief Source du lecteur de fichier.msh issu de GMSH et parseur pour la résolution
 * \version 1.0
 */

#include "../include/maillage.h"
#include <string>
#include <sstream>
#include <iostream>

int type2nnodes(int) ;

Maillage::Maillage(ifstream& fd) {


    // POUR LES FICHIERS VERSION 2.2 SEULEMENT
    string line ;
    // on lit dans le fichier fd
    if(fd.is_open() ) {
        // on lit les lignes une par une et on avise en fonction des mots cle
        while ( fd.good() ) {
            getline (fd,line);
            // on regarde si on a un mot cle
            if( line[0] == '$' ) {

                if( line.compare("$MeshFormat")==0 ) {
                    // je lis le format mais je n'en fais rien
                    getline(fd,line) ;
                    // je lis le mot cle de fin
                    getline(fd,line) ;
                    if( ! line.compare("$EndMeshFormat")==0 ) {
                        cout << "Probleme de lecture du fichier maillage" << endl ;
                        cout << "Mot cle MeshFormat non suivi de EndMeshFormat" << endl ;
                        abort() ;
                    }
                }
                else if( line.compare("$Nodes")==0 ) {
                    // lecture du nombre de noeuds
                    fd >> n_nodes ;
                    // on alloue le tableau de coord des noeuds
                    nodes_coords = new double[3*n_nodes] ;
                    nodes_ref = new int[n_nodes] ;

                    // on lit les coords
                    for (int i = 0; i< n_nodes ; i++)
                    {
                        // variable dummy pour le numero de noeud
                        double dummy ;
                        fd >> dummy >> nodes_coords[3*i] >> nodes_coords[3*i+1] >> nodes_coords[3*i+2] ;
                        nodes_ref[i] = 0 ;
                    }
                    // on passe a la ligne (je ne sais pas pourquoi c'est necessaire
                    getline(fd,line) ;
                    // on lit le mot cle de fin
                    getline(fd,line) ;
                    if( ! line.compare("$EndNodes")==0 ) {
                        cout << "Probleme de lecture du fichier maillage" << endl ;
                        cout << "Mot cle Nodes non suivi de EndNodes" << endl ;
                        abort() ;
                    }

                }
                else if( line.compare("$Elements")==0 ) {
                    // on lit le nb d'elements
                    fd >> n_elems ;
                    // on passe a la ligne pour la suite
                    getline(fd,line) ;
                    // allocation
                    elems_ref = new int[n_elems] ;
                    n_partition = new int[n_elems] ;
                    partition_ref = new int[n_elems] ;
                    elems_type = new int[2*n_elems] ;
                    this->n_triangles=0;
                    nb_partitions=0;

                    // je considere qu'on n'a que des segments et des triangles et je fais un tableau simple
                    elems_sommets = new int[3*n_elems] ;

                    for( int i=0; i<n_elems ; i++)
                    {


                        // on ne sais pas combien il y a de parametres sur chaque ligne : on lit ligne par ligne
                        getline(fd,line) ;
                        // on compte les espaces pour savoir combien on a de parametres
                        int nspace=0;
                        for(unsigned int j=0; j!=line.size(); ++j)
                            nspace+=( line.at(j)==' ');
                        int ntmp = nspace+1 ;
                        int* tmp = new int[ntmp] ;
                        // on met les entiers dans les cases
                        size_t spos ;
                        size_t isize ;
                        size_t sdeb = 0 ;
                        int ind = 0 ;
                        do {
                            spos = line.find(' ', sdeb) ;
                            isize = spos - sdeb ;
                            istringstream(line.substr(sdeb, isize)) >> tmp[ind] ;
                            ind ++;
                            sdeb = spos+1 ;
                        } while (spos!=string::npos ) ;
                        // on a tous les entiers dans tmp
                        // on remplit les membres de maillage
                        elems_type[2*i] = tmp[1] ; // type d'elements
                        elems_type[2*i+1] = type2nnodes(tmp[1]);  // nb de sommets pour cet element


                        // GMSH met les elements du bord en 1er
                        // la ref du 1er element correspond donc a la ref du bord
                        int n_tags = tmp[2] ;
                        if (n_tags > 2) {
                            elems_ref[i] = tmp[3] ;
                            n_partition[i] = tmp[5] ;
                            partition_ref[i] = tmp[6] ;
                            if (tmp[6]>nb_partitions)
                            {
                                nb_partitions=tmp[6];
                            }
                        }

                        for (int j=0;j<elems_type[2*i+1] ; j++) {
                            elems_sommets[3*i+j] = tmp[3+n_tags+j] ;
                            if (tmp[1]==2)
                            {
                            }
                            // on remplit nodes_ref : le -1 est la pour recoller a la numerotation 0..n-1
                            if(nodes_ref[ tmp[3+n_tags+j]-1 ] == 0)
                                nodes_ref[ tmp[3+n_tags+j]-1 ] = elems_ref[i] ;
                        }

                        if (elems_type[2*i+1]==3){
                            n_triangles++;
                        }

                        // TODO : Remplir partition_map et partition_map_voisins
                        /*
                         *Methode : lors de la lecture des élément par le parser, les tags donnent une partition principale
                         *dans le cas ou le point d'indice i n'est pas mappé i.e. partition_map[i].empty vrai,on associe les
                         *autres points de l'élément à la partition principale.
                         *
                         *Si il est mappé , on ajoute les points voisins non-encore existants automatiquement.
                         *Si un point voisin est déja mappé, mais avec une valeur différente, cela veut dire qu'il est sur l'interface
                         *
                         *Une fois la boucle effectuée sur tous les points, on peut merger toutes les maps et on obtient une map*<int,int>
                         *donnant le (indice,partition) global [
                        */

                        /*initialiser partition_map_voisins
                        for (int j=0;j<n_elems;j++) {
                            if (partition_map_voisins[tmp[3]]->empty())

                        }*/


                        delete [] tmp ;
                    }
                    getline(fd,line) ;
                    if( ! line.compare("$EndElements")==0 ) {
                        cout << "Probleme de lecture du fichier maillage" << endl ;
                        cout << "Mot cle Nodes non suivi de EndNodes" << endl ;
                    }

                    triangles_sommets = new int[3*n_triangles];
                    for (int k=0;k<n_triangles;k++)
                    {
                        for (int l=0;l<3;l++)
                        {


                            triangles_sommets[3*k+l]=elems_sommets[3*k+l+3*(n_elems-n_triangles)];
                            int a=elems_sommets[3*k+l+3*(n_elems-n_triangles)];
                            //std::cout<<a<<" | ";


                            //attention on suppose qu'il n'y a que des segments et des triangles, et que les segments sont tous avant les triangles dans le .msh
                        }
                        //std::cout<<std::endl;
                    }
                }
            }
        }
    }
}

Maillage::~Maillage()
{
  /*
    delete [] triangles_sommets;
    delete [] nodes_coords;
    delete [] nodes_ref;
    delete [] elems_ref;
    delete [] elems_sommets;
    delete [] elems_type;
    delete [] n_partition;
    delete [] partition_ref;
    */
}

// correspondance entre le type d'element et le nb de sommets de cet element
int type2nnodes(int type) {

    int nnodes[31] = {2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 4, 5, 6, 20, 35, 56} ;
    return nnodes[type-1] ;

}

//// Méthdes SET et GET pour l'encapsulation
  //** Méthodes SET 
  
  void Maillage::Set_n_nodes (int _n_nodes) {
    n_nodes=_n_nodes;
  }
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  void Maillage::Set_n_elems (int _n_elems) {
    n_elems=_n_elems;
  }

  void Maillage::Set_n_triangles (int _n_triangles) {
    n_triangles=_n_triangles;
  }

  void Maillage::Set_nb_partitions (int _nb_partitions) {
    nb_partitions=_nb_partitions;
  }

  void Maillage::Set_nodes_coords (double* _nodes_coords) {
    nodes_coords=_nodes_coords;
  }

  void Maillage::Set_nodes_ref (int* _nodes_ref) {
    nodes_ref=_nodes_ref;
  }

  void Maillage::Set_elems_type (int* _elems_type) {
    elems_type=_elems_type;
  }

  void Maillage::Set_elems_ref (int* _elems_ref) {
    elems_ref=_elems_ref;
  }

  void Maillage::Set_elems_sommets (int* _elems_sommets) {
    elems_sommets=_elems_sommets;
  }

  void Maillage::Set_triangles_ref (int* _triangles_ref) {
    triangles_ref=_triangles_ref;
  }

  void Maillage::Set_triangles_sommets (int* _triangles_sommets) {
    triangles_sommets=_triangles_sommets;
  }

  void Maillage::Set_n_partition (int* _n_partition) {
    n_partition=_n_partition;
  }

  void Maillage::Set_partition_ref (int* _partition_ref) {
    partition_ref=_partition_ref;
  }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  //** Méthodes GET
  int Maillage::Get_n_nodes () {
    return n_nodes;
  }
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  int Maillage::Get_n_elems () {
    return n_elems;
  }

  int Maillage::Get_n_triangles() {
    return n_triangles;
  }

  int Maillage::Get_nb_partitions() {
    return nb_partitions;
  }

  double* Maillage::Get_nodes_coords() {
    return nodes_coords;
  }

  int* Maillage::Get_nodes_ref() {
    return nodes_ref;
  }

  int* Maillage::Get_elems_type() {
    return elems_type;
  }

  int* Maillage::Get_elems_ref() {
    return elems_ref;
  }

  int* Maillage::Get_elems_sommets() {
    return elems_sommets;
  }

  int* Maillage::Get_triangles_ref() {
    return triangles_ref;
  }

  int* Maillage::Get_triangles_sommets() {
    return triangles_sommets;
  }

  int* Maillage::Get_n_partition() {
    return n_partition;
  }

  int* Maillage::Get_partition_ref() {
    return partition_ref;
  }

#endif /* DOXYGEN_SHOULD_SKIP_THIS */