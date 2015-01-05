/**
 * \file main.cpp
 * \brief Programme de résolution parallélisée de problèmes linéaire par la méthode de Jacobi par points.
 * \author Mehdi Ennaïme & Pierre Fournier
 * \version 1.0
 * \date 5 Janvier 2015
 *
 * Programme utilisant le standard MPI pour paralléliser la résolution d'un problème linéaire par la méthode de Jacobi sur un maillage 
 * non structuré de type éléments finis.
 * 
 */

 #include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "../include/maillage.h"
#include "../include/probleme.h"
#include "../include/nonParallel.h"


int main(int argc, char *argv[])
{
    VectorXd u, u_nouveau, u_avant_com, second_membre,u_exact;
    Eigen::SparseMatrix<double> mat_rigidite;
    Eigen::MatrixXd d_mat_rigidite_local, diag_inv_local;
    Eigen::SparseMatrix<double> diag,diag_inv;
    int it;
    bool convergence;
    double t1,t2;
    double diffnorm;


    Eigen::DiagonalMatrix<double, Eigen::Dynamic> diagonale;
    int rang;
    int nb_procs;

    /* Initialisation de MPI */
    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rang);
    MPI_Comm_size( MPI_COMM_WORLD, &nb_procs);

    /*lecture du fichier .msh*/
    ifstream FILE;
    FILE.open("./fichierTest/testpart.msh", ios::in);

    if (FILE.fail())
    {
        std::cout<<"Task : "<<rang<< " erreur lors de l'ouverture de FILE"<<std::endl;
        return -1;
    }

    /* Creation des donnees de maillage a partir du fichier lu */
    Maillage mon_maillage=Maillage(FILE);

    /* Calcul des matrices elements finis et des voisinages pour les communciations */
    Probleme mon_probleme(mon_maillage, rang);

    u = *(mon_probleme.Get_u());
    mat_rigidite = *(mon_probleme.Get_p_K());
    second_membre = *(mon_probleme.Get_felim());
    u_exact = *(mon_probleme.Get_uexa());
    diag_inv= *(mon_probleme.Get_diag());

    Eigen::VectorXd second_membre_global(mon_maillage.Get_n_nodes());
    Eigen::MatrixXd K_total(mon_maillage.Get_n_nodes(),mon_maillage.Get_n_nodes());
    Eigen::MatrixXd diag_global(mon_maillage.Get_n_nodes(),mon_maillage.Get_n_nodes());


    if (PARALLELE)
      {
        std::cout<<" on est bien en parallele"<<std::endl;
      MPI_Allreduce(second_membre.data(),second_membre_global.data(),second_membre.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      d_mat_rigidite_local = Eigen::MatrixXd(mat_rigidite);
      MPI_Allreduce(d_mat_rigidite_local.data(),K_total.data(),d_mat_rigidite_local.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      diag_inv_local = Eigen::MatrixXd(diag_inv);
      MPI_Allreduce(diag_inv_local.data(),diag_global.data(),diag_inv_local.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
    else
      {
        second_membre_global = second_membre;
        K_total = mat_rigidite;
        diag_global = diag_inv;
      }

    K_total+=diag_global;

    //Inversion de diag en dehors de la boucle de calcul
    for(int i=0;i<mon_maillage.Get_n_nodes();i++)
    {
        if (!PARALLELE || (mon_probleme.Get_partition_noeud())[i]==rang)
        {
            diag_inv.coeffRef(i,i)=1/diag_inv.coeffRef(i,i);
        }
        else
        {
            diag_inv.coeffRef(i,i)=0;
        }
    }

    

    

    vector<vector<int> > voisins_partition = mon_probleme.Get_voisins_partition();
    vector<vector<int> > voisins_interface = mon_probleme.Get_voisins_interface();


    /* Schema iteratif en temps */
    it = 0;
    convergence = false;

    /* Mesure du temps en seconde dans la boucle en temps */
    t1 = MPI_Wtime();

    while ( !(convergence) && (it < it_max) )
    {   
        it = it+1;
        u_avant_com = u;

        /* Echange des points aux interfaces pour u a l'iteration n */
        if (PARALLELE)
        {
            vector<double> valeurs_a_envoyer;
            if (rang!=0)
            {
                int taille = voisins_interface[rang-1].size();
                for (unsigned int i=0;i<taille;i++)
                {
                   valeurs_a_envoyer.push_back(u.coeffRef(voisins_interface[rang-1][i]-1,0));
                }
                MPI_Send(&valeurs_a_envoyer[0],taille,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            }
            else
            {
                for (int num_proc=1;num_proc<nb_procs;num_proc++)
                {
                    int taille = voisins_interface[num_proc-1].size();
                    vector<double> temp(taille);
                    MPI_Recv(&temp[0],taille,MPI_DOUBLE,num_proc,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    for (int ind=0;ind<taille;ind++)
                    {
                        u.coeffRef(voisins_interface[num_proc-1][ind]-1,0)=temp[ind];
                    }
                }
            }

            if (rang==0)
            {
                for (int num_proc=1;num_proc<nb_procs;num_proc++)
                {
                    int taille = voisins_partition[num_proc-1].size();
                    for(int i=0;i<taille;i++)
                    {
                        valeurs_a_envoyer.push_back(u.coeffRef(voisins_partition[num_proc-1][i]-1,0));
                    }
                    MPI_Send(&valeurs_a_envoyer[0],taille,MPI_DOUBLE,num_proc,0,MPI_COMM_WORLD);
                }
            }
            else
            {
                int taille = voisins_partition[rang-1].size();
                vector<double> temp(taille);
                MPI_Recv(&temp[0],taille,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for (int ind=0;ind<taille;ind++)
                {
                    u.coeffRef(voisins_partition[rang-1][ind]-1,0)=temp[ind];
                }
            }
        }

        /* Calcul de u a l'iteration n+1 */
        u_nouveau = diag_inv * ((-mat_rigidite) * u + second_membre);

        /* Calcul de l'erreur locale */
        double erreur_locale, diffnorm;
        erreur_locale = 0;
        for (int iter=0; iter<u_avant_com.size(); iter++)
        {
            double temp = fabs( u_avant_com[iter] - u_nouveau[iter] );
            if (erreur_locale < temp) 
            {
                erreur_locale = temp;
            }
        }

        /* Calcul de l'erreur sur tous les sous-domaines */
        if (PARALLELE)
        {
            MPI_Allreduce( &erreur_locale, &diffnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
        else
        {
            diffnorm = erreur_locale;
        }

        u = u_nouveau ;

        /* Arret du programme si on a atteint la precision machine obtenu */
        convergence = ( diffnorm < eps );

        /* Affichage pour le processus 0 de la difference */
        if (rang == 0&& it%100==0) 
        {
            cout<<"Iteration "<<it<<" : Erreur_globale = "<<diffnorm<<endl;
        }
    }

    Eigen::VectorXd u_global(mon_maillage.Get_n_nodes());

    if (PARALLELE)
    {
        MPI_Allreduce(u.data(),u_global.data(),u.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    else
    {
        u_global = u;
    }

    /* Mesure du temps a la sortie de la boucle */
    t2 = MPI_Wtime();
    
    /* Calcul de l'erreur par rapport à la solution exacte */
    double erreur_exa=0;
    for (int iter=0;iter<u.size();iter++)
    {
        double temp = fabs( u_global.coeffRef(iter,0) - u_exact.coeffRef(iter,0) );
        if (erreur_exa < temp) 
        {
            erreur_exa = temp;
        }
    }

    if (rang ==  0) 
    {
        cout<<"voici K_total*u_global-second_membre_global qui doit etre nul si l'inversion du probleme lineaire a fonctionné"<<endl;
        affichVector(K_total*u_global-second_membre_global);

        output_vector(u_global,"u_global");
        output_vector(u_exact,"u_exact");
    
        /* Affichage du temps de convergence par le processus 0 */
        cout<<"Convergence apres "<<it<<" iterations en "<< t2-t1<<" secs"<<endl;
        cout<<"L'erreur par rapport a la solution exacte vaut "<<erreur_exa<<endl;

      /* Comparaison de la solution calculee et de la solution exacte
       * sur le processus 0 */
    }

    MPI_Finalize();
    return 0;
}


// Le contenu qui suit sert uniqument à la génération de la page d'index du documentation Doxygen
/*! \mainpage Accueil
 *
 * \section intro_sec Introduction
 *
 * L'objectif de ce projet est la réalisation d'un solveur de système linéaires
 * utilisant la méthode de Jacobi capable de travailler sur un problème à géométrie non strucurée, et ce
 * en utilisant une architecture multiprocesseur.
 *
 * L'outil utilisé pour la génération du problème non structuré est GMSH, qui permet de
 * générer un maillage spatial.
 *
 * Pour la gestion Multiprocesseur, le standard MPI est utilisé, ainsi que la bibliothèque Open MPI (v 1.8.2).
 *  Le compilateur C++ utlilisé est mpic++
 *
 * La résolution du problème posé est faite à l'aide de la méthode de Jacobi par points
 *
 * \section install_sec Contenu
 *
 * \subsection sub1 Fonctionnement du programme 
 *
 *Lorsque le programme se lance, pour chaque tâche MPI, le main.cpp instancie l'objet Probleme correspondant à la
 * tache MPI qui s'execute. 
 *
 * Chaque Probleme instancié initialise ses attributs, dont le Maillage correspondant au problème complet, mais les 
 * matrices de rigidités construites via le parcourt des éléments sont des matrices partielles correspondant aux éléments uniquement liés à la partition étudiée dans cette tache.
 *
 * Une fois toutes les matrices des partition calculées, il est possible de résoudre itérativement le système linéaire en effectuant les communications MPI
 * appropriées selon la partition considérée.
 * 
 * Une fois le critère d'arrêt atteint pour la convergence, le vecteur solution est récupéré, reconstruit, puis sorti dans un fichier texte "u_global".
 *
 *\subsection sub2 Visualisation des résultats
 *
 * Une fois le progamme éxécuté, celui-ci à écrit deux fichiers de sortie, l'un u_exact contient la solution exacte du problème pour la géométrie considérée
 * l'autre u_global contient le résultat de la résolution du système par Jacobi.
 * 
 * A l'aide du script MATLAB Viewing.m ,on peut visualiser facilement les deux solutions obtenues et l'erreur absolue, ce qui donne une
 * idée précise de la validité de la solution proposée pour le problème.
 *
 * ![Rendu du script MATLAB pour un problème à 4 partitions et 1500 points](rendu.png)
 *
 * \subsection sources_sec Sources et fichiers
 *  - C++ source files
 *      -# main.cpp
 *      -# maillage.cpp
 *      -# probleme.cpp
 *      -# nonParallel.cpp
 *
 *  - C++ header files
 *      -# maillage.h
 *      -# probleme.h
 *      -# nonParallel.h
 *
 *  - GMSH mesh files
 *      -# testpart.msh
 *      -# testpart_basic.msh
 *      -# testSimple.msh
 *
 *  - Other files
 *      -# submit_gin.qsub
 *      -# Viewing.m
 *      -# makefile
 */