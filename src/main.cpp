#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "../include/parallel.h"
#include "../include/maillage.h"
#include "../include/probleme.h"
#include "../include/nonParallel.h"
#include "../include/param.h"

int main(int argc, char *argv[])
{
    VectorXd u, u_nouveau, second_membre;
    Eigen::SparseMatrix<double> mat_rigidite;
    Eigen::SparseMatrix<double> diag;
    int it;
    bool convergence;
    double t1,t2;
    double diffnorm;


    Eigen::DiagonalMatrix<double, Eigen::Dynamic> diagonale;
    extern int rang;
    extern int nb_procs;

    /* Initialisation de MPI */
    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rang);
    MPI_Comm_size( MPI_COMM_WORLD, &nb_procs);

    /*lecture du fichier .msh*/
    ifstream FILE;
    FILE.open("./fichierTest/testSimple.msh", ios::in);

    if (FILE.fail())
    {
        std::cout<<"Task : "<<rang<< " erreur lors de l'ouverture de FILE"<<std::endl;
        return -1;
    }

    /* Creation des donnees de maillage a partir du fichier lu */
    Maillage mon_maillage=Maillage(FILE);

    //cout << "Task : "<<rang<< " creation du maillage réussie" << endl;

    /* Calcul des matrices elements finis et des voisinages pour les communciations */
    Probleme mon_probleme=Probleme(mon_maillage, rang);

    //cout << "Task : "<<rang<< " creation du probleme reussie" << endl;

    u = *(mon_probleme.Get_u());
    cout<<"affichage du vecteur solution récupéré depuis probleme :"<<endl;
    affichVector(u);
    second_membre = *(mon_probleme.Get_felim());
    cout<<"afficahge du vecteur second membre récupéré depuis probleme :"<<endl;
    affichVector(second_membre);
    mat_rigidite = *(mon_probleme.Get_p_K());
    diag = *(mon_probleme.Get_diag());
    cout<<"affichage de la matrice de rigidite finale obtenue dans probleme :"<<endl;   
    affich(mat_rigidite);

    //cout << "Task : "<<rang<< " Initialisation des variables" << endl;

    /* Schema iteratif en temps */
    it = 0;
    convergence = false;

    /* Mesure du temps en seconde dans la boucle en temps */
    t1 = MPI_Wtime();
    //cout<<"Task : "<<rang<< " Temps horloge avant le calcul mémorisé t1 = "<<t1<<endl;

    vector<vector<int> > voisins_partition = mon_probleme.Get_voisins_partition();
    vector<vector<int> > voisins_interface = mon_probleme.Get_voisins_interface();

    //cout<<"affichage de voisins_partition :"<<endl;
    //affiche_vector(voisins_partition);
    //cout<<"affichage de voisins_interface :"<<endl;
    //affiche_vector(voisins_interface);

    while ( !(convergence) && (it < it_max) )
    {   
        cout<<"_________________________"<<endl;
        cout<<"Task : "<<rang<< " Iteration "<<it<<endl;
        it = it+1;

        cout<<"Task : "<<rang<< " Echange des valeurs entre interface et partitions"<<endl;
        /* Echange des points aux interfaces pour u a l'iteration n */
        //communication(u, voisins_partition, voisins_interface);
        cout<<"voici le vecteur u dans le proc "<<rang<<" avant communication"<<endl;
        affichVector(u);

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

        cout<<"voici le vecteur u dans le proc "<<rang<<" après communication"<<endl;
        affichVector(u);

        cout << "Task : "<<rang<< " operation de communication terminee" << endl;

        /* Calcul de u a l'iteration n+1 */
        //calcul(u,  u_nouveau, mat_rigidite, diagonale, second_membre);
        for(int i=0;i<mon_maillage.Get_n_nodes();i++)
        {
            if (mon_probleme.Get_partition_noeud()[i]=rang)
            {
                diag.coeffRef(i,i)=1/diag.coeffRef(i,i);
            }

        }
        affich(diag);

        VectorXd vecteur_interm = mat_rigidite * u + second_membre;
        u_nouveau = diag * vecteur_interm;

        /* Calcul de l'erreur globale */
        diffnorm =  erreur_entre_etapes (u, u_nouveau);

        u = u_nouveau ;

        /* Arret du programme si on a atteint la precision machine obtenu */
        convergence = ( diffnorm < eps );

        /* Affichage pour le processus 0 de la difference */
        if ( (rang == 0) && ( (it % 100) == 0) )
        {
            printf("Iteration %d erreur_globale = %g\n", it, diffnorm);
        }
    }

    /* Mesure du temps a la sortie de la boucle */
    t2 = MPI_Wtime();

    if (rang ==  0) {
      /* Affichage du temps de convergence par le processus 0 */
      printf("Convergence apres %d iterations en %f secs\n", it, t2-t1);

      /* Comparaison de la solution calculee et de la solution exacte
       * sur le processus 0 */
    }

    finalisation_mpi();
    return 0;
}


