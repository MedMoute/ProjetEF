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
    extern int rang;
    extern int nb_procs;

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
    
    second_membre = *(mon_probleme.Get_felim());
    Eigen::VectorXd second_membre_global(mon_maillage.Get_n_nodes());
    MPI_Allreduce(second_membre.data(),second_membre_global.data(),second_membre.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    mat_rigidite = *(mon_probleme.Get_p_K());
    d_mat_rigidite_local = Eigen::MatrixXd(mat_rigidite);
    Eigen::MatrixXd K_total(mon_maillage.Get_n_nodes(),mon_maillage.Get_n_nodes());
    MPI_Allreduce(d_mat_rigidite_local.data(),K_total.data(),d_mat_rigidite_local.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    u_exact = *(mon_probleme.Get_uexa());

    diag_inv= *(mon_probleme.Get_diag());
    diag_inv_local = Eigen::MatrixXd(diag_inv);
    Eigen::MatrixXd diag_global(mon_maillage.Get_n_nodes(),mon_maillage.Get_n_nodes());
    MPI_Allreduce(diag_inv_local.data(),diag_global.data(),diag_inv_local.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    K_total+=diag_global;

    //Inversion de diag en dehors de la boucle de calcul
    for(int i=0;i<mon_maillage.Get_n_nodes();i++)
    {
        if ((mon_probleme.Get_partition_noeud())[i]==rang)
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
        if (rang == 0) 
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
    
        /* Affichage du temps de convergence par le processus 0 */
        cout<<"Convergence apres "<<it<<" iterations en "<< t2-t1<<" secs"<<endl;
        cout<<"L'erreur par rapport a la solution exacte vaut "<<erreur_exa<<endl;

      /* Comparaison de la solution calculee et de la solution exacte
       * sur le processus 0 */
    }

    MPI_Finalize();
    return 0;
}


