#include <mpi.h>
#include "../include/parallel.h"
#include "../include/nonParallel.h"
int rang;
int nb_procs;

/*
 * Initialisation pour chaque processus de son rang et du
 * nombre total de processus nb_procs
 */


void communication(VectorXd u, vector<vector<int> > voisins_partition, vector<vector<int> > voisins_interface)
{
    if (rang==0)
    {
    cout << "Task : "<<rang<< " Est entrée dans la fonction de communication de l'interface" << endl;

        const int etiquette = 100;
        MPI_Status statut;
        vector<vector<double> > valeurs_a_envoyer;
        vector<vector<double> > valeurs_a_recevoir;

            //affiche_vector(voisins_partition);
            cout<<u.coeffRef(voisins_partition[1][1])<<endl;

        //valeurs_a_envoyer.resize(nb_procs-1);
        //valeurs_a_recevoir.resize(nb_procs-1);
        cout<<"Task : "<<rang<< " sur "<<nb_procs<<endl;
        for (int i=1;i<nb_procs;i++)
        {
            cout<<"Task : "<<rang<< " On regarde ce qu'on va envoyer au proc "<<i<<" a savoir les "<<voisins_partition[i].size();
            cout<<" valeurs de l'interface qui sont voisines de la partition "<<i<<endl;
            for (unsigned int j=0;j<voisins_partition[i].size();j++)
            {
                cout<<"Test1"<<endl;
                valeurs_a_envoyer[i].push_back(u.coeffRef(voisins_partition[i][j],0));
                cout<<"Test3"<<endl;
            }
            cout<<endl;
            MPI_Send(&valeurs_a_envoyer[i],valeurs_a_envoyer[i].size(),MPI_DOUBLE,i,etiquette,MPI_COMM_WORLD);
            MPI_Recv(&valeurs_a_recevoir[i],voisins_interface[i].size(),MPI_DOUBLE,i,etiquette,MPI_COMM_WORLD,&statut);
        }
    }
    else
    {
    cout << "Task : "<<rang<< " Est entrée dans la fonction de communication de la non-interface" << endl;
    affiche_vector(voisins_interface);
        const int etiquette = 200;
        MPI_Status statut;
        vector<double> valeurs_a_envoyer;
        vector<double> valeurs_a_recevoir;
        cout<<"Task : "<<rang<< " sur "<<nb_procs<<endl;
        for (unsigned int i=0;i<voisins_interface[rang].size();i++)
        {
            valeurs_a_envoyer.push_back(u.coeffRef(voisins_interface[rang][i],0));
        }
        MPI_Send(&valeurs_a_envoyer,valeurs_a_envoyer.size(),MPI_DOUBLE,0,etiquette,MPI_COMM_WORLD);
        MPI_Recv(&valeurs_a_recevoir,voisins_partition[rang].size(),MPI_DOUBLE,0,etiquette,MPI_COMM_WORLD,&statut);
    }
    return ;
}

double erreur_entre_etapes (VectorXd u, VectorXd u_nouveau)
{
    double erreur_locale, diffnorm;

    erreur_locale = 0;
    for (int iter=0; iter<u.size(); iter++)
    {
        double temp = fabs( u[iter] - u_nouveau[iter] );
        if (erreur_locale < temp) erreur_locale = temp;
    }

    /* Calcul de l'erreur sur tous les sous-domaines */
    MPI_Allreduce( &erreur_locale, &diffnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return diffnorm;
}

void finalisation_mpi()
{
  /* Desactivation de MPI */
  MPI_Finalize();
}
