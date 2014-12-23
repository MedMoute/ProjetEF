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
        vector<vector<double> > valeurs_a_envoyer;
        vector<vector<double> > valeurs_a_recevoir;
        valeurs_a_envoyer.resize(nb_procs-1);
        valeurs_a_recevoir.resize(nb_procs-1);
        for (int i=1;i<nb_procs;i++)
        {
            for (unsigned int j=0;j<voisins_partition[i].size();j++)
            {
                valeurs_a_envoyer[i-1].push_back(u.coeffRef(voisins_partition[i-1][j]-1,0));
            }            
            MPI_Send(&valeurs_a_envoyer[i-1][0],valeurs_a_envoyer[i-1].size(),MPI_DOUBLE,i,0,MPI_COMM_WORLD);
            MPI_Recv(&valeurs_a_recevoir[i-1][0],voisins_interface[i-1].size(),MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for (int k=0;k<valeurs_a_recevoir[i-1].size();k++)
            {
                u.coeffRef(voisins_interface[i-1][k],0)=valeurs_a_recevoir[i-1][k];
            }
        }
    }
    else
    {
        cout << "Task : "<<rang<< " Est entrée dans la fonction de communication de la non-interface" << endl;
        affiche_vector(voisins_interface);
        vector<double> valeurs_a_envoyer;
        vector<double> valeurs_a_recevoir;
        for (unsigned int i=0;i<voisins_interface[rang-1].size();i++)
        {
            valeurs_a_envoyer.push_back(u.coeffRef(voisins_interface[rang-1][i]-1,0));
        }
        MPI_Recv(&valeurs_a_recevoir[0],voisins_partition[rang-1].size(),MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(&valeurs_a_envoyer[0],valeurs_a_envoyer.size(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        for (int i=0;i<valeurs_a_recevoir.size();i++)
        {
            u.coeffRef(voisins_partition[rang-1][i],0)=valeurs_a_recevoir[i];
        }
    }
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
