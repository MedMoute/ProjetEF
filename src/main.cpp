#include <iostream>
#include <fstream>
#include <mpi.h>
#include "../include/parallel.h"
#include "../include/maillage.h"
#include "../include/probleme.h"
#include "../include/nonParallel.h"
#include "../include/param.h"

int main(int argc, char *argv[])
{
    VectorXd u, u_nouveau, second_membre;
    int it,convergence;
    double t1,t2;
    double diffnorm;
    Eigen::SparseMatrix<double> mat_rigidite;

    /* Initialisation de MPI */
    initialisation_mpi( argc, argv);

    /*lecture du fichier .msh*/
    ifstream FILE;
    FILE.open("./fichierTest/testPart.msh", ios::in);

    if (FILE.fail())
    {
        std::cout<<"erreur lors de l'ouverture de FILE"<<std::endl;
        return -1;
    }

    /* Creation des donnees de maillage a partir du fichier lu */
    Maillage mon_maillage=Maillage(FILE);

    /* Calcul des matrices elements finis et des voisinages pour les communciations */
    Probleme mon_probleme=Probleme(mon_maillage, rang);

    u = *(mon_probleme.u);
    mat_rigidite = *(mon_probleme.p_Kelim);
    second_membre = *(mon_probleme.felim);

    /* Schema iteratif en temps */
    it = 0;
    convergence = faux;

    /* Mesure du temps en seconde dans la boucle en temps */
    t1 = MPI_Wtime();

    while ( !(convergence) && (it < it_max) )
    {
        it = it+1;

        /* Echange des points aux interfaces pour u a l'iteration n */
        communication(u, mon_probleme.voisins_partition, mon_probleme.voisins_interface);

        /* Calcul de u a l'iteration n+1 */
        calcul(u,  u_nouveau, mat_rigidite, second_membre);

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


