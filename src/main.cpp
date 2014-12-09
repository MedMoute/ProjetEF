#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include "../include/parallel.h"
#include "../include/calcul.h"

using namespace std;

int main(int argc, int *argv[])
{
    int rang;
    static int nb_procs;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    MPI_Comm_size(MPI_COMM_WOLRD, &nb_procs);

    ifstream FILE;
    FILE.open("./fichierTest/testPart.msh", ios::in);

    if (FILE.fail())
    {
        std::cout<<"ERREUR LORS DE L\'OUVERTURE DU FICHIER MSH"<<std::endl;
        std::cout<<"LE PROGRAMME VA MAINTENANT S\'INTERROMPRE"<<std::endl;
        std::exit (EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    if (rang == 0)
    {
        std::cout<<"Execution code poisson avec "<<nb_procs<<" processus MPI\n"<<std::endl;
    }

     Maillage mon_maillage=Maillage(FILE);
     Probleme mon_probleme=Probleme(mon_maillage);



}


