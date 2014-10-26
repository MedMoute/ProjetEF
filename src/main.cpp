#include "../include/probleme.h"

int main(int argc, char* argv[])
{
    ifstream FILE;
    FILE.open("./fichiertest/testpart.msh", ios::in);

    if (FILE.fail())
    {
        std::cout<<"ERREUR LORS DE L\'OUVERTURE DU FICHIER MSH"<<std::endl;
        std::cout<<"LE PROGRAMME VA MAINTENANT S\'INTERROMPRE"<<std::endl;

        std::exit (EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    Maillage mon_maillage=Maillage(FILE);
    Probleme mon_probleme=Probleme(mon_maillage);

    mon_probleme.affich(*(mon_probleme.p_K));
    mon_probleme.~Probleme();
    return EXIT_SUCCESS;


}



