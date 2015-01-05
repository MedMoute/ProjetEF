/*!
 * \file nonParallel.cpp
 * \brief Source des fonctions n'exploitant pas la paraléllisation dans le solveur
 * \version 1.0
 */
#include <fstream>
#include <iostream> 
#include "../include/nonParallel.h"


void affichVector(VectorXd V)
{

    for(int i=0;i<V.rows();i++)
    {
        std::cout.precision(2);
        std::cout<<V(i,0)<<"|";
    }
    std::cout<<"__________"<<std::endl;
}
void affich(Eigen::SparseMatrix<double> _mat)
{
    std::cout<<std::endl<<"-------------------------------------"<<std::endl;

    if (_mat.size()!=0)
    {
        Eigen::MatrixXd dMat;

        dMat=Eigen::MatrixXd(_mat);

        for(int i=0;i<dMat.rows();i++)
        {
            for(int j=0;j<dMat.cols();j++)
            {
                std::cout.precision(2);

                std::cout<<dMat(i,j)<<"|";
            }
            std::cout<<std::endl;
        }
    }


    else
    {
        std::cout<<"Matrice de taille zero."<<std::endl;
    }
}

void affiche_vector(vector<vector<int> > v)
{
    for (int int_ligne=0;int_ligne<v.size();int_ligne++)
    {
        for (int int_col=0;int_col<v[int_ligne].size();int_col++)
        {
            std::cout<<v[int_ligne][int_col]<<" | ";
        }
        std::cout<<std::endl;
    }
}


void output_vector(VectorXd u,string s)
{
    ofstream OuputFile;
    OuputFile.open (s.c_str());

    for(int i=0;i<u.rows();i++)
    {
    OuputFile << u(i,0);
    OuputFile << "\n";
    }
    OuputFile.close();
}