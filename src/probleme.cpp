#include "../include/probleme.h"
#include <stdio.h>
#include <math.h>

Probleme::Probleme(Maillage & monMaillage)
{
    uexa = new VectorXd;
    maillage = &monMaillage;
    g = new VectorXd;
    u = new VectorXd;
    felim = new VectorXd;
    
    p_K = new Eigen::SparseMatrix<double>;
    p_Kelim = new Eigen::SparseMatrix<double>;
    
    int nb_triangles = maillage->n_triangles;
    int *tab_local_global;
    tab_local_global = new int[3*nb_triangles];
    tab_local_global = maillage->triangles_sommets;
    
    double *coorneu = maillage->nodes_coords;

    p_K->resize(maillage->n_nodes,maillage->n_nodes);
    
    for(int ind_triangle=0;ind_triangle<nb_triangles;ind_triangle++)
    {
        double *p_K_elem;
        p_K_elem = new double[9];

        mat_K_elem(p_K_elem,coorneu,tab_local_global,ind_triangle);
        assemblage(*p_K,p_K_elem,tab_local_global,ind_triangle);
        delete p_K_elem;
        p_K_elem=0;
        
    }
    
    
}

Probleme::~Probleme()
{
    delete p_K;
    delete uexa ;
    delete g ;
    delete u ;
    delete felim ;
}

void Probleme::affich(Eigen::SparseMatrix<double> _mat)
{

    std::cout<<std::endl<<"Affichage de la matrice de rigidite :";
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

void Probleme::assemblage(Eigen::SparseMatrix<double> & mat, double* mat_elem, int* tab,int n){

    for(int i=0;i<3;i++)
    {
        int ind_global_1=tab[3*n+i];

        for(int j=0;j<3;j++)
        {
            int ind_global_2=tab[3*n+j];
            //Debug : Affiche les sommets correspondant à la matrice élémentaire en construction
            //std::cout<<" Sommets : "<<ind_global_1<<" | " <<ind_global_2<< std::endl;
            double AddedCoeff = mat_elem[3*i+j];

            mat.coeffRef(ind_global_1-1,ind_global_2-1)+=AddedCoeff;

        }
    }
}

void Probleme::mat_K_elem(double *mat_elem,double *coorneu, int *tab_local_global, int ind_triangle){

    int ind_pt1=tab_local_global[3*ind_triangle];
    int ind_pt2=tab_local_global[3*ind_triangle+1];
    int ind_pt3=tab_local_global[3*ind_triangle+2];

    // Les -3 sont dus au fait que les valeurs commencent à l'indice 0 or min(ind_pt)=1
    double x1=coorneu[0+ind_pt1*3-3];
    double y1=coorneu[1+ind_pt1*3-3];
    double z1=coorneu[2+ind_pt1*3-3];
    double x2=coorneu[0+ind_pt2*3-3];
    double y2=coorneu[1+ind_pt2*3-3];
    double z2=coorneu[2+ind_pt2*3-3];
    double x3=coorneu[0+ind_pt3*3-3];
    double y3=coorneu[1+ind_pt3*3-3];
    double z3=coorneu[2+ind_pt3*3-3];

    double x12 = x2 - x1;
    double y12 = y2 - y1;
    double z12 = z2 - z1;
    double x13 = x3 - x1;
    double y13 = y3 - y1;
    double z13 = z3 - z1;
    double x23 = x3 - x2;
    double y23 = y3 - y2;
    double z23 = z3 - z2;

    double tab_interm[] = {-y23,x23,y13,-x13,-y12,x12};

    // Le calcul de la méthode qui suit semble être faux...
    //double aire_triangle = (1/2)*sqrt((y12*z13 - z12*y13)*(y12*z13 - z12*y13)+((z12*x13-x12*z13)*(z12*x13-x12*z13))+((x12*y13-y12*x13)*(x12*y13-y12*x13)));

    //Formule de héron
    double a=sqrt((x12)*(x12)+(y12)*(y12)+(z12)*(z12));
    double b=sqrt((x13)*(x13)+(y13)*(y13)+(z13)*(z13));
    double c=sqrt((x23)*(x23)+(y23)*(y23)+(z23)*(z23));

    double s=(a+b+c)/2;

    double aire_triangle = sqrt(s*(s-a)*(s-b)*(s-c));

    //std::cout<<std::endl<< "Aire du triangle : "<<aire_triangle<<std::endl;

    //rajouter une erreur si l'aire du triangle est nulle
    //A partir de maintenant on considere que la troisieme coord est toujours nulle

    //std::cout <<std::endl<< "On affiche la matrice elementaire"<<std::endl;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            mat_elem[3*i+j]=(1/(4*aire_triangle))*(tab_interm[2*i]*tab_interm[2*j]+tab_interm[2*i+1]*tab_interm[2*j+1]);
            //   std::cout<<mat_elem[3*i+j]<<" | ";
        }
        // std::cout<<std::endl;
    }
}



