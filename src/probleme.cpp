
#include "../include/probleme.h"
#include "../include/parallel.h"
#include <stdio.h>
#include <math.h>
#define PI 3.14159
#define PARALLELE true

Probleme::Probleme(Maillage & monMaillage, int rang)
{
    uexa = new VectorXd;
    uexa->resize(maillage->n_nodes,1);

    maillage = &monMaillage;

    g = new VectorXd;
    g->resize(maillage->n_nodes,1);

    u = new VectorXd;

    felim = new VectorXd;
    felim->resize(maillage->n_nodes,1);

    p_K = new Eigen::SparseMatrix<double> (maillage->n_nodes,maillage->n_nodes);
    p_Kelim = new Eigen::SparseMatrix<double> (maillage->n_nodes,maillage->n_nodes);

    partition_noeud = new int[maillage->n_nodes];

    vector<vector<int> > voisins_interface;
    voisins_interface.resize(maillage->nb_partitions-1);
    vector<vector<int> > voisins_partition;
    voisins_partition.resize(maillage->nb_partitions-1);

    calcul_voisins();

    /* Calcul de la solution exacte */

    for(int ind_node=0;ind_node<maillage->n_nodes;ind_node++)
    {
        uexa->coeffRef(ind_node,0)+=calcul_uexa(maillage->nodes_coords[3*ind_node],maillage->nodes_coords[3*ind_node+1]);
    }

    /* Initialisation des conditions au bord, dans un premier temps a une constante*/

    for(int ind_node=0;ind_node<maillage->n_nodes;ind_node++)
    {
        if (maillage->nodes_ref[ind_node]!=0)
        {
            double x=maillage->nodes_coords[3*ind_node+0];
            double y=maillage->nodes_coords[3*ind_node+1];
            double addedCoeff = calcul_g(x,y);
            g->coeffRef(ind_node,0)+=addedCoeff;
        }
    }

    /* Assemblage de la matrice de rigiditÃ© par parcours de tous les triangles */

    assemblage(rang);

    /* On garde en mÃ©moire la matrice assemblÃ©e avant pseudo Ã©limination pour calculer l'erreur H1 plus tard */

    Eigen::SparseMatrix<double> K_err = *p_K;

    /* Le second membre total, prenant en compte les f et g de la formulation variationnelle
     * felim a ete obtenu par formules de quadrature. Le second membre en g est obtenu par interpolation
     */

    *felim=*felim-(*p_K) * (*g);

    /* Mise en oeuvre de la pseudo elimination */

    for(int i=0;i<maillage->n_nodes;i++)
    {
        /* seuls les elements du bord sur felim doivent etre changes */

        if (maillage->nodes_ref[i]==1)
        {
            double nouvCoeff = p_K->coeffRef(i,i)*g->coeffRef(i,0);
            felim->coeffRef(i,0)=nouvCoeff;
            for(int j=0;j<maillage->n_nodes;j++)
            {
                if (j!=i)
                {
                    /* Mise a zero des lignes et colonne (hors diagonale) dont les indices sont sur le bord */

                    p_K->coeffRef(i,j)=0;
                    p_K->coeffRef(j,i)=0;
                }
            }
        }
    }
}


void Probleme::calcul_voisins()
{
    for (int ind_triangle=0;ind_triangle<maillage->n_triangles;ind_triangle++)
    {
        int ind_pt1=maillage->triangles_sommets[3*ind_triangle];
        int ind_pt2=maillage->triangles_sommets[3*ind_triangle+1];
        int ind_pt3=maillage->triangles_sommets[3*ind_triangle+2];

        int numero_partition_triangle=maillage->partition_ref[ind_triangle];

        /* Un noeud Ã  l'interface est reconnu par son appartenance Ã  deux partitions distinctes. */

        if (partition_noeud[ind_pt1]==0)
        {
            partition_noeud[ind_pt1]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt1]==numero_partition_triangle)
        {
        }
        else
        {
            partition_noeud[ind_pt1]=0;
        }

        if (partition_noeud[ind_pt2]==0)
        {
            partition_noeud[ind_pt2]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt2]==numero_partition_triangle)
        {
        }
        else
        {
            partition_noeud[ind_pt2]=0;
        }

        if (partition_noeud[ind_pt3]==0)
        {
            partition_noeud[ind_pt3]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt3]==numero_partition_triangle)
        {
        }
        else
        {
            partition_noeud[ind_pt3]=0;
        }
    }

    /* On sait maintenant quel noeud appartient Ã  quel triangle. En parcourant Ã  nouveau les triangles,
     * en s'arretant sur ceux qui sont au contact de l'interface sans y etre totalement inclus, on peut
     * donc dÃ©terminer les vecteurs voisins_interface et voisins_partition
     */

    for (int ind_triangle=0;ind_triangle<maillage->n_triangles;ind_triangle++)
    {
        int ind_pt1=maillage->triangles_sommets[3*ind_triangle];
        int ind_pt2=maillage->triangles_sommets[3*ind_triangle+1];
        int ind_pt3=maillage->triangles_sommets[3*ind_triangle+2];

        int numero_partition_triangle=maillage->partition_ref[ind_triangle];

        if ((partition_noeud[ind_pt1]==0 || partition_noeud[ind_pt2]==0 || partition_noeud[ind_pt3]==0) &&
                (!(partition_noeud[ind_pt1]==0 && partition_noeud[ind_pt2]==0 && partition_noeud[ind_pt3]==0)))
        {
            if (partition_noeud[ind_pt2]!=0)
            {
                /* Si le noeud n'est pas sur l'interface, il en est voisin */
                voisins_interface[partition_noeud[ind_pt2]].push_back(ind_pt2);
            }
            else
            {
                /* Si le noeud est sur l'interface, il est voisin de la partition du triangle auquel il appartient */
                voisins_partition[numero_partition_triangle].push_back(ind_pt2);
            }
            if (partition_noeud[ind_pt3]!=0)
            {
                voisins_interface[partition_noeud[ind_pt3]].push_back(ind_pt3);
            }
            else
            {
                voisins_partition[numero_partition_triangle].push_back(ind_pt3);
            }
            if (partition_noeud[ind_pt1]!=0)
            {
                voisins_interface[partition_noeud[ind_pt1]].push_back(ind_pt1);
            }
            else
            {
                voisins_partition[numero_partition_triangle].push_back(ind_pt1);
            }
        }
     }

    /* On trie chaque liste de noeuds des vecteurs voisins_partition et voisins_interface par ordre
     * croissant afin de pouvoir les utiliser dans les communications
     * */

    for (int i=0;i<maillage->nb_partitions-1;i++)
    {
        std::sort (voisins_partition[i].begin(),voisins_partition[i].end());
        std::sort (voisins_interface[i].begin(),voisins_interface[i].end());
    }

}




void Probleme::assemblage(int rang)
{
    /* On definit les cinq coordonnees utiles dans la formule */

    double s_0=1./3;
    double s_1=(6-sqrt(15))/21;
    double s_2=(6+sqrt(15))/21;
    double s_3=(9+2*sqrt(15))/21;
    double s_4=(9-2*sqrt(15))/21;

    /* et les trois poids */

    double eta_0=9./80;
    double eta_1=(155-sqrt(15))/2400;
    double eta_2=(155+sqrt(15))/2400;

    /* puis le vecteur pour boucler */

    double tab[] = {eta_0,s_0,s_0,eta_1,s_1,s_1,eta_1,s_1,s_3,eta_1,s_3,s_1,eta_2,s_2,s_2,eta_2,s_2,s_4,eta_2,s_4,s_2};

    for(int ind_triangle=0;ind_triangle<maillage->n_triangles;ind_triangle++)
    {
        int ind_pt1=maillage->triangles_sommets[3*ind_triangle];
        int ind_pt2=maillage->triangles_sommets[3*ind_triangle+1];
        int ind_pt3=maillage->triangles_sommets[3*ind_triangle+2];

        /* Les -3 sont dus au fait que les valeurs commencent Ã  l'indice 0 or min(ind_pt)=1 */

        double x1=maillage->nodes_coords[0+ind_pt1*3-3];
        double y1=maillage->nodes_coords[1+ind_pt1*3-3];
        double z1=maillage->nodes_coords[2+ind_pt1*3-3];
        double x2=maillage->nodes_coords[0+ind_pt2*3-3];
        double y2=maillage->nodes_coords[1+ind_pt2*3-3];
        double z2=maillage->nodes_coords[2+ind_pt2*3-3];
        double x3=maillage->nodes_coords[0+ind_pt3*3-3];
        double y3=maillage->nodes_coords[1+ind_pt3*3-3];
        double z3=maillage->nodes_coords[2+ind_pt3*3-3];

        /* Calcul des composantes des vecteurs formant les aretes du triangle */

        double x12 = x2 - x1;
        double y12 = y2 - y1;
        double z12 = z2 - z1;
        double x13 = x3 - x1;
        double y13 = y3 - y1;
        double z13 = z3 - z1;
        double x23 = x3 - x2;
        double y23 = y3 - y2;
        double z23 = z3 - z2;

        /* Calcul de l'aire du triangle avec la formule de hÃ©ron */

        double a=sqrt((x12)*(x12)+(y12)*(y12)+(z12)*(z12));
        double b=sqrt((x13)*(x13)+(y13)*(y13)+(z13)*(z13));
        double c=sqrt((x23)*(x23)+(y23)*(y23)+(z23)*(z23));

        double s=(a+b+c)/2;
        double aire_triangle = sqrt(s*(s-a)*(s-b)*(s-c));

        /* Assemblage de felim */

        assemblage_felim(tab, x12, x13, x1, y12, y13, y1, aire_triangle, ind_triangle, rang);

        /* Calcul des matrices Ã©lÃ©mentaires */

        /* Creation de la matrice elementaire pour le triangle de la boucle */

        double *p_K_elem;
        p_K_elem = new double[9];

        double tab_interm[] = {-y23,x23,y13,-x13,-y12,x12};
        mat_K_elem(tab_interm, p_K_elem, aire_triangle);

        /* Assemblage de p_K Ã  partir de chaque matrice Ã©lÃ©mentaire */
        assemblage_pKelem(ind_triangle, p_K_elem, rang);

        /* Suppression des matrices Ã©lÃ©mentaires provisoires */

        delete p_K_elem;
        p_K_elem=0;
    }
}

void Probleme::assemblage_felim(double* tab, double x12, double x13, double x1, double y12, double y13, double y1, double
                                aire_triangle, int ind_triangle, int rang)
{
    double integ_interpolee=0;
    for(int j=0;j<3;j++)
    {
        int ind_global=maillage->triangles_sommets[3*ind_triangle+j];
        if (PARALLELE && partition_noeud[ind_global]!=rang)
        {
        }
        else
        {
            for (int i=0;i<7;i++)
            {
                /* L'integrale sur le triangle de la boucle se ramene au triangle de base par chgmt de var.
                 * (s1,s2) est le nouveau point ou l'on calcule f apres changement de var */

                double s1=x12*tab[3*i+1]+x13*tab[3*i+2]+x1;
                double s2=y12*tab[3*i+1]+y13*tab[3*i+2]+y1;
                /* base_loc(j,x,y) calcule la valeur de la fonction barycentrique lambda_j en (x,y) */
                double term_som=tab[3*i]*base_loc(j,tab[3*i+1],tab[3*i+2])*calcul_f(s1,s2);
                integ_interpolee+=term_som;
            }
            double addedCoeff = 2*aire_triangle*integ_interpolee;
            /* Assemblage de felim */
            felim->coeffRef(ind_global-1,0)+=addedCoeff;
        }
    }
}

void Probleme::mat_K_elem(double* tab_interm, double* p_K_elem, double aire_triangle)
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            /* Remplissage de la matrice elementaire par le calcul de l integrale des fonctions barycentriques */
            p_K_elem[3*i+j]=(1/(4*aire_triangle))*(tab_interm[2*i]*tab_interm[2*j]+tab_interm[2*i+1]*tab_interm[2*j+1]);
        }
    }
}

void Probleme::assemblage_pKelem(int ind_triangle, double* p_K_elem, int rang)
{
    for(int i=0;i<3;i++)
    {
        int ind_global_1=maillage->triangles_sommets[3*ind_triangle+i];
        for(int j=0;j<3;j++)
        {
            int ind_global_2=maillage->triangles_sommets[3*ind_triangle+j];
            if (PARALLELE && partition_noeud[ind_global_1]!=rang && partition_noeud[ind_global_2]!=rang)
            {
            }
            else
            {
                double AddedCoeff = p_K_elem[3*i+j];
                p_K->coeffRef(ind_global_1-1,ind_global_2-1)+=AddedCoeff;
            }
        }
    }
}

void Probleme::affichVector(VectorXd V)
{

    for(int i=0;i<V.rows();i++)
    {
        std::cout.precision(2);
        std::cout<<V(i,0);
        std::cout<<std::endl;
    }
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


double Probleme::base_loc(int j, double coor_1, double coor_2)
{
    return((j==0)*(1-coor_1-coor_2)+(j==1)*(coor_1)+(j==2)*(coor_2));
}

double Probleme::calcul_f(double coor_1, double coor_2)
{
    return(2*pow(PI,2)*sin(PI*coor_1)*sin(PI*coor_2));
}

double Probleme::calcul_uexa(double coor_1, double coor_2)
{
    return(sin(PI*coor_1)*sin(coor_2*PI));
}

double Probleme::calcul_g(double coor_1, double coor_2)
{
    return(0*coor_1+0*coor_2);
}

Probleme::~Probleme()
{
    delete uexa;
    uexa=0;
    delete g;
    g=0;
    delete felim;
    felim=0;
    delete p_K;
    p_K=0;
    delete p_Kelim;
    p_Kelim=0;
    delete u;
    u=0;
}

