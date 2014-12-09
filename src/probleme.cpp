
#include "../include/probleme.h"
#include <stdio.h>
#include <math.h>
#define PI 3.14159

Probleme::Probleme(Maillage & monMaillage)
{
    uexa = new VectorXd;
    maillage = &monMaillage;
    g = new VectorXd;
    u = new VectorXd;
    felim = new VectorXd;

    p_K = new Eigen::SparseMatrix<double> (maillage->n_nodes,maillage->n_nodes);
    p_Kelim = new Eigen::SparseMatrix<double> (maillage->n_nodes,maillage->n_nodes);
/*
    //vecteur qui a en position i la partition a laquelle appartient le noeud i
    int* partition_noeud = new int[maillage->n_nodes];

    //vecteur

    for (ind_triangle=0;ind_triangle<maillage->n_triangles;ind_triangle++)
    {
        int ind_pt1=maillage->triangles_sommets[3*ind_triangle];
        int ind_pt2=maillage->triangles_sommets[3*ind_triangle+1];
        int ind_pt3=maillage->triangles_sommets[3*ind_triangle+2];

        numero_partition_triangle=maillage->partition_ref[ind_triangle];

        if (partition_noeud[ind_pt1]=0)
        {
            partition_noeud[ind_pt1]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt1]=numero_partition_triangle)
        {
        }
        else
        {
            partition_noeud[ind_pt1]=3;
        }

        if (partition_noeud[ind_pt2]=0)
        {
            partition_noeud[ind_pt2]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt2]=numero_partition_triangle)
        {
        }
        else
        {
            partition_noeud[ind_pt2]=3;
        }

        if (partition_noeud[ind_pt3]=0)
        {
            partition_noeud[ind_pt3]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt3]=numero_partition_triangle)
        {
        }
        else
        {
            partition_noeud[ind_pt3]=3;
        }

    }

    voisins_interface = new int[maillage->]

    uexa->resize(maillage->n_nodes,1);

    for(int ind_node=0;ind_node<maillage->n_nodes;ind_node++){
        uexa->coeffRef(ind_node,0)+=calcul_uexa(maillage->nodes_coords[3*ind_node],maillage->nodes_coords[3*ind_node+1]);
    }
    */
    g->resize(maillage->n_nodes,1);

    //on initialise les conditions au bord, dans un premier temps a une constante

    for(int ind_node=0;ind_node<maillage->n_nodes;ind_node++){
        if (maillage->nodes_ref[ind_node]!=0){
            double x=maillage->nodes_coords[3*ind_node+0];
            double y=maillage->nodes_coords[3*ind_node+1];
            double addedCoeff = calcul_g(x,y);
            g->coeffRef(ind_node,0)+=addedCoeff;
        }
    }
    std::cout<<"vecteur conditions aux limites"<<endl;
    affichVector(*g);

    felim->resize(maillage->n_nodes,1);

    for(int ind_triangle=0;ind_triangle<maillage->n_triangles;ind_triangle++){
        int ind_pt1=maillage->triangles_sommets[3*ind_triangle];
        int ind_pt2=maillage->triangles_sommets[3*ind_triangle+1];
        int ind_pt3=maillage->triangles_sommets[3*ind_triangle+2];

        // Les -3 sont dus au fait que les valeurs commencent à l'indice 0 or min(ind_pt)=1
        double x1=maillage->nodes_coords[0+ind_pt1*3-3];
        double y1=maillage->nodes_coords[1+ind_pt1*3-3];
        double z1=maillage->nodes_coords[2+ind_pt1*3-3];
        double x2=maillage->nodes_coords[0+ind_pt2*3-3];
        double y2=maillage->nodes_coords[1+ind_pt2*3-3];
        double z2=maillage->nodes_coords[2+ind_pt2*3-3];
        double x3=maillage->nodes_coords[0+ind_pt3*3-3];
        double y3=maillage->nodes_coords[1+ind_pt3*3-3];
        double z3=maillage->nodes_coords[2+ind_pt3*3-3];

        //calcul des composantes des vecteurs formant les aretes du triangle
        double x12 = x2 - x1;
        double y12 = y2 - y1;
        double z12 = z2 - z1;
        double x13 = x3 - x1;
        double y13 = y3 - y1;
        double z13 = z3 - z1;
        double x23 = x3 - x2;
        double y23 = y3 - y2;
        double z23 = z3 - z2;

        //Calcul de l'aire du triangle avec la formule de héron

        double a=sqrt((x12)*(x12)+(y12)*(y12)+(z12)*(z12));
        double b=sqrt((x13)*(x13)+(y13)*(y13)+(z13)*(z13));
        double c=sqrt((x23)*(x23)+(y23)*(y23)+(z23)*(z23));

        double s=(a+b+c)/2;

        double aire_triangle = sqrt(s*(s-a)*(s-b)*(s-c));

        //Assemblage du second membre F par la formule de quadrature Gauss Lobatto

        //On definit les cinq coordonnees utiles dans la formule

        double s_0=1./3;
        double s_1=(6-sqrt(15))/21;
        double s_2=(6+sqrt(15))/21;
        double s_3=(9+2*sqrt(15))/21;
        double s_4=(9-2*sqrt(15))/21;

        //et les trois poids

        double eta_0=9./80;
        double eta_1=(155-sqrt(15))/2400;
        double eta_2=(155+sqrt(15))/2400;

        //chaque sommet du triangle de la boucle a trois contributions au vecteur felim

        for(int j=0;j<3;j++){
            //la liste des poids et points pour l'interpolation
            double tab[] = {eta_0,s_0,s_0,eta_1,s_1,s_1,eta_1,s_1,s_3,eta_1,s_3,s_1,eta_2,s_2,s_2,eta_2,s_2,s_4,eta_2,s_4,s_2};
            double integ_interpolee=0;
            for (int i=0;i<7;i++){
                //l'integrale sur le triangle de la boucle se ramene au triangle de base par chgmt de var.
                //(s1,s2) est le nouveau point ou l'on calcule f apres changement de var
                double s1=x12*tab[3*i+1]+x13*tab[3*i+2]+x1;
                double s2=y12*tab[3*i+1]+y13*tab[3*i+2]+y1;
                //base_loc(j,x,y) calcule la valeur de la fonction barycentrique lambda_j en (x,y)
                double term_som=tab[3*i]*base_loc(j,tab[3*i+1],tab[3*i+2])*calcul_f(s1,s2);
                integ_interpolee+=term_som;
            }

            double addedCoeff = 2*aire_triangle*integ_interpolee;

            //l'assemblage a proprement parler

            int ind_global=maillage->triangles_sommets[3*ind_triangle+j];
            felim->coeffRef(ind_global-1,0)+=addedCoeff;
        }

        double tab_interm[] = {-y23,x23,y13,-x13,-y12,x12};

        //creation de la matrice elementaire pour le triangle de la boucle
        double *p_K_elem;
        p_K_elem = new double[9];

        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                //remplissage de la matrice elementaire par le calcul de l integrale des fonctions barycentriques
                p_K_elem[3*i+j]=(1/(4*aire_triangle))*(tab_interm[2*i]*tab_interm[2*j]+tab_interm[2*i+1]*tab_interm[2*j+1]);
            }
        }
        //l'assemblage de p_K a proprement parler
        for(int i=0;i<3;i++){
            int ind_global_1=maillage->triangles_sommets[3*ind_triangle+i];

            for(int j=0;j<3;j++){
                int ind_global_2=maillage->triangles_sommets[3*ind_triangle+j];
                //Debug : Affiche les sommets correspondant à la matrice élémentaire en construction
                //std::cout<<" Sommets : "<<ind_global_1<<" | " <<ind_global_2<< std::endl;
                double AddedCoeff = p_K_elem[3*i+j];

                p_K->coeffRef(ind_global_1-1,ind_global_2-1)+=AddedCoeff;

            }
        }

        //Je mets en commentaire cette partie ou on decide si on fait en sequentiel ou pas
       /* if (maillage->nb_partitions==2)
        {
            //Résolution séquentielle et assemblage "simple"
            mat_K_elem(p_K_elem,coorneu,tab_local_global,ind_triangle);
            assemblage(*p_K,p_K_elem,tab_local_global,ind_triangle);
        }
        else
        {
            mat_K_elem_par(p_K_elem,coorneu,tab_local_global,ind_triangle);
            assemblage_par(*p_K,p_K_elem,tab_local_global,ind_triangle);
        }
        */

        delete p_K_elem;
        p_K_elem=0;
    }

    //std::cout<<"second membre avant pseudo elimination"<<endl;
    //affichVector(*felim);



    //affich(*p_K);

    Eigen::SparseMatrix<double> K_err = *p_K;

    //le second membre total, prenant en compte les f et g de la formulation variationnelle
    //felim a ete obtenu par formules de quadrature, le produit p_K*g montre que le second membre en g est obtenu par interpolation
    *felim=*felim-(*p_K) * (*g);
    std::cout<<"second membre dans la reso du systeme avant pseudi elimin"<<endl;
    affichVector(*felim);
    //algorithme de pseudo elimination
    for(int i=0;i<maillage->n_nodes;i++)
    {
        if (maillage->nodes_ref[i]==1)
        {
            //std::cout<<"le noeud  "<<i<<"est sur le bord"<<endl;
            //seuls les elements du bord sur felim doivent etre changes
            double nouvCoeff = p_K->coeffRef(i,i)*g->coeffRef(i,0);
            felim->coeffRef(i,0)=nouvCoeff;
            for(int j=0;j<maillage->n_nodes;j++)
            {
                if (j!=i)
                {
                    //on met a zero les lignes et colonne (hors diagonale) dont les indices sont sur le bord
                    //affichUnElem(*p_K,i,j);
                    p_K->coeffRef(i,j)=0;
                    p_K->coeffRef(j,i)=0;
                }
            }
        }
    }

    affich(*p_K);
    //std::cout<<"second membre dans la resolution du systeme lineaire"<<endl;
    //affichVector(*felim);

    /* Schema iteratif en temps */

    /*
    int it = 0;
    convergence = faux;


    while ( !(convergence) && (it < it_max) ) {
        it = it+1;

        temp = u; u = u_nouveau; u_nouveau = temp;

        communication( u );

        calcul( u,  u_nouveau);

        diffnorm =  erreur_globale (u, u_nouveau);

        convergence = ( diffnorm < eps );

        if ( (rang == 0) && ( (it % 100) == 0) )
          printf("Iteration %d erreur_globale = %g\n", it, diffnorm);

    std::cout<<"la solution finale"<<endl;
    affichVector(*u);

    std::cout<<"l'erreur"<<endl;
    affichVector(*u-*uexa);

    double erreurH1=((*u-*uexa).dot(K_err*((*u-*uexa))));
    cout<<"l'erreur H1 vaut "<<erreurH1<<endl;
    */
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

/*void Probleme::assemblage_par(Eigen::SparseMatrix<double> & mat, double* mat_elem, int* tab,int n)
{
    /*Version parallélisée de l'assemblage de la matrice

     * On suppose avoir N+1 coeurs de calculs pour N partitions de tailles N_i noeuds chacunes
     *(le N-ième coeur correspond au noeuds de l'interface entre partitions)

     * La matrice de rigidité est globalement creuse, pour une partition, on peut définir
     *une matrice de rigidité locale qui contiendra :
     *                   -La matrice des noeuds internes à la partition pour les N premieres coeurs de travail
     * avec size(P_K_loc(i))=N_i*N_i
     *                   -La matrice des noeuds présents sur l'interface pour le N+1 ème coeur
     * avec size(P_K_loc(i))=N*N
     *
     *A partir d'un objet de la classe maillage créé par la méthode Maillage::Maillage(ifstream& fd)
     *chaque process va créer sa matrice de rigidité locale.


}


void Probleme::mat_K_elem_par(double *mat_elem,double *coorneu, int *tab_local_global, int ind_triangle)
{
    //TODO Aussi
}
*/
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
