
#include "../include/probleme.h"
#include "../include/parallel.h"
#include "../include/nonParallel.h"
#include <stdio.h>
#include <math.h>
#define PI 3.14159
#define PARALLELE true

Probleme::Probleme(Maillage monMaillage, int rang)
{
    Set_maillage(&monMaillage);

    uexa = new VectorXd;
    uexa->resize(maillage->Get_n_nodes(),1);

    g = new VectorXd;
    g->resize(maillage->Get_n_nodes(),1);

    u = new VectorXd;
    u->resize(maillage->Get_n_nodes(),1);

    felim = new VectorXd;
    felim->resize(maillage->Get_n_nodes(),1);

    p_K = new Eigen::SparseMatrix<double> (maillage->Get_n_nodes(),maillage->Get_n_nodes());
    //p_Kelim = new Eigen::SparseMatrix<double> (maillage->n_nodes,maillage->n_nodes);

    diag = new Eigen::DiagonalMatrix<double, Eigen::Dynamic>;

    partition_noeud = new int[maillage->Get_n_nodes()];
    for (int i=0;i<maillage->Get_n_nodes();i++)
    {
        partition_noeud[i]=-1;
    }

    voisins_interface.resize(maillage->Get_nb_partitions());
    voisins_partition.resize(maillage->Get_nb_partitions());

    //cout<<"Task : "<<rang<< " voisins_interface contient "<<voisins_interface.size()<<" lignes"<<endl;

    //cout <<"Task : "<<rang<< " Initialisation de tous les parametres de la classe" << endl;
    //cout <<"Task : "<<rang<< " remplissage des vecteurs voisins_partition et voisins_interface"<< endl;


    calcul_voisins();

    //cout<<"affichage de voisins_interface"<<endl;
    //affiche_vector(voisins_interface);

    //cout<<"affichage de voisins_partition"<<endl;
    //affiche_vector(voisins_partition);


    /* Calcul de la solution exacte */

    for(int ind_node=0;ind_node<maillage->Get_n_nodes();ind_node++)
    {
        uexa->coeffRef(ind_node,0)+=calcul_uexa(maillage->Get_nodes_coords()[3*ind_node],maillage->Get_nodes_coords()[3*ind_node+1]);
    }

    //cout <<"Task : "<<rang<< " calcul de la solution exacte"<<endl;

    /* Initialisation des conditions au bord, dans un premier temps a une constante*/

    for(int ind_node=0;ind_node<maillage->Get_n_nodes();ind_node++)
    {
        if (maillage->Get_nodes_ref()[ind_node]!=0)
        {
            double x=maillage->Get_nodes_coords()[3*ind_node+0];
            double y=maillage->Get_nodes_coords()[3*ind_node+1];
            double addedCoeff = calcul_g(x,y);
            g->coeffRef(ind_node,0)+=addedCoeff;
        }
    }

    for(int ind_node=0;ind_node<maillage->Get_n_nodes();ind_node++)
    {
        if (partition_noeud[ind_node]==rang)
        {
            u->coeffRef(ind_node,0)+=partition_noeud[ind_node];
        }
        else
        {
            u->coeffRef(ind_node,0)+=-1;
        }
    }

    //cout<<"Task : "<<rang<< " calcul des conditions au bords"<<endl;

    /* Assemblage de la matrice de rigiditÃ© par parcours de tous les triangles */

    assemblage(rang);

    //cout<<"Task : "<<rang<< " assemblage de la matrice de rigidite reussie."<<endl;

    //affich(*p_K);

    /* On garde en mÃ©moire la matrice assemblÃ©e avant pseudo Ã©limination pour calculer l'erreur H1 plus tard */

    Eigen::SparseMatrix<double> K_err = *p_K;

    /* Le second membre total, prenant en compte les f et g de la formulation variationnelle
     * felim a ete obtenu par formules de quadrature. Le second membre en g est obtenu par interpolation
     */

    *felim=*felim-(*p_K) * (*g);

    /* Mise en oeuvre de la pseudo elimination */

    for(int i=0;i<maillage->Get_n_nodes();i++)
    {
        /* seuls les elements du bord sur felim doivent etre changes */

        if (maillage->Get_nodes_ref()[i]==1)
        {
            double nouvCoeff = p_K->coeffRef(i,i)*g->coeffRef(i,0);
            felim->coeffRef(i,0)=nouvCoeff;
            for(int j=0;j<maillage->Get_n_nodes();j++)
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

    //cout <<"Task : "<<rang<< " etape de pseudo elimination terminee."<<endl;
    //affich(*p_K);

    /* On sotcke la diagonale de la matrice de rigidité pour les itérations */

    (*diag).diagonal() = (*p_K).diagonal();

    /* On retire cette diagonale de la matrice de rigidité car le produit matriciel au cours
     * des itérations ne les fait pas intervenir ; une boucle pourrait etre évitée en utilisant
     * les fonctions diagonal de eigen mais la version installee empeche d'ajouter des matrices creuses
     * et des matrices denses
     * */

    for (int i=0;i<maillage->Get_n_nodes();i++)
    {
        p_K->coeffRef(i,i)=0;
    }

    //cout<<"voici la matrice de rigidite sans diagonale, telle qu'elle intervient dans les calculs"<<endl;
    //affich(*p_K);

}



void Probleme::calcul_voisins()
{
    for (int ind_triangle=0;ind_triangle<maillage->Get_n_triangles();ind_triangle++)
    {
        int ind_pt1=(maillage->Get_triangles_sommets())[3*ind_triangle]-1;
        int ind_pt2=maillage->Get_triangles_sommets()[3*ind_triangle+1]-1;
        int ind_pt3=maillage->Get_triangles_sommets()[3*ind_triangle+2]-1;

        int numero_partition_triangle=maillage->Get_partition_ref()[maillage->Get_n_elems()-maillage->Get_n_triangles()+ind_triangle];

        //cout<<"________________________________"<<endl;
        //cout<<"on regarde le triangle de noeuds "<<ind_pt1+1<<" "<<ind_pt2+1<<" "<<ind_pt3+1<<endl;
        //cout<<"c'est l element "<<maillage->Get_n_elems()-maillage->Get_n_triangles()+ind_triangle+1<<endl;;


        /* Un noeud Ã  l'interface est reconnu par son appartenance Ã  deux partitions distinctes. */

        //cout<<"le sommet "<<ind_pt1+1<<" est pour l'instant sur la partition "<<partition_noeud[ind_pt1]<<endl;

        if (partition_noeud[ind_pt1]==-1)
        {
            //cout<<"il est non initialise donc sa partition est maintenant celle de son triangle"<<endl;
            partition_noeud[ind_pt1]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt1]==numero_partition_triangle)
        {
            //cout<<"il est deja sur la partition de son triangle"<<endl;
        }
        else
        {
            //cout<<"il est initialise sur une partition differente du triangle parcouru"<<endl;
            partition_noeud[ind_pt1]=0;
        }

        //cout<<"le sommet "<<ind_pt2+1<<" est pour l'instant sur la partition "<<partition_noeud[ind_pt2]<<endl;

        if (partition_noeud[ind_pt2]==-1)
        {
            //cout<<"il est non initialise donc sa partition est maintenant celle de son triangle"<<endl;
            partition_noeud[ind_pt2]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt2]==numero_partition_triangle)
        {
            //cout<<"il est deja sur la partition de son triangle"<<endl;
        }
        else
        {
            //cout<<"il est initialise sur une partition differente du triangle parcouru"<<endl;
            partition_noeud[ind_pt2]=0;
        }

        //cout<<"le sommet "<<ind_pt3+1<<" est pour l'instant sur la partition "<<partition_noeud[ind_pt3]<<endl;

        if (partition_noeud[ind_pt3]==-1)
        {
            //cout<<"il est non initialise donc sa partition est maintenant celle de son triangle"<<endl;
            partition_noeud[ind_pt3]=numero_partition_triangle;
        }
        else if(partition_noeud[ind_pt3]==numero_partition_triangle)
        {
            //cout<<"il est deja sur la partition de son triangle"<<endl;
        }
        else
        {
            //cout<<"il est initialise sur une partition differente du triangle parcouru"<<endl;
            partition_noeud[ind_pt3]=0;
        }

        //cout<<"apres parcours de l element "<<maillage->Get_n_elems()-maillage->Get_n_triangles()+ind_triangle+1<<" les noeuds "<<ind_pt1+1<<", "<<ind_pt2+1<<" et "<<ind_pt3+1<<endl
        //<<" sont sur les partitions "<<partition_noeud[ind_pt1]<<", "<<partition_noeud[ind_pt2]<<" et "<<partition_noeud[ind_pt3]<<endl;
    }

    /* On sait maintenant quel noeud appartient Ã  quel triangle. En parcourant Ã  nouveau les triangles,
     * en s'arretant sur ceux qui sont au contact de l'interface sans y etre totalement inclus, on peut
     * donc dÃ©terminer les vecteurs voisins_interface et voisins_partition
     */

    //cout<<"________________________________"<<endl;
    //cout<<"PASSAGE AU REMPLISSAGE DES VOISINS"<<endl;
    //cout<<"voisins_interface contient "<<voisins_interface.size()<<" lignes"<<endl;


    for (int ind_triangle=0;ind_triangle<maillage->Get_n_triangles();ind_triangle++)
    {
        int ind_pt1=maillage->Get_triangles_sommets()[3*ind_triangle]-1; 
        int ind_pt2=maillage->Get_triangles_sommets()[3*ind_triangle+1]-1;
        int ind_pt3=maillage->Get_triangles_sommets()[3*ind_triangle+2]-1;

        //cout<<"on regarde le triangle de noeuds "<<ind_pt1+1<<" "<<ind_pt2+1<<" "<<ind_pt3+1<<endl;
        //cout<<"c'est l element "<<maillage->Get_n_elems()-maillage->Get_n_triangles()+ind_triangle+1<<endl;;

        int numero_partition_triangle=maillage->Get_partition_ref()[maillage->Get_n_elems()-maillage->Get_n_triangles()+ind_triangle];
        //cout<<"il est sur la partition "<<numero_partition_triangle<<endl;

        if ((partition_noeud[ind_pt1]==0 || partition_noeud[ind_pt2]==0 || partition_noeud[ind_pt3]==0) &&
                (!(partition_noeud[ind_pt1]==0 && partition_noeud[ind_pt2]==0 && partition_noeud[ind_pt3]==0)))
        {

            //cout<<"l'element "<<maillage->Get_n_elems()-maillage->Get_n_triangles()+ind_triangle+1<<" est au contact de l'interface"<<endl;
            if (partition_noeud[ind_pt2]!=0)
            {
                /* Si le noeud n'est pas sur l'interface, il en est voisin */
                voisins_interface[partition_noeud[ind_pt2]-1].push_back(ind_pt2+1);
                //cout<<"le noeud "<<ind_pt2+1<<" n'est pas sur l'interface et appartient a la partition "<<partition_noeud[ind_pt2]<<endl;
                //cout<<"vecteur voisins_interface"<<endl;
                //affiche_vector(voisins_interface);
            }
            else
            {
                /* Si le noeud est sur l'interface, il est voisin de la partition du triangle auquel il appartient */
                voisins_partition[numero_partition_triangle-1].push_back(ind_pt2+1);
                //cout<<"le noeud "<<ind_pt2+1<<" est sur l'interface et appartient au triangle dont la partition est "<<numero_partition_triangle<<endl;
                //cout<<"vecteur voisins_partition"<<endl;
                //affiche_vector(voisins_partition);
            }
            if (partition_noeud[ind_pt3]!=0)
            {
                voisins_interface[partition_noeud[ind_pt3]-1].push_back(ind_pt3+1);
                //cout<<"le noeud "<<ind_pt3+1<<" n'est pas sur l'interface et appartient a la partition "<<partition_noeud[ind_pt3]<<endl;
                //cout<<"vecteur voisins_interface"<<endl;
                //affiche_vector(voisins_interface);
            }
            else
            {
                voisins_partition[numero_partition_triangle-1].push_back(ind_pt3+1);
                //cout<<"le noeud "<<ind_pt3+1<<" est sur l'interface et appartient au triangle dont la partition est "<<numero_partition_triangle<<endl;
                //cout<<"vecteur voisins_partition"<<endl;
                //affiche_vector(voisins_partition);
            }
            if (partition_noeud[ind_pt1]!=0)
            {
                voisins_interface[partition_noeud[ind_pt1]-1].push_back(ind_pt1+1);
                //cout<<"le noeud "<<ind_pt1+1<<" n'est pas sur l'interface et appartient a la partition "<<partition_noeud[ind_pt1]<<endl;
                //cout<<"vecteur voisins_interface"<<endl;
                //affiche_vector(voisins_interface);
            }
            else
            {
                voisins_partition[numero_partition_triangle-1].push_back(ind_pt1+1);
                //cout<<"le noeud "<<ind_pt1+1<<" est sur l'interface et appartient au triangle dont la partition est "<<numero_partition_triangle<<endl;
                //cout<<"vecteur voisins_partition"<<endl;
                //affiche_vector(voisins_partition);
            }
        }
        else
        {
            //cout<<"Il n'est pas au contact de l'interface"<<endl;
        }

     }

    /* On trie chaque liste de noeuds des vecteurs voisins_partition et voisins_interface par ordre
     * croissant afin de pouvoir les utiliser dans les communications
     * */

    vector<vector<int> >::iterator vect_it; 

    for (vect_it=voisins_partition.begin();vect_it!=voisins_partition.end();vect_it++)
    {
        std::sort((*vect_it).begin(),(*vect_it).end());
        
        vector<int>::iterator it;

        it = std::unique((*vect_it).begin(),(*vect_it).end());
        (*vect_it).resize(std::distance((*vect_it).begin(),it));

    }

    for (vect_it=voisins_interface.begin();vect_it!=voisins_interface.end();vect_it++)
    {
        std::sort((*vect_it).begin(),(*vect_it).end());

        vector<int>::iterator it;

        it = std::unique((*vect_it).begin(),(*vect_it).end());
        (*vect_it).resize(std::distance((*vect_it).begin(),it));
    }



    //cout<<"voici le vecteur voisins_partition, contenant sur la i-eme ligne la liste ordonnée des noeuds de l'interface voisins de la partition i"<<endl;
    //affiche_vector(voisins_partition);

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

    for(int ind_triangle=0;ind_triangle<maillage->Get_n_triangles();ind_triangle++)
    {
        int ind_pt1=maillage->Get_triangles_sommets()[3*ind_triangle];
        int ind_pt2=maillage->Get_triangles_sommets()[3*ind_triangle+1];
        int ind_pt3=maillage->Get_triangles_sommets()[3*ind_triangle+2];

        /* Les -3 sont dus au fait que les valeurs commencent Ã  l'indice 0 or min(ind_pt)=1 */

        double x1=maillage->Get_nodes_coords()[0+ind_pt1*3-3];
        double y1=maillage->Get_nodes_coords()[1+ind_pt1*3-3];
        double z1=maillage->Get_nodes_coords()[2+ind_pt1*3-3];
        double x2=maillage->Get_nodes_coords()[0+ind_pt2*3-3];
        double y2=maillage->Get_nodes_coords()[1+ind_pt2*3-3];
        double z2=maillage->Get_nodes_coords()[2+ind_pt2*3-3];
        double x3=maillage->Get_nodes_coords()[0+ind_pt3*3-3];
        double y3=maillage->Get_nodes_coords()[1+ind_pt3*3-3];
        double z3=maillage->Get_nodes_coords()[2+ind_pt3*3-3];

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
        int ind_global=maillage->Get_triangles_sommets()[3*ind_triangle+j];
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
        int ind_global_1=maillage->Get_triangles_sommets()[3*ind_triangle+i];
        for(int j=0;j<3;j++)
        {
            int ind_global_2=maillage->Get_triangles_sommets()[3*ind_triangle+j];
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
    //delete p_Kelim;
    //p_Kelim=0;
    delete u;
    u=0;
}


////// Méthdoes GET et SET pour l'encapsulation
    //** Méthodes SET

    void Probleme::Set_maillage (Maillage* _maillage) {
        maillage=_maillage;
    }

    void Probleme::Set_uexa (VectorXd* _uexa) {
        uexa=_uexa;
    }

    void Probleme::Set_g (VectorXd* _g) {
        g=_g;
    }

    void Probleme::Set_u (VectorXd* _u) {
        u=_u;
    }

    void Probleme::Set_felim (VectorXd* _felim) {
        felim=_felim;
    }

    void Probleme::Set_partition_noeud (int* _partition_noeud) {
        partition_noeud = _partition_noeud;
    }

    void Probleme::Set_voisins_interface (vector<vector<int> > _voisins_interface)  {
        voisins_interface=_voisins_interface;
    }

    void Probleme::Set_voisins_partition (vector<vector<int> > _voisins_partition) {
        voisins_partition=_voisins_partition;
    }

    void Probleme::Set_p_K (Eigen::SparseMatrix<double>* _p_K) {
        p_K=_p_K;
    }

    void Probleme::Set_p_Kelim (Eigen::SparseMatrix<double>* _p_Kelim) {
        p_Kelim=_p_Kelim;
    }

    void Probleme::Set_diag (Eigen::DiagonalMatrix<double,Eigen::Dynamic>* _diag) {
        diag=_diag;
    }

    //** Méthodes GET
    
    Maillage* Probleme::Get_maillage (){
        return maillage;
    }

    VectorXd* Probleme::Get_uexa (){
        return uexa;
    }

    VectorXd* Probleme::Get_g (){
        return g;
    }

    VectorXd* Probleme::Get_u (){
        return u;
    }

    VectorXd* Probleme::Get_felim () {
        return felim;
    }

    int* Probleme::Get_partition_noeud (){
        return partition_noeud;
    }

    vector<vector<int> > Probleme::Get_voisins_partition (){
        return voisins_partition;
    }

    vector<vector<int> > Probleme::Get_voisins_interface (){
        return voisins_interface;
    }

    Eigen::SparseMatrix<double>* Probleme::Get_p_K (){
        return p_K;
    }

    Eigen::SparseMatrix<double>* Probleme::Get_p_Kelim () {
        return p_Kelim;
    }

    Eigen::DiagonalMatrix<double,Eigen::Dynamic>* Probleme::Get_diag () {
        return diag;
    }
