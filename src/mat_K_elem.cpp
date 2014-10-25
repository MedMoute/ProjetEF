#include "../include/mat_K_elem.h"

void mat_K_elem(double *mat_elem,double *coorneu, int *tab_local_global, int ind_triangle){
  

  double x1=coorneu[0+tab_local_global[3*ind_triangle+0]];
  double y1=coorneu[1+tab_local_global[3*ind_triangle+0]];
  double z1=coorneu[2+tab_local_global[3*ind_triangle+0]];
  double x2=coorneu[0+tab_local_global[3*ind_triangle+1]];
  double y2=coorneu[1+tab_local_global[3*ind_triangle+1]];
  double z2=coorneu[2+tab_local_global[3*ind_triangle+1]];
  double x3=coorneu[0+tab_local_global[3*ind_triangle+2]];
  double y3=coorneu[1+tab_local_global[3*ind_triangle+2]];
  double z3=coorneu[2+tab_local_global[3*ind_triangle+2]];

  double x12 = x2 - x1;
  double y12 = y2 - y1;
  double z12 = z2 - z1;
  double x13 = x3 - x1;
  double y13 = y3 - y1;
  double z13 = z3 - z1;
  double x23 = x3 - x2;
  double y23 = y3 - y2;
  //double z23 = z3 - z2;

  double tab_interm[] = {-y23,x23,y13,-x13,-y12,x12};

  double aire_triangle = (1/2)*sqrt(pow((y12*z13 - z12*y13),2)+pow((z12*x13-x12*z13),2)+pow((x12*y13-y12*x13),2));
  //rajouter une erreur si l'aire du triangle est nulle
  //a partir de maintenant on considere que la troisieme coord est toujours nulle
  
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat_elem[3*i+j]=(1/(4*aire_triangle))*(tab_interm[2*i]*tab_interm[2*j]+tab_interm[2*i+1]*tab_interm[2*j+1]);
    }
  }
}
    

  
  
