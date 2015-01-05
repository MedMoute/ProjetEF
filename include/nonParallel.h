/*!
 * \file nonParallel.h
 * \brief Header des fonctions n'exploitant pas la paraléllisation dans le solveur
 *
 * Contient principalement des fonctions d'affichage et d'écriture
 *
 * \version 1.0
 */

#ifndef NONPARALLEL_H
#define NONPARALLEL_H

#include <Eigen/Sparse>
#include "probleme.h"


//void calcul(VectorXd, VectorXd, Eigen::SparseMatrix<double>, Eigen::DiagonalMatrix<double, Eigen::Dynamic>, VectorXd);

 /*!
 * \fn void affich(Eigen::SparseMatrix<double> _mat)
 * \brief Fonction d'affichage de matrice
 * 
 * Fonction d'affichage d'une matrice creuse (SparseMatrix) définie dans la doc d'Eigen 
 *
 * \param _mat Matrice à afficher dans la sortie standard (à ne pas utiliser sur de gros maillages)
 */
void affich(Eigen::SparseMatrix<double>);

 /*!
 * \fn void affichVector(VectorXd V)
 * \brief Fonction d'affichage de vecteur
 *
 * Fonction d'affichage d'un vecteur présent sous forme de VectorXd
 *
 * \param V vecteur à afficher dans la sortie standard (à ne pas utiliser sur de gros maillages)
 */

void affichVector(VectorXd);

 /*!
 * \fn void affiche_vector(vector<vector<int> > v)
 * \brief Fonction d'affichage de vecteur/matrice
 *
 * Fonction d'affichage d'un vecteur ou d'une matrice présent sous forme de vector<vector<int>>
 *
 * \param v vecteur à afficher dans la sortie standard (à ne pas utiliser sur de gros maillages)
 * Vecteur est un vector<vector<>>, et peut être vu comme un tableau de vecteurs
 */
void affiche_vector(vector<vector<int> >);

/*!
* \fn void output_vector(VectorXd u, string s)
* \brief fonction d'écriture d'un vecteur dans un fichier externe
*
* Fonction écrivant le contenu du VectorXd passé en premier argument dans un fichier nommé par la chaîne de caractère 
*passée en deuxième argument
*
*\param u Vecteur à écrire dans le fichier de sortie
*\param s Nom du fichier de sortie 
*/
void output_vector(VectorXd, string);

#endif // NONPARALLEL_H
