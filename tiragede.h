#ifndef _TIRAGEDE_H_
#define _TIRAGEDE_H_

#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <vector>

using namespace std;

#define MAX 1//Max du tirage aléatoire
#define MIN 0//Min du tirage aléatoire

class Tiragede
{
  unsigned int size;//taille du message à envoyer
  vector<unsigned int> I;//vecteur de bits à envoyer
 public:
  //Constructeur de ma classe tirage
  Tiragede(unsigned int _size): size(_size)
  {
    //On change la taille du vecteur en fonction de la taille souhaité pour l'envoi de données
    I.resize(size);
  }
  std:: vector<unsigned int> gettirage(unsigned int _size)
  {
    //Générateur de nombre aléatoire binaire
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(MIN, MAX);
    for(unsigned int i=0; i<I.size(); i++)//on affecte nos valeurs de '0' ou '1' dans le vecteur
    {
      I[i]=distr(eng);
    }
    return I;
  }
  void afficher(std:: vector<unsigned int> I)
  {
    printf("Le message à encoder est ");
      for (unsigned int i=0; i<size; i++)
	{
	  printf("%d", I[i]);
	}
      printf("\n");
  }

  
};
#endif

