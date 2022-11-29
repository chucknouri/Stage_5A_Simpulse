#ifndef _AJOUTERR_H_
#define _AJOUTERR_H_

#include <math.h>
#include <bitset>
#include <vector>
#include <iterator>
#include <random>

static unsigned int nberreel=0;
class AjoutErr
{
  std:: vector<unsigned int> messaerr;//vecteur contenant notre message encodé avec les erreurs
  std:: vector<unsigned int> poserr;//Vecteur contanant la position de(s) erreur(s) ajoutée(s)
  //unsigned int nberreel;//nombre d'erreur dans le message
 public:
  std:: vector<unsigned int> ajouterr(std:: vector<unsigned int> messaencode, unsigned int errcormax)
  {
    printf("errcormax=%d\n\n\n", errcormax);
    messaerr.resize(messaencode.size());
    //nberreel=0;
    unsigned int k=0;
    //On génére aléatoirement un nombre entier entre 0 et 4
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(0,4);
    //On copie notre message à encoder dans le message erreur
    for (unsigned int j=0; j<messaerr.size(); j++)
    {
	messaerr[j]=messaencode[j];
    }
    //On commence à ajouter 7 erreurs (le nb max d'erreur coorigeable)
    for(unsigned int i=0; i<messaerr.size(); i++)
    {
      //Si le nombre tiré au hasard est 1 ou si le nombre max d'erreur n'est pas encore généré
      if ((distr(eng)==1)&&(nberreel<errcormax+1))
      {
	nberreel++;//On incrémente le nombre d'erreur 
	messaerr[i]=messaencode[i]+1;//On change le bit de valeur en ajoutant 1
	//On resize pour ajouter une case de plus et y stocker la valeur de la position d'erreur
	poserr.resize(k+1);
	//On stocke i qu'est la position de l'erreur
	poserr[k]=i;
	//On incrémente k pour préparer la case suivante
	k++;
	//On tient compte du mod(2)
	if (messaerr[i]==2)
	{
	  messaerr[i]=0;
	}
      }
    }
    return messaerr;
  }
  void afficher(std:: vector<unsigned int> messagerreur)
  {
    printf("Le message erronée est ");
    
    for (unsigned int j=0; j<messagerreur.size(); j++)
    {
      printf("%d", messagerreur[messagerreur.size()-1-j]);
    }
    printf("\n\n");
    printf("Il y a %d erreurs dans notre message en position(s) ", nberreel);
    for(unsigned int j=0; j<poserr.size(); j++)
    {
      printf("%d ", poserr[j]);
    }
    printf("\n\n");
  } 
};//Fin de classe 

#endif
