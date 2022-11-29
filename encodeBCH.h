#ifndef _ENCODEBCH_H_
#define _ENCODEBCH_H_

#include <math.h>
#include <bitset>
#include <vector>
#include <iterator>

static unsigned int errcormax=0; //Nombre d'erreur maximum à introduire
  
class Encodebch
{
  unsigned int n;//longueur du mot encodé (g(x)*m(x))
  unsigned int size;//Taille du message m(x) à transmettre

  std:: vector<unsigned int> mesencode;//vecteur de bits contenant le message encodé
  std:: vector<unsigned int> g; //Tableau contenant le polynome générateur g
  
 public:
  Encodebch(std:: vector<unsigned int> message, unsigned int degre)
  {
    //On définit la longueur du message
    n=pow(2,degre)-1;
    //On resize la longueur du polynome generateur
    g.resize(n-1);
    //Selon le degré on affiche les bonnes informations
    if (degre==5)
    {
      errcormax=7;
      printf("Le nombre maximum d'erreurs corrrigeable pour un message de longueur %d est %d\n", n, errcormax);
      printf("Le polynome générateur associé est ");
      g=octal_binary(313365047);
      printf("\n");
      printf("La distance de hamming minimum entre tous les mots du code (ensemble des mots divisible par g(x)) sera alors 15\n");
    }
    if (degre==4)
    {
      errcormax=3;
      printf("Le nombre maximum d'erreurs corrrigeable pour un message de longueur %d est %d\n", n, errcormax);
      printf("Le polynome générateur associé est ");
      g=octal_binary(2467);
      printf("\n");
      printf("La distance de hamming minimum entre tous les mots du code (ensemble des mots divisible par g(x)) sera alors 8\n");
    }
    if(degre==6)
      {
	printf("nul\n\n\n");
	 g=octal_binary(313365047);
      }
  }
  
  /* Fonction pour convertir un nombre  octal en binaire.*/
  std:: vector<unsigned int> octal_binary(unsigned int octal)
  {
    unsigned int decimal=0, i=0, j=0;
    std:: vector<unsigned int> h(n-2,0);
    h[0]=0;
    while (octal!=0)
    {
        decimal+=(octal%10)*pow(8,i);
        ++i;
        octal/=10;
    }
    /* At this point, the decimal variable contains corresponding decimal value of that octal number. */
    i=1;
    while(decimal!=0)
    {
        h[j]+=(decimal%2);
        decimal/=2;
	j++;
    }
    for (unsigned int k=0; k<n-2; k++)
    {
      printf("%d", h[n-3-k]);
    }
    printf(" (en binaire)\n");
    return h;
  }

  std:: vector<unsigned int> encodeBCH (std:: vector<unsigned int> a)
  {
    //vecteur contenant notre polynome générateur translaté d'un certain de 0 vers la gauche
    std:: vector<unsigned int> resultint;
    //Resize de la taille du message encodé suivant la logique de multiplication de 2 vecteurs a et g
    //(taille max du vecteur = taille de a + taille de g
    mesencode.resize(a.size()+g.size());
    unsigned int k=mesencode.size();
    //On initialise nos valeurs dans le vecteur qui contiendra le  message encodé à 0
    for (unsigned int i=0; i<mesencode.size(); i++)
    {
	mesencode[i]=0;
    }
    //On parcourt les élément de notre vecteur a (ici le message à encoder)
    for (unsigned int i=0; i<a.size(); i++)
      {
	//Dès qu'un bit est à 1 dans a
	if(a[i]==1)
	  {
	    //On resize le résultat intermédiaire avec décalage de i correspondant à la position du 1 trouvé dans a
	    //On prépare donc le vecteur resultint a accueillir le polynome générateur translaté de i
	    resultint.resize(g.size()+a.size()-i);
	    //On remplit les premiers bits de notre vecteur avec des 0 suivant les regles de multiplications
	    //de mots binaires, on met autant de 0 que la position du 1 trouvé
	    for (unsigned int j=0; j<a.size()-1-i; j++)
		  {
		    resultint[j]=0;
		  }
	    //On recopie à la suite des 0 mis précédemment, le polynome générateur. 
	    for (unsigned int j=a.size()-1-i; j<resultint.size(); j++)
		  {
		    resultint[j]=g[j-a.size()+1+i];
		  }
	    for (unsigned int j=0; j<resultint.size(); j++)
	      {
		//On ajoute dans notre vecteur contenant le message encodé, les resultats intermédiaire
		//qui ne sont autres que tous les vecteurs générateur décalé de la position i
		//(position ou se trouve les 1 dans le vecteur message (a))
		mesencode[j]+=resultint[j];
		//On s'assure de respecter mod(x^n - 1)
		if((resultint.size()-1-j>=n)&&(resultint[resultint.size()-1-j]==1))
		  {
		    mesencode[resultint.size()-1-j]=0;
		    mesencode[resultint.size()-1-j-n]+=1;	    
		  }
		//On s'assure de respecter le mod(2)
		if(mesencode[j]==2)
		  {
		    mesencode[j]=0;
		  }
	      }
	  }//if(a[i]...)
      }//for
    //On resize notre vecteur de bits encodé sans tenir compte des 0 résiduel en fin de mot.
    //En même temps on s'assure que le cas où le mot binaire à encoder est 0...00
    //soit aussi traiter par notre programme (seconde condition du while)
    while((mesencode[k-1]==0)&&(k-1!=0))
    {   
      k--;
    }
    mesencode.resize(k);
    return mesencode;
  }//méthode encodeBCH
  void afficher(std:: vector<unsigned int> c)
  {
    printf("Le message encodé est ");
    for (unsigned int j=0; j<c.size(); j++)
    {
      printf("%d", c[c.size()-1-j]);
    }
    printf("\n\n");
  }//méthode afficher

};//Fin de classe 

#endif
