#ifndef _CALCALPHA_H_
#define _CALCALPHA_H_

#include <math.h>
#include <bitset>
#include <vector>

class Alpha
{
  unsigned int deg;//degrée du polynome (longueur du message)
  unsigned int nblig; //nb de ligne de notre tableau 2D
  unsigned int nbpair; //nombre de 1 pair pour éliminer les polynomes de la liste des potentielles primitives
  unsigned int cptlig; //compteur pour savoir le nombre de ligne du tableau allpolyredu (pour la libération de ressources)
  vector<unsigned int> elem;//vecteur qui contiendra tous les élements généré à partir de l'élément primitif
  vector<unsigned int> elem0;//vecteur qui contiendra tous les élements généré à partir de l'élément primitif
  vector<unsigned int> result;
  vector<unsigned int> polyminimal;//Vecteur contenant le polynome minimal
  unsigned int** poly; //Tableau 2D contenant tous les polynomes possible de degrées "deg"
  unsigned int** irredu; //Tableau 2D intermediaire contenant les polynomes possiblement irreductible et donc possiblement primitifs
  std:: vector<unsigned int> alpha; //Tableau contenant les valeurs construite à partir de l'élément primitifs
  
  unsigned int  binary, remainder, product ;
 
 public:

  //Cette méthode créée un tableau qui contiendra tous les polynome possible du degré souhaité en binaire
  //(Ex : x²+1 de degré 2 => (101) en binaire) 
  void toutpoly(unsigned int nbdeg)
  {
    deg=nbdeg+1;
     //Si on a polynome de degrée 5 alors il y a 2⁵-1 polynome possible (avec des coefficients 1).
    //On l'initialise à 2⁵ car les boucles for employé ensuite vont de 0 à 2⁵-1
    nblig=pow(2,deg);
    //printf("nblig=%d\n\n", nblig);
    //Initialisation du tableau qui va contenir tous les polynomes avec leurs coefficient sur chaque case  
    poly = new unsigned int*[nblig];
    for(unsigned int i=0; i<nblig; i++)
    {
     poly[i] = new unsigned int[deg];
    }
    // Remplissage du tableau 2D avec tous les polynome possible avec un coefficient unitaire ou égale à 0.
    for(unsigned int i=0;i<nblig;i++)
    {
      //Génération de tous les nombres entre 0 et 2^deg-1, en binaire
      //Qui seront alors tous les polynomes possibles
      unsigned int j=0;
      unsigned int decimal =i;
      while (decimal != 0 && j<deg) {
	  binary = 0;
	  product = 1;
	  remainder = decimal % 2;
	  binary = binary + (remainder * product);
	  poly[i][deg-j-1]=binary;
	  j++;
	  decimal = decimal / 2;
	  product *= 10;
      }//while		
    }//for 
  }//méthode toutpoly

  //Cette méthode va "élaguer" notre tableau de polynome obtenu avec toutpoly
  //Suivant 3 critères simples définis plus bas.
  void polyprim()
  {
    //Voir les attributs de la classe
    nbpair=0;
    cptlig=0;
    //Définition de k pour parcourir + librement les éléments du tableau irredu plus bas
    unsigned int k=0;
    
    //Initialisation du tableau qui va contenir tous les polynomes avec leurs coefficient sur chaque case  
    irredu = new unsigned int*[nblig];
    for(unsigned int i=0; i<nblig; i++)
    {
      irredu[i] = new unsigned int[deg];
    }
    //Initialisation d'alpha qui va contenir tous les nombres générables par le polynome possiblement primitifs  
    alpha.resize(nblig);
    //On parcourt les lignes de notre tableau allpoly qui représente chacun un polynome possible
    for(unsigned int i =0; i<nblig; i++)
    {
      //On verifie si le dernier coefficient du polynome est bien égale à 1
      //sinon ca veut dire qu'on a un polynome de la forme x^n+x^(n-1)+...+x
      //qui est factorisable par x
      //En même temps on s'assure que le coeffcient de plus haut degré
      //ne soit pas à 0 car cela signifierait que le polynome serait alors de
	//degré deg-1 qui serait donc définis dans un espace a deg-1 dimensions
	//il serait alors impossible de fabriquer des éléments de dimension supérieures
	//(Si un point est défini sur x,y alors on ne peut construire une composante en 3D)
      if((poly[i][0]!=0)&&(poly[i][deg-1]==1))
      {
	//On parcourt les colonnes de notre tableau qui sont les coefficients de notre polynome
	for(unsigned int j=0; j<deg; j++)
	{
	  //On va compter le nombre de coeffcient à 1
	  if(poly[i][j]==1)
	  {
	    nbpair++;
	  }
	}
	//Si le nombre de coeffcient unitaire est paire
	//alors cela signifie que le polynome est reductible
	//En effet, notre corps contient 2 éléments {0,1} qui nous permettrons de fabriquer l'ensemble des valeur de 0 à 1 ((taille du message)^n)-1
	//Donc notre ensemble de travail n'est composé que de 0 et 1
	//Ainsi, si un polynome admet une racine dans son ensemble de définitions (0 ou 1)
	//alors il est reductible et donc impossible qu'il soit le polynome minimal
	//Par exemple, x⁴+x³+x²+1 => (11101) en binaire et il a bien un nombre de 1 paire
	//Or comme nous travaillons dans un corps de galois (= corps finis)  à deux élements
	//notre inconnue a soit 0 ou 1 comme valeur possible
	//Si on remplace, p(1)=1⁴+1³+1²+1=4
	//Cependant il faut noter que dans un corps de galois à 2 élements nous travaillons en modulo 2
	//donc 4(mod 2)=0 => notre polynome est bien reductible dans le corps de galois à 2 éléments car il existe une racine dans le corps de définition
	//Le polynome est donc factorisable/reductible=> il ne peut être primitif.

	if ((nbpair%2==0)||(nbpair==1))
	{
	  nbpair=0;
	}
	else
	{
	  //Sinon si toutes les précédentes conditions ne sont pas validé alors le polynomes peut-être primitifs
	  //On stocke alors le polynome dans un nouveau tableau avant de passer à la ligne/polynome suivant.
	  for(unsigned int j=0; j<deg; j++)
	  {
	    irredu[k][j]=poly[i][j];
	  }
	  k++;
	  cptlig++;
	  nbpair=0;
	}//else
      }//if
    }//for(i)
  }//méthode polyprim

  //Cette nous permet d'obtenir le polynome minimal
  //(polynome primitif avec les plus faibles degrés pour ses coeff)
  std:: vector<unsigned int> polymin()
  { 
    unsigned int dec=0;
   
    //On initialise la taille de nos vecteurs
    elem.resize(deg-1);
    elem0.resize(deg-1);
    polyminimal.resize(deg);
    result.resize(deg);
    //Tout polynome primitif admet un élément primitifs 
    for (unsigned int i= 0; i<cptlig; i++)//On parcourt nos potentiels primitifs (Ex : alpha⁴+alpha+1=0 avec alpha l'élément primitif et (10011) en binaire)
    {
      unsigned int ind=0;
       unsigned int cptid=0;
      //On stocke le premier polynome dans elem privé du bit lié au coeff dominant
       //En gros, on isole alpha⁴=-alpha-1=alpha+1=(0011) (1+1=0=>1=-1)
      for (unsigned int j=1; j<deg; j++)
      {
	elem0[j-1]=irredu[i][j];
	elem[j-1]=irredu[i][j];
      }
      //On parcourt les 2^deg de notre corps finis
      for(unsigned int l=0;l<pow(2,deg-1)-1;l++)
      {
	//On multiplie alpha⁴(=0011) par alpha(0010) pour avoir alpha⁵(ici result)
	//Qui nous donne 00110
	//Ce qui revient à opère un décalage sur alpha⁴ (ici elem)
	result[deg-1]=0;
	for (unsigned int k=0; k<deg-1;k++)
	{
	  result[k]=elem[k];
	}
	//On convertit alpha⁴ (elem en décimal)
	dec=binarytodec(elem);
       
	
	//On vérifie que alpha⁵ ne soit pas egale à une somme de alpha avec un degré 4
	//car il faudrait veiller à remplacer alpha⁴ par alpha+1
	if (result[0]==1)
	{
	  //Si on a bien un terme alpha⁴
	  //On parcourt notre vecteur de bit
	  for (unsigned int k=1; k<deg;k++)
	  {
	    //Et on remplace alpha⁴ par (alpha+1)
	    //Ce qui revient à additionné alpha+1 (=0010)
	    //à notre vecteur result qui contient le décalage opéré précédemment (00110)
	    //
	    elem[k-1]=result[k]+elem0[k-1];
	    //On tient compte du modulo 2 pour l'addition
	    if(elem[k-1]==2)
	    {
	      elem[k-1]=0;
	    }
	  }
	}//if
	//Sinon c'est que alpha⁴ n'est pas présent dans l'expression d'alpha⁵ (par exemple)
	else
	{
	  for (unsigned int k=1; k<deg;k++)
	  {
	    elem[k-1]=result[k];
	  }
	}//else
	//On convertit le vecteur de alpha^n en décimal
	alpha[ind]=dec;
	//dec=binarytodec(elem);
	//printf("alpha[%d] =%d\n",ind, alpha[ind]);
	//On le stocke dans le tableau des valeurs de alpha

	ind++;
      }//for de tous nos éléments
      //On va parcourir toutes les valeurs possible de notre CG(2^n) 
      for(unsigned int p=0; p<(nblig/2)-1; p++)
      {
	unsigned int deuxfois=0;
	//On va parcourir notre tableau avec les valeurs decimales lié aux différentes expression de alpha^n 
	for(unsigned int n=0; n<(nblig/2)-1; n++)
	{
	  //Si on trouve une équivalence
	  if(alpha[n]==p)
	  {
	    //On incrémente ce compteur qui s'assure qu'il n'y ait pas de doublon
	    deuxfois++;
	    //On incrémente ce compteur qui vérifie le nombre d'équivalence trouvé
	    cptid++;
	  }
	  //Si on a un doublon c'est que alpha^n a généré deux valeur décimales (et donc deux polynomes) identiques
	 
	  if (deuxfois==2)
	  {
	    goto exit;
	  }
	 }//for n
      }//for p
      //Si le nombre d'équivalence comptait est bien égale au cardinal de notre corps de galois
      //Alors le polynome est bien primitifs et comme les polynomes étaient bien classées par ordre croissant de degré
      //dans toutpoly => implique que le premier polynome primitif trouvé est bien le polynome minimal (avec les degrés les moins haut) 
      if (cptid==(nblig/2)-2)
      {
	for (unsigned int j=0; j<deg; j++)
	{
	  polyminimal[j]=irredu[i][j];
	}
	break;
      }
      //Le polynome ne peut-être primitifs
      //On casse la boucle et on passe au polynome suivant
     exit:
      printf ("Le polynome ");
      for (unsigned int j=0; j<deg; j++)
      {
	printf("%d",irredu[i][j]);
      }
      printf (" n'est pas primitifs\n\n");
    }//troisieme for des polynomes primitifs
    return alpha;
  }//Méthode polymin
  
  //Méthode pour convertir du binaire en décimal
  unsigned int  binarytodec (std::vector<unsigned int> vec)
  {
    unsigned int dec=0;
    for(unsigned int i=1; i<deg+1; i++)
      {
	if(vec[deg-i]==1)
	  {
	    dec+=pow(2,i-1);
	  }
      }
    dec=dec/2;
    return dec;
  }
  
 void afficher()
 {
   //Verif pour être sur de ne pas avoir un tableau de taille nulle
   if(deg == 0)
   {
     printf("Erreur : taille nulle");
     return;
   }
   else
   {
     /* printf("Tous les polynomes possible de degrès %d sont \n\n\n", deg-1);
     for (unsigned int i=0; i<nblig; i++)
     {
       for (unsigned int j=0; j<deg; j++)
       {
	 printf("%d", poly[i][j]);
       }
       printf("\n");
       }*/
     /*printf("Les polynomes qui sont peut-être primitifs sont \n\n\n");
     for (unsigned int i=0; i<cptlig; i++)
     {
       for (unsigned int j=0; j<deg; j++)
       {
	 printf("%d", irredu[i][j]);
       }
       printf("\n");
       }*/
     /*unsigned l=0;
     printf("Les valeurs de alpha sont :\n");
     for(unsigned int n=0; n<(nblig/2)-1; n++)
     {
       printf("%d pour alpha degré %d\n",alpha[n], deg+l-1);
       l++;
     }
     printf("\n\n");*/
     printf ("Le polynome qui s'annule pour l'élement primitif alpha  de degre %d est\n", deg-1);
     for (unsigned int j=0; j<deg; j++)
     {
       printf("%d",polyminimal[j]);
       }
     printf ("\nIl est à noter que ce polynome s'annule aussi pour les élément alpha², alpha⁴, alpha⁸, alpha¹⁶, alpha³²...jusqu'à une puissance inf à %d \n", deg-1);
   }//else
 }//Méthode

 ~Alpha()
 {
   //On libère la ressource allouée 
   if(deg != 0)
   {
     for(unsigned int i=0; i<deg; i++)
     {
        delete[] poly[i];
	delete[] irredu[i];
     }
     delete[] poly;
     delete[] irredu;
   }
 }
};//Fin de classe 

#endif
