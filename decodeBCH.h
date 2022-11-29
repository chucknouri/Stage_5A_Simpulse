#ifndef _DECODEBCH_H_
#define _DECODEBCH_H_

#include <math.h>
#include <bitset>
#include <vector>
#include <iterator>
#include <random>

static unsigned int u=0;



class decodeBCH
{
  std:: vector<unsigned int> puiss;//Vecteur qui va contenir les différentes puissance (en binaire) des alphas présents dans notre expression erreur(alpha^i)
  std:: vector<unsigned int> alphaord;//vecteur qui va contenir les différentes puissances de alpha de notre CG (Corps de Galois)
                                      //et leur équivalence avec des puissances simple (qui vont de alpha⁰ à alpha^(degré-1)) et ainsi pouvoir les ajouter.
                                      //Actuellement notre tableau alpha en arguments de la fonction syndrome commence à alpha^degre mais en ordonnant il commencera
                                      // à alpha⁰ ce qui simplifiera la recherche d'équivalence par la suite
  std:: vector<unsigned int> syndrun;//sert pour les calculs intermédiaire du syndrome
  std::vector<unsigned int>syndr; //vecteur contenant les syndromes finaux
  unsigned int** dico;//Tableau en 2D avec sur la première colonne la puissance de alpha et sur la deuxieme, la valeur décimale correspondante au mot binaire de l'expression équivalente avec des puissances simples (mentionné dans le commentaire de alphaord). Ex : dico[0][0]=[0 pour alpha puissance 0] et dico[0][1]=[1 car alpha⁰=1]
                                                                                 //dico[1][0]=[1 pour alpha puissance 1] et dico[1][1]=[2 car alpha en binaire = 0010]
  unsigned int nblig;//entier contenant le degré max pour le mod(x^nblig-1)
  std::vector<std::vector<unsigned int> > Matsyndrinit;//Notre matrice des syndromes
  
  std::vector<std::vector<unsigned int> > Matsyndr;
  unsigned int nberr;//Le nombre d'erreur détecter

  std::vector<unsigned int> sigma;
  
  unsigned int allerr=0;
  
 public:
  void syndromes(std:: vector<unsigned int> messaerr, std:: vector <unsigned int> alpha, unsigned int deg)
  {
    unsigned int l=0;
    unsigned int j=0;
    unsigned int m=0;
    unsigned int p=0;
    unsigned int q=0;
    
    nblig=pow(2,deg)-1;
    alpha.resize(nblig);
    alphaord.resize(alpha.size());
    //On devrait avoir au moins 2*nberreurs+1 syndromes pour positionner le(s) erreur(s)
    //Mais dans le cas pratique on ignore le nombre d'erreurs qu'il peut y avoir
    syndr.resize(2*errcormax);
    puiss.resize(alpha.size());
    syndrun.resize(2*errcormax);
    
    //Initialisation du "dictionnaire" des alpha pour le CG(2^degre)
    dico=new unsigned int*[alphaord.size()];
    for (unsigned int i = 0;i<alphaord.size(); i++)
      {
	dico[i]=new unsigned int[2];
      }
    //On va parcourir le vecteur alpha obtenu dans calcalpha.
    for(unsigned int i=0; i<alpha.size(); i++)
    {
      //Dans le tableau alpha, c'est à partir de nblig qu'on a alpha⁰, alpha¹, alpha²...alpha^degre-1
      //et le tab alpha commence à alpha^degre donc on va tenté de mettre tout ce qu'il après nblig dans le tab alpha
      //en premier élément
      if(deg+l+1>nblig)
	{
	  alphaord[j]=alpha[i];
	  
	  j++;
	}
      l++;
    }
    //Une les premiers éléments correctement placé dans alphaord
    //On copie le reste des éléments de alpha (aprés alpha^degre-1) dans alphaord.
    for(unsigned int k=deg; k<alphaord.size(); k++)
      {
	alphaord[k]=alpha[k-deg];
      } 
    //On recopie ensuite les éléments de alphaord dans la deuxième colonne
    //qui sont les valeurs décimales des équivalences
    //Et dans la première colonne du dico, la puissance associée
    for(unsigned int i=0; i<alphaord.size(); i++)
      {
	dico[i][0]=m;
	dico[i][1]=alphaord[m];
	m++;
      }
    //On parcourt syndr de la position 1 à la position 2*t pour bien avoir syndrome 1...syndrome 2*t en termes de positions
    //(C'est du détail mais c'est surtout pour faciliter la lecture lors du print des syndromes)
    for (unsigned int k=1; k<=syndr.size(); k++)
      {
	//On remet à 0 syndrun
	for(unsigned int i=0; i<syndrun.size(); i++)
	  {
	    syndrun[i]=0;
	  }
	//On parcourt notre message encodé avec les erreurs
	for(unsigned int i=0; i<messaerr.size(); i++)
	  {
	    //Dès qu'un bit est à un (présence d'un alpha^i dans l'expression)
	    if(messaerr[i]==1)
	      {
		//On applique l'opérateur "rond" c'est à dire que si l'on cherche le syndrome 3, on remplace
		//tous les alpha par alpha³, donc toutes les puissances sont multiplié par 3.
		p=k*i;
		//On applique mod(x^nblig-1)
		if (p>=nblig)
		  {
		    if((p%nblig==0)&&(p!=0))
		      {
			p=0;
		      }
		    while(p>nblig)
		      {
			p-=nblig;
		      }
		  }
		//On va chercher dans le vecteur alpha le nombre décimal correspondant à la puissance de alpha détécté
		//et on convertit la valeur en binaire 
		puiss=toBinary(alphaord[p], deg);
		for(unsigned int j=0; j<deg; j++)
		  {
		    //On l'ajoute syndrun qui vaut 0
		    syndrun[j]+=puiss[j];
		    //On tient compte du mod(2)
		    if (syndrun[j]==2)
		      {
			syndrun[j]=0;
		      }//if
		  }//for j
		}//if messaerr
	    //On continue de détecter les autres puissance de alpha de notre premier syndrome et les ajouter ensemble dans syndrun en tenant
	    //compte du mod(2), il est à préciser qu'on est obliger de convertir en bianire avant d'ajouter car le mod(2) ne s'applique pas aux
	    //valeure décimales, par conséquent cela donne des résultats erronnés.
	  }//for i
	
	//On repasse le binaire obtenus par l'ajout des équivalences de toutes les puissances détectés en valeur décimale
	unsigned int decimal=todec(syndrun,deg);
	
	//On s'assure que le syndrome soit pas nul, ce qui ne figure pas dans le dico, donc on traite ce cas ici
	if(decimal==0)
	  {
	    syndr[q]=0;
	  }
	//Puisqu'on stocke les puissances de alpha dans  notre matricesyndrome, il faut distinguer alpha¹ de 1
	//1=alpha⁰ mais le probleme reste le même pour distinguer alpha⁰ de 0
	//Donc on met une étiquette caracteristique pour aider le programme à comprendre que alpha¹!=1
	else if(decimal == 1)
	  {
	    syndr[q]=1111;
	  }
	//Sinon on cherche dans le dico à quelle puissance de alpha correspond ce nombre décimal.
	else
	  {
	    for(unsigned int i=0; i<alphaord.size(); i++)
	      {
		if(dico[i][1]==decimal)
		  {
		   syndr[q]=dico[i][0];
		  }
	      }//for i
	  }//else
	q++;//On passe au prochain syndrome qui sera stocké dans la case suivante de syndr
      }//for k
  }
  
  //Méthode qui va calculer le nombre d'erreur dans le message
  unsigned int findnberr(unsigned int degre)
  {
    unsigned int determinant;//Le determinant final de la matrice
    nberr=errcormax;//On stocke le nombre d'erreur corrigeable maximum 
    Matsyndrinit.resize(nberr, vector <unsigned int>(nberr));//On resize à la bonne dimension
    Matsyndr.resize(nberr, vector <unsigned int>(nberr));//On resize à la bonne dimension
    //On stocke nos syndromes
    for(unsigned int i=0; i<errcormax; i++)//colonne
      {
	unsigned int k=i;
	for(unsigned int j=0; j<errcormax; j++)//ligne
	  {
	    Matsyndrinit[j][i]=syndr[j+k];
	  }//for j
      }//for i
    //On sauvegarde notre matrice des syndromes pour l'afficher
    for(unsigned int i=0; i<errcormax; i++)//colonne
      {
	for(unsigned int j=0; j<errcormax; j++)//ligne
	  {
	    Matsyndr[j][i]=Matsyndrinit[j][i];
	  }//for j
      }//for i
    /* printf("Notre matrice des syndromes est :\n");
    for(unsigned int i=0; i<nberr; i++)
      {
	for (unsigned int j=0; j<nberr; j++)
	  {
	    printf(" %d", Matsyndr[j][i]);
	  }
	printf("\n");
      }
      printf("\n");*/
     
    //On calcule le determinant de notre matrice
    determinant= determinantOfMatrix(Matsyndr,nberr, degre);
    //printf("Determinant=%d\n\n", determinant);
    //Si le determinant vaut 0, c'est qu'il n'y a pas nberr erreur
    //On passe à la verification pour nberr-1
    if (determinant == 0)
      {
	//printf("Il n'y a pas %d erreurs\n\n", nberr);
       //On réitère le calcul du déterminant jusqu'à obtenir un determinant !=0
       //Indiquant le bon nombre d'erreur
       while((determinant==0)&&(nberr>0))
	 {
	   nberr--;
	   if(nberr==0)
	     {
	       return nberr;
	     }
	   else if((nberr==1)&&(Matsyndr[0][0]!=0))
	     {
	       return nberr;
	     }

	   else
	     {
	       Matsyndr.resize(nberr, vector <unsigned int>(nberr));
	       for(unsigned int i=0; i<nberr; i++)//colonne
		 {
		   unsigned int k=i;
		   for(unsigned int j=0; j<nberr; j++)//ligne
		     {
		       Matsyndr[j][i]=syndr[j+k];
		     }//for j
		 }//for i
	       determinant= determinantOfMatrix(Matsyndr,nberr, degre);
	       if (determinant != 0)
		 {
		   return nberr;
		 }
	     }//else
	 }//while
      }
    //Sinon c'est qu'il y a au moins nberr
    else 
      {
	//printf("Il y au moins %d erreurs\n", nberr);
	return nberr;
      }
    return nberr;
  }//Méthode

  //méthode qui determine les sous matrices d'une
  // EX les sous matrices de la matrice
  //|1 2 3| 
  //|4 5 6|
  //|7 8 9|
  // sont
  //|5 6| et |4 6| et |4 5|
  //|8 9|    |7 9|    |7 8|
 std::vector<std::vector<unsigned int> > subMatrix(std::vector<std::vector<unsigned int> > matrice,std::vector<std::vector<unsigned int> > temp, unsigned int p, unsigned int q, unsigned int size) {
   unsigned int i = 0, j = 0;
   // filling the sub matrix
   for (unsigned int row = 0; row < size; row++) {
      for (unsigned int col = 0; col < size; col++) {
         // skipping if the current row or column is not equal to the current
         // element row and column
	if ((row != p) && (col != q) ){
            temp[i][j++] = matrice[row][col];
            if (j == size - 1) {
               j = 0;
               i++;
            }//if j
	}//if row p col q
      }//for col      
   }//for row
   return temp;
}
  //Méthode qui calcule le determinant d'une matrice
  unsigned int determinantOfMatrix(std::vector<std::vector<unsigned int> > matrix, unsigned int size, unsigned int deg) {
    //On sauvegarde la taille initiale de notre matrice syndrome 
  unsigned int sizevect=deg;
  unsigned int determinant = 0;
  bool flagz=false;
  bool flago=false;
  //Pterm et sterm correpondent aux 2 produits de coefficient des sous matrices 2x2
  // EX de la matrice
  //|1 2 3| 
  //|4 5 6|
  //|7 8 9|
  // Les sous matrices sont :
  //|5 6| et |4 6| et |4 5|
  //|8 9|    |7 9|    |7 8|
  //pour la premiere sous matrice pterm = 45 qui est 5*9
  // et sterm = 8*6 = 48
  //En fait notre matrice stockes les puissance dans les coefficients
  //Comme on ne peut pas stocker des puissances de alpha directement
  //On doit pas multiplier les coeffients entre eux mais les ajouter
  //car alpha⁴*alpha⁵ != 5*4
  //Donc on fait 5+4 = à la puissance de alpha⁹
  unsigned int pterm=160;
  unsigned int sterm=160;
  //Addterm est l'addition de pterm et sterm
  unsigned int addterm;
  //Conversion décimale des puissance pterm et sterm avant l'addition via le dictionnaire 
  unsigned int ptermdec;
  unsigned int stermdec;
  //Spterm corresponds aux coefficient multiplicatifs de nos sous matrices
  // EX de la matrice
  //|1 2 3| 
  //|4 5 6|
  //|7 8 9|
  // Les sous matrices sont :
  //Spterm 1 * |5 6| et spterm2*|4 6| et spterm3*|4 5|
  //           |8 9|            |7 9|            |7 8|
  //spterm1 = 1, spterm2 = 2 , spterm3 = 3
  unsigned int spterm;
  //ssterm contiendra le determinant de nos sous matrices
  unsigned int ssterm;
  //Conversion en  binaire de nos puissance pterm et sterm
  std::vector<unsigned int> ptermbin;
  std::vector<unsigned int> stermbin;
  //Addition binaire de nos ptermbin et stermbin
  std::vector<unsigned int> detbin;
  //Conversion binaire de spterm et ssterm
  std::vector<unsigned int> sptermbin;
  std::vector<unsigned int> sstermbin;
  std::vector<unsigned int> sdetbin;
  //Contiendra tous nos determinants des sous matrices avant de les sommer
  std::vector<unsigned int> listdetinter;
  //Conversion en binaire pour ajout 
  std::vector<unsigned int> listdetinterbin;
  //Vecteur binaire contenant les ajouts des determinants intermédiaire
  std::vector<unsigned int> somdeterminant(sizevect,0);

  //Resize nos vecteurs
  ptermbin.resize(sizevect);
  stermbin.resize(sizevect);
  detbin.resize(sizevect);
  sptermbin.resize(sizevect);
  sstermbin.resize(sizevect);
  sdetbin.resize(sizevect);
  listdetinter.resize(factorial(errcormax+1));
  listdetinterbin.resize(sizevect);
  somdeterminant.resize(sizevect);
  
  //Si on a une matrice à un coefficient le determinant est la valeur de son unique coefficient
  if (size == 1) {
    return matrix[0][0];
  }
  //Si la taille de la matrice est 2
  if (size == 2) {
    //On controle les valeurs des 4 cases de cette matrices
    /*if((matrix[size-2][size-2]==0)||(matrix[size-1][size-1]==0))
      {
	pterm=0;
      }
    else if((matrix[size-1][size-2]==0)||(matrix[size-2][size-1]==0))
      {
	sterm=0;
      }
    else
    {*/
	switch(matrix[size-2][size-2])
	  {
	    //Si la valeur vaut 0 alors pterm = 0*matrix[size-1][size-1] =0
	  case 0:
	    pterm=0;
	    //printf("dans case 0 pterm size-2 size-2\n");	
	    break;
	    //Si la valeur vaut 1111 = 1 alors pterm = 1*matrix[size-1][size-1]
	  case 1111:
	    pterm=matrix[size-1][size-1];
	    flago=true;
	    //printf("dans case 1111 pterm size-2 size-2\n");
	    break;
	    //Sinon si elle vaut tout autre valeur 
	  default:
	    //On analyse la valeur de pterm pour savoir si elle n'a pas déjà été calculé dans un autre switch case
	    switch(pterm)
	      {
		//Si elle est nulle, elle le reste
	      case 0:
		pterm=0;
		break;
		//Si elle vaut 1 alors on est déjà rentré dans un case où matrix[size-1][size-1] = 1
	      case 1111:
		pterm= matrix[size-2][size-2];
		break;
		//Sinon elle vaut encore 160 (valeur initiale) alors pterm est encore vierge et non calculé
		//Donc on peut ajouter le coefficient non nulle et non unitaire de matrix[size-2][size-2] (qui est la puissance de alpha)
	      case 160:
		//printf("cas défaut\n");
		pterm=0;
		pterm+=matrix[size-2][size-2];
		break;
		//Sinon c'est que pterm n'est pas vierge et que matrix[size-2][size-2] est non unitaire
		//En fait la condition if est la au cas ou matrix[size-2][size-2]=1111 et que pterm = 1111
		//Si matrix[size-2][size-2] on rentre dans le cas ou pterm=matrix[size-1][size-1] mais si matrix[size-1][size-1]==1111
		//alors pterm = 1111 et donc on va tomber dans le cas pterm=matrix[size-2][size-2]=1111
		//Donc pour éviter d'ajouter deux fois 1111 qui donnerait 0 on vérifie que pterm soit pas déjà unitaire pour ajouter matrix[size-2][size-2]
	      default:
		if((pterm != 1111)&&(flago==false))
		  {
		    pterm+=matrix[size-2][size-2];
		    //printf("ajout1\n\n");
		  }
		else
		  {
		    flago=false;
		  }
		break;
	    
	      }	
	  }
	//Idem que précédemment
	switch(matrix[size-1][size-1])
	  {
	  case 0:
	    pterm=0;
	    //printf("dans case 0 pterm size-1 size-1\n");
	    break;
	  case 1111:
	    pterm=matrix[size-2][size-2];
	    flago=true;
	    // printf("dans case 1111 pterm size-1 size-1\n");
	    break;
	  default:
	    switch(pterm)
	      {
	      case 0:
		pterm=0;
		break;
	      case 1111:
		pterm= matrix[size-1][size-1];
		break;
	      case 160:
		//printf("cas défaut\n");
		pterm=0;
		pterm+=matrix[size-1][size-1];
		break;
	      default:
		if((pterm != 1111)&&(flago==false))
		  {
		    pterm+=matrix[size-1][size-1];
		    //printf("ajout2\n\n");
		  }
		else
		  {
		    flago=true;
		  }
		break;   
	      }
	  }
	switch(matrix[size-1][size-2])
	  {
	  case 0:
	    sterm=0;
	    //printf("dans case 0 sterm size-1 size-2\n");
	    break;
	  case 1111:
	    sterm=matrix[size-2][size-1];
	    flagz=true;
	    //printf("dans case 1111 sterm size-1 size-2\n");
	    break;
	  default:
	    switch(sterm)
	      {
	      case 0:
		sterm=0;
		break;
	      case 1111:
		sterm= matrix[size-1][size-2];
		break;
	      case 160:
		//printf("cas défaut\n");
		sterm=0;
		sterm+=matrix[size-1][size-2];
		break;
	      default:
		if((sterm != 1111)&&(flagz==false))
		  {
		    sterm+=matrix[size-1][size-2];
		    //printf("ajout3\n\n");
		  }
		else
		  {
		    flagz=false;
		  }
		break;	    
	      }
	  }
	switch(matrix[size-2][size-1])
	  {
	  case 0:
	    sterm=0;
	    //printf("dans case 0 sterm size-2 size-1\n");
	    break;
	  case 1111:
	    sterm=matrix[size-1][size-2];
	    flagz=true;
	    //printf("dans case 1111 sterm size-2 size-1\n");
	    break;
	  default:
	    switch(sterm)
	      {
	      case 0:
		sterm=0;
		break;
	      case 1111:
		sterm= matrix[size-2][size-1];
		break;
	      case 160:
		//printf("cas défaut\n");
		sterm=0;
		sterm+=matrix[size-2][size-1];
		break;
	      default:
		if((sterm !=1111)&&(flagz==false))
		  {
		    sterm+=matrix[size-2][size-1];
		    //printf("ajout4\n\n");
		  }
		else
		  {
		    flagz=false;
		  }	    
		break;   
	      }
	  }
	// }
      //Si 
	if(((matrix[size-2][size-2]!=0)&&(matrix[size-2][size-2]!=1111))&&((matrix[size-1][size-1]!=0)&&(matrix[size-1][size-1]!=1111))&&((matrix[size-1][size-2]!=0)&&(matrix[size-1][size-2]!=1111))&&((matrix[size-2][size-1]!=0)&&(matrix[size-2][size-1]!=1111)))
	{
	  pterm=matrix[size-2][size-2] + matrix[size-1][size-1];
	  sterm=matrix[size-2][size-1] + matrix[size-1][size-2];
	}
      //printf("Les derniers termes sont\n matrix[%d][%d]=%d\nmatrix[1][1]=%d\nmatrix[0][1]=%d\nmatrix[1][0]=%d\n",size-2, size-2, matrix[size-2][size-2],matrix[size-1][size-1],matrix[size-2][size-1],matrix[size-1][size-2]);
      //printf("pterm avant=%d\n", pterm);
      //printf("sterm avant=%d\n", sterm);
      //On applique le modulo correspondant au degré
      if ((pterm>=nblig)||(sterm>=nblig))
      {
	//Si la puissance de alpha est multiple de (2^deg)-1
	//Alors elle vaut 1 modulo
	//On est obligé d'inclure pterm != 0 car 0 divisé par n'importe quel nombre (hormis 0)
	// donne un reste nul car le if précédent est un ou, donc pterm peut valoir 0 et sterm >= nblig...
	if(((pterm%nblig==0)&&(pterm!=0))||(pterm==1111))
	  {
	    //printf("je rentre 1\n");
	    pterm=1111;
	  }
	else
	  {
	    while(pterm>nblig)
	      {
		pterm-=nblig;
	      }
	  }
	if(((sterm%nblig==0)&&(sterm!=0))||(sterm==1111))
	  {
	    //printf("je rentre 2\n");
	    sterm=1111;
	  }
	else
	  {
	    while(sterm>nblig)
	      {
		sterm-=nblig;
	      }
	  }
      }//if pterm sterm
      //printf("pterm apres=%d\n", pterm);
      //printf("sterm apres=%d\n", sterm);
      //On traite le cas ou pterm est nul
      if((pterm==0)&&(sterm!=0))
       {
	 //ptermbin va valoir 0
	  ptermbin=toBinary(0,sizevect);
	  //On s'assure qu'sterm soit ou non unitaire
	  if(sterm==1111)
	    {
	      stermbin=toBinary(1,sizevect);
	    }
	  else
	    {
	      //Si sterm n'est pas unitaire alors on peut sans problème aller chercher sa valeur décimale dans le dictionnaire
	      for(unsigned int i=0; i<alphaord.size(); i++)
		{
		  if(dico[i][0]==sterm)
		    {
		      stermdec=dico[i][1];
		      //printf("stermdec=%d\n\n", stermdec);
		    } 
		}//for i
	      //Convertir cette valeur en binaire qui donnera l'équivalent en puissance de alpha "exploitable"
	      stermbin=toBinary(stermdec,sizevect);
	    }
       }
      //On traite de la même le cas inverse ou seulement sterm est nul 
      else if((pterm!=0)&&(sterm==0))
	{
	  if(pterm==1111)
	    {
	      ptermbin=toBinary(1,sizevect);
	    }
	  else
	    {
	      for(unsigned int i=0; i<alphaord.size(); i++)
		{
		  if(dico[i][0]==pterm)
		    {
		      ptermdec=dico[i][1];
		      //printf("ptermdec=%d\n\n", ptermdec);
		    }
		}//for i
	      ptermbin=toBinary(ptermdec,sizevect);
	    }
	  stermbin=toBinary(0,sizevect);
	}
      //On traite le cas ou pterm est unitaire
      else if((pterm==1111)&&(sterm!=1111))
	{
	  ptermbin=toBinary(1,sizevect);
	  for(unsigned int i=0; i<alphaord.size(); i++)
	    {
	      if(dico[i][0]==sterm)
		{
		  stermdec=dico[i][1];
		  //printf("stermdec=%d\n\n", stermdec);
		} 
	    }//for i
	  stermbin=toBinary(stermdec,sizevect);
	}
      else if((sterm==1111)&&(pterm!=1111))
	{
	  for(unsigned int i=0; i<alphaord.size(); i++)
	    {
	      if(dico[i][0]==pterm)
		{
		  ptermdec=dico[i][1];
		  //printf("ptermdec=%d\n\n", ptermdec);
		}
	    }//for i
	  ptermbin=toBinary(ptermdec,sizevect);
	  stermbin=toBinary(1,sizevect);
	}
      else if((sterm==1111)&&(pterm==1111))
	{
	  ptermbin=toBinary(1,sizevect);
	  stermbin=toBinary(1,sizevect);
	}
      else if((pterm==0)&&(sterm==0))
	{
	  ptermbin=toBinary(0,sizevect);
	  stermbin=toBinary(0,sizevect);
	}
      //Si sterm et pterm sont non null et non unitaire
      //On peut aller chercher leurs valeurs respectives dans le dans le dico
      else
	{
	  for(unsigned int i=0; i<alphaord.size(); i++)
	    {
	      if(dico[i][0]==pterm)
		{
		  ptermdec=dico[i][1];
		  //printf("ptermdec=%d\n\n", ptermdec);
		}
	      if(dico[i][0]==sterm)
		{
		  stermdec=dico[i][1];
		  //printf("stermdec=%d\n\n", stermdec);
		} 
	    }//for i
	  ptermbin=toBinary(ptermdec,sizevect);
	  stermbin=toBinary(stermdec,sizevect);
	}
      //On additionne les valeurs binaire de pterm et sterm
      for(unsigned int i=0; i<sizevect; i++)
	{
	  detbin[i]=ptermbin[i]+stermbin[i];
	  if(detbin[i]==2)
	    {
	      detbin[i]=0;
	    }
	}
      addterm=todec(detbin,sizevect);
      //pour ne pas confondre addterm=0 de addterm=1, qui donnent tous les deux lieux au même déterminant numérique
      //il faut distinguer 0 de alpha⁰
      //On définis donc 1112 qui est notre étiquette pour un détérminant nul
      if(addterm==0)
	{
	  determinant=1112;
	}
      else
	{
	  //On parcourt le dictionnaire en veillant à bien interpréter alpha⁰=1
	  for(unsigned int i=0; i<alphaord.size(); i++)
	    {
	      if(dico[i][1]==addterm)
		{
		  determinant=dico[i][0];

		}
	      if (determinant == 0)
		{
		  determinant = 1111;
		}
	    }
	}
      //Si le determinant vaut 1112 alors on l'interprète comme 0 et non alpha⁰
      if(determinant==1112)
	{
	  determinant=0;
	}
      return determinant;
  }//fin if(size==2)
  
  std::vector<std::vector<unsigned int> > temp;
  temp.resize(size, vector <unsigned int>(size));
  std::vector<std::vector<unsigned int> > ntemp;
  ntemp.resize(size, vector <unsigned int>(size));
  // On va calculer autant de sous determinants que la taille de la matrice 
  for (unsigned int i = 0; i < size; i++) {
    ntemp=subMatrix(matrix, temp, 0, i, size);
    /*for(unsigned int t=0; t<ntemp.size(); t++)
      {
	for(unsigned int v=0; v<ntemp.size(); v++)
	  {
	    printf(" %d", ntemp[t][v]);
	  }
	printf("\n");
      }
      printf("\n");*/
    //On stocke le premier spterm de la sous matrice ntemp
    spterm=matrix[0][i];
    //printf("spterm=%d\n\n", spterm);
    //On recalcule le determinant de la première sous matrice ntemp en réàactualisant la size
    //ssterm stockera donc le determinant de la sous matrice associé à son coefficient multiplicatif (spterm)
    ssterm= determinantOfMatrix(ntemp, size - 1, deg);
    //printf("ssterm=%d\n\n", ssterm);
    //On controle les cas particulier unitaire et nulle
    if((ssterm==0)||(spterm==0))
      {
	determinant=0;
      }
    else if(spterm==1111)
      {
	determinant=ssterm;
      }
    else if(ssterm==1111)
      {
	determinant=spterm;
      }
    //S'il n'y a pas de résultats unitaire ou nul
    //on additionne nos deux coeffcients 
    else
      {
	determinant=ssterm+spterm;
      }
    //On applique le modulo ((x^deg)-1)
    if (determinant>=nblig)
      {
	if(determinant%nblig==0)
	  {
	    determinant=1111;
	  }
	if(determinant!=1111)
	  {
	    while(determinant>nblig)
	      {
		determinant-=nblig;
	      }
	  }
      }//if pterm sterm
    //printf("Le determinant en dehors des if est alpha puissance = %d\n\n", determinant);
    //On stocke les determinants intermédiaires calculés 
    listdetinter[u]=determinant;
    u++;
  }
  /* printf("La liste des determinants est \n");
  for(unsigned int s=0; s<listdetinter.size(); s++)
    {
      if(listdetinter[s]!=0)
	{
	  printf("%d ", listdetinter[s]);
	}
    }
    printf("\n\n");*/
  //On additionne nos determinants intermédiaires
  for(unsigned int v=0; v<listdetinter.size();v++)
    {
      if(listdetinter[v]==0)
	{
	  listdetinterbin=toBinary(0,sizevect);
	  for(unsigned int w=0; w<somdeterminant.size(); w++)
	    {
	      somdeterminant[w]+=listdetinterbin[w];
	      if(somdeterminant[w]==2)
		{
		  somdeterminant[w]=0;
		}
	    }//for w
	}//if
      else if (listdetinter[v]==1111)
	{
	  listdetinterbin=toBinary(1,sizevect);
	  for(unsigned int w=0; w<somdeterminant.size(); w++)
	    {
	      somdeterminant[w]+=listdetinterbin[w];
	      if(somdeterminant[w]==2)
		{
		  somdeterminant[w]=0;
		}
	    }//for w
	}
      else
	{
	  for(unsigned int x=0; x<alphaord.size(); x++)
	    {
	      if(dico[x][0]==listdetinter[v])
		{
		  listdetinterbin=toBinary(dico[x][1], sizevect);
		  break;
		}
	    }//for i
	  for(unsigned int w=0; w<somdeterminant.size(); w++)
	    {
	      somdeterminant[w]+=listdetinterbin[w];
	      if(somdeterminant[w]==2)
		{
		  somdeterminant[w]=0;
		}
	    }
	   
	}//else
    }//for v
  /*printf("La somme des determinant en binaire est ");
  for(unsigned int y=0; y<somdeterminant.size(); y++)
    {
      printf("%d", somdeterminant[y]);
    }
    printf("\n\n");*/
  determinant=todec(somdeterminant, sizevect);
  //printf("Le determinant final en decimal est égale à %d\n", determinant);
  if(determinant!=0)
    {
      for(unsigned int y=0; y<alphaord.size(); y++)
	{
	  if(dico[y][1]==determinant)
	    {
	      if(y==0)
		{
		  determinant=1111;
		}
	      else
		{
		  determinant=dico[y][0];
		  break;
		}
	    }//if dico
	}//for y
    }//if determinant
  //printf("Le determinant final est égale à %d\n", determinant);
  return determinant;
  }
  
  void coeff(unsigned int degre)
  {
    std::vector<std::vector<unsigned int> > Matsyndrinter;//Notre matrice des syndromes intermédiare
    Matsyndrinter.resize(nberr, vector <unsigned int>(nberr));
    unsigned int det;
    unsigned int determinantf;
    unsigned int val=0;
    std::vector<std::vector<unsigned int> > termes;//Notre matrice des syndromes intermédiare
    std::vector<unsigned int> termesbinu;//Notre matrice des syndromes intermédiare
  
    //unsigned int poserreur;    
     
    termes.resize(nblig, vector <unsigned int>(degre));
    termesbinu.resize(degre);
    sigma.resize(nberr);
    //Mettre en place les matrices nécessaire pour le calcul des t coefficients
    //du polynome localisateur d'erreurs
    for (unsigned int i=0; i<nberr; i++)
      {  
	for (unsigned int j=0; j<nberr; j++)
	  {
	    if(i==j)
	      {
		for (unsigned int k=0; k<nberr; k++)
		  {
		    Matsyndrinter[k][j]=syndr[nberr+k];
		  }
	      }
	    else
	      {
		for (unsigned int k=0; k<nberr; k++)
		  {
		    Matsyndrinter[k][j]=Matsyndrinit[j][k];
		  }
	      }
	  }
	for (unsigned int t=0; t<nblig; t++)
	  {
	    for (unsigned int u=0; u<degre; u++)
	      {
		termes[t][u]=0;
	      }
	  }
	/*printf("La matrice des syndromes pour le coefficient %d est \n", i+1);
	for (unsigned int l=0; l<nberr; l++)
	  {  
	    for (unsigned int m=0; m<nberr; m++)
	      {
		printf(" %d", Matsyndrinter[l][m]);    
	      }
	    printf("\n");	 
	  }
	  printf("\n");*/
	det= determinantOfMatrix(Matsyndrinter, nberr, degre);
	//printf("Le determinant %d de la matrice associé est %d\n", i+1, det);
	determinantf= determinantOfMatrix(Matsyndrinit, nberr, degre);
	//printf("Le determinant de la matrice initial est %d\n", determinantf);  
	if((determinantf>det)&&(determinantf!=1111))
	  { 
	    sigma[nberr-i]=det-determinantf;
	    sigma[nberr-i]+=nblig;
	  }
	else if((determinantf<det)&&(det!=1111))
	  {
	    sigma[nberr-i]=det-determinantf;
	  }
	else if(determinantf==1111)
	  {
	    sigma[nberr-i]=det;
	  }
	else if((det==1111)&&(determinantf!=1111))
	  {
	    sigma[nberr-i]=-determinantf+nblig;
	  }
	else
	  {
	    sigma[nberr-i]=1111;
	  }
	//intf("Le sigma %d de la matrice associé est %d\n", nberr-i, sigma[nberr-i]);

      }//for i
    //trouver la position de l'erreur avec la recherche de Chien
    if(nberr==1)
      {
	printf("L'erreur est à la position %d\n\n\n", sigma[nberr]);
      }
    else
      {
	for(unsigned int j=0; j<nblig; j++)
	  {
	    for(unsigned int i=0; i<=nberr; i++)
	      {
		if(i==nberr)
		  {
		    if(j==0)
		      {
			val=1111;
		      }
		    else
		      {
			val=i*dico[j][0];
		      }
		    // printf(" val dans nberr = %d\n", val);
		    if(val==1111)
		      {
			termesbinu=toBinary(1, degre);
		      }
		    else
		      {
			while(val>nblig)
			  {
			    val-=nblig;
			  }
			if((val%nblig==0)&&(val!=0))
			  {
			    termesbinu=toBinary(1, degre);
			  }
			else if(val==0)
			  {
			    termesbinu=toBinary(0, degre);
			  }
			else
			  {
			    termesbinu=toBinary(dico[val][1], degre);
			  }
		      }
		    /* printf("Le termes en binaire est égale à \n");
		    for(unsigned int z=0; z<degre; z++)
		      {
			printf("%d", termesbinu[z]);
		      }
		      printf("\n");*/
		  }
		else
		  {
		    if(sigma[nberr-i]==1111)
		      {
			val=dico[j][0]*i;
		      }
		    else 
		      {
			val=dico[j][0]*i+sigma[nberr-i];
		      }
		    // printf(" val ailleurs = %d\n", val);
		    while(val>nblig)
		      {
			val-=nblig;
		      }
		    if((val%nblig==0)&&(val!=0))
		      {
			val=1;
		      }
		    if(val==0)
		      {
			termesbinu=toBinary(0, degre);
		      }
		    else
		      {
			termesbinu=toBinary(dico[val][1], degre);
		      }
		    /* printf("Le termes en binaire est égale à \n");
		    for(unsigned int z=0; z<degre; z++)
		      {
			printf("%d", termesbinu[z]);
		      }
		      printf("\n");*/
		  }
		for(unsigned int k=0; k<degre; k++)
		  {
		    termes[j][k]+=termesbinu[k];
		    if(termes[j][k]==2)
		      {
			termes[j][k]=0;
		      }
		  }
		/*printf("L'addition en binaire est égale à \n");
		for(unsigned int z=0; z<degre; z++)
		  {
		    printf("%d", termes[j][z]);
		  }
		  printf("\n");*/
	      }//for i
	    /*if((todec(termes[j],degre))!=0)
	      {
		printf("aucune erreur en position %d\n\n\n", j);
		}*///if dico
	   if((todec(termes[j],degre))==0)
	      {
		//printf("il y a erreur en position %d\n\n\n", j);
		allerr++; 
	      }
	  }//for j
	/*if(allerr!=nberr)
	  {
	    printf("Le message n'est pas corrigeable\n\n\n");
	  }
	else
	  {
	    printf("Toutes les erreurs ont été trouvée(s)\n\n\n");
	    }*/
      }
  }
  
  //Méthode pour convertir un nombre décimal en vecteur binaire
  std:: vector <unsigned int> toBinary(int nombre, int n)
  {
    std:: vector <unsigned int> r;
    r.resize(n);
    unsigned int i=0;
    unsigned int reste=0;
    
    if (nombre==1)
      {
	for(unsigned int i=0; i<r.size()-1; i++)
	  {
	    r[i]=0;
	    }
	r[r.size()-1]=1;
	}
    else
      {
	while (nombre != 0){
	  reste=nombre%2;
	  if(reste==0)
	    {
	      r[r.size()-1-i] =0;
	    }
	  else
	    {
	      r[r.size()-1-i]=1;
	    }
	  nombre /= 2;
	  i++;
	}
      }
    return r;
  }

  //Méthode pour convertir un vecteur binaire en décimal
  unsigned int todec(std:: vector<unsigned int> syndrome, int n)
  {
    unsigned int decimal=0;
    syndrome.resize(n);
    for (unsigned int i=0; i<syndrome.size(); i++)
      {
	if(syndrome[syndrome.size()-1-i]==1)
	  {
	    decimal += pow(2, i);
	  }
      }
    return decimal;
  }
  unsigned int factorial(unsigned int n)
  {
    unsigned int ret = 1;
    while (n > 1)
        ret *= n--;
    return ret;
  }
 
  void afficher()
  {
     //On print les éléments du dico
    printf("\nLes valeurs des puissances de alpha sont :\n");
    for(unsigned int i=0; i<alphaord.size(); i++)
    {
      printf("Pour alpha puissance %d la valeur décimale du polynome associé est %d\n", dico[i][0], dico[i][1]);
    }
    printf("\n");
    //On print nos syndromes
    for(unsigned int i=0; i<syndr.size(); i++)
    {
      if(syndr[i]==0)
	{
	   printf("Le syndrome %d est ", i+1);
	   printf("%d\n", 0);
	}
      else if(syndr[i]==1111)
	{
	  printf("Le syndrome %d est ", i+1);
	  printf("%d\n", 1);
	}	
      else
	{
	  printf("Le syndrome %d est alpha puissance ", i+1);
	  printf("%d\n", syndr[i]);
	}
    }
    //printf("alphor=%ld\n", alphaord.size());
    printf("\n\n");
    printf("Notre matrice des syndromes est :\n");
    for(unsigned int i=0; i<errcormax; i++)
    {  
      for (unsigned int j=0; j<errcormax; j++)
	{	  
	  printf(" %d", Matsyndrinit[j][i]);
	}
       printf("\n");
    }
    printf("\n");
    if(nberr==errcormax)
      {
	printf("Il y a au moins %d erreurs\n", nberr);
      }
    else
      {
	printf("Il y a %d erreurs\n", nberr);
      }
    printf("Il y avait %d erreurs\n", nberreel);
    for(unsigned int i=0; i<nberr; i++)
      {
	printf("Le coefficient %d du polynôme localisateur d'erreur est %d\n", nberr-i, sigma[nberr-i]);
      }
    printf("La recherche de Chien donne %d racine(s)\n", allerr);
  }
  ~decodeBCH()
  {
    for(unsigned int i=0; i<alphaord.size(); i++)
      {
	delete[] dico[i];
      }
    delete[] dico;
  }
};//Fin de classe 

#endif
