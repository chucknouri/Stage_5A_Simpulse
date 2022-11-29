#include "tiragede.h"
#include "encodeBCH.h"
#include "ajouterr.h"
#include "calcalpha.h"
#include "decodeBCH.h"
#include <stdio.h>
#include <stdlib.h>

static unsigned int degre=7; //Degré du message (sous forme polynomial) à envoyer

std:: vector<unsigned int> testtirage()
{
  Tiragede Ip(degre+1);//Instanciation de l'objet Ip
  vector<unsigned int> Ipp=Ip.gettirage(degre+1);//Stocke le tirage dans le vecteur Ipp pour l'afficher (impossible de faire Ip.afficher(Ip.gettiragede..) Ca va relancer un nouveau tirage.
  Ip.afficher(Ipp);
  return Ipp;
}
std:: vector<unsigned int> encodeBCH(std:: vector<unsigned int> vectencode)
{
  //On va compter le nombre de 0 à la postion des coefficient dominants dans le mot binaire à encoder 
  //Ex dans "010001" il y a 1 zero pour le coefficient dominant donc le polynome associé doit être x⁴+0*x³+0*x²+x+1
  //et non 0*x⁵+x⁴+0*x³+0*x²+x+1, sinon ça fause le calcul de l'encodage
   unsigned int s=0;
   while(vectencode[s]==0)
   {
     s++;
   }
   //Une fois qu'on a notre nombre de 0
   if(s==0)
     //Si s vaut 0 alors le nombre du coefficient dominant (de plus haut degré) vaut 1
     {
       //On ajuste le degre du polynome
       degre=vectencode.size()-1;
     }
   else
     {
       //On va décaler le vecteur d'autant de position que de 0 dominant
       //Ex 010010 va donner 100100
       for(unsigned int i=0; i<vectencode.size()-s; i++)
	 {
	   vectencode[i]=vectencode[i+s];
	 }
       for(unsigned int i=vectencode.size()-s; i<vectencode.size()-1; i++)
	 {
	   vectencode[i]=0;
	 }
       //On ajuste le degré et la taille du vecteur
       degre=vectencode.size()-1-s;
       vectencode.resize(degre+1);
     }
  std:: vector<unsigned int> messageencode;
  Encodebch enc(vectencode, degre);
  messageencode = enc.encodeBCH(vectencode);
  enc.afficher(messageencode);
  return messageencode;
}
std:: vector<unsigned int> ajouterre(std:: vector<unsigned int> vectencode)
{
  std:: vector<unsigned int> messagerre;
  AjoutErr err;
  messagerre = err.ajouterr(vectencode,errcormax);
  err.afficher(messagerre);
  return messagerre;
}
void decodebch(std:: vector<unsigned int> vecterr)
{
  std:: vector <unsigned int> tabalpha;
  Alpha alpha;
  alpha.toutpoly(degre);
  alpha.polyprim();
  tabalpha=alpha.polymin();
  alpha.afficher();
  decodeBCH dec;
  dec.syndromes(vecterr, tabalpha, degre);
  dec.findnberr(degre);
  dec.coeff(degre);
  dec.afficher();
 

}
int main()
{
  vector<unsigned int> I(degre+1);
  vector<unsigned int> encode;
  vector<unsigned int> encoderr;
  vector <unsigned int> mess;
  mess.resize(degre+1);
  //On génére un message m=4bits
  mess[0]=0;
  mess[1]=0;
  mess[2]=0;
  mess[3]=1;
  mess[4]=0;
  mess[5]=0;
  mess[6]=1;
  mess[7]=0;
  //Ou on peut générer un message aléatoire
  //(mettre la variable degré à 4) et remplacer mess
  //par I
  
  //I=testtirage(); 	
  encode=encodeBCH(mess);
  encoderr=ajouterre(encode);
  decodebch(encoderr);
  
  return 0;
}



