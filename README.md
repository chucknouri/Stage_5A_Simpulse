# Stage_5A_Simpulse
[![forthebadge](http://forthebadge.com/images/badges/built-with-love.svg)](http://forthebadge.com)  [![forthebadge](http://forthebadge.com/images/badges/powered-by-electricity.svg)](http://forthebadge.com)

Ce projet est le codage en C++ de l'algorithme de PGZ.

## Contenus du code
Le code permet de générer un mot de longueur m aléatoirement avec tiragede.h.
Il permet dans un second temps d'encoder en BCH, par multiplication de g(x), le mot binaire m avec encodeBCH.h.
Ensuite, de bruiter ce mots avec ajouterr.h
Puis de trouver les polynômes minimaux pour générer nos équivalences entre les éléments primitifs 𝜶^i et les puissance (𝜶¹, 𝜶²,...,𝜶^(m-1)).
Une fois ces équivalence obtenus les opérations implémenté permettent le calculs des syndrômes d'erreurs.
Enfin, ce code va permettre de décoder le mot pour compter le nombre d'erreur, et ensuite trouver les racines (decodeBCH.h).

## Problèmes du code
Tout d'abord, la méthode d'encodage implémentée n'est pas celle préconisé par l'ETSI. Il s'agit de la multiplication par g(x) qui présente l'inconvénient de ne pas clairement distinguer BBFrame et BCHFEC. De plus, l'algorithme de décodage n'est pas non plus bon puisqu'il s'agit de l'algorithme de Peterson-Gorenstein-Ziegler. Cette algorithme travaille avec des déterminants de matrice (t*t) pour obtenir le nombre et les positions des erreurs. Or, pour des tailles t très grandes, le calcul n'est pas viable en termes de temps et complexité. 

### Pré-requis
Le décodeur PGZ se divise en 2 parties : Comptage du nombre d'erreurs et localisation des erreurs. Contrairement à BM qui combine les 2 étapes. De plus, l'algorithme n'est pas totalement au point, en particulier sur la recherche de Chien. C'est à peu près au moment de l'implémenter que je me suis rendu compte que le PGZ n'irait pas avec les applications Simpulse.

### Installation


Il suffit simplement de recopier le répertoire en local.

## Démarrage

Lancer le main avec les commandes g++ -Wall main.cpp -o opt et ensuite ./opt comme suit :

![30 11 2022_01 43 52_REC](https://user-images.githubusercontent.com/87069145/204685720-ffd633fb-564a-411d-b948-59e93d9d8cee.gif)

## Auteurs
* **Nouri Huseynov** _alias_ [@chucknouri](https://github.com/chucknouri)
