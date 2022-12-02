# Stage_5A_Simpulse
[![forthebadge](http://forthebadge.com/images/badges/built-with-love.svg)](http://forthebadge.com)  [![forthebadge](http://forthebadge.com/images/badges/powered-by-electricity.svg)](http://forthebadge.com)

Ce projet est le codage en C++ de l'algorithme de PGZ.

## Contenus du code
Le code permet de générer un mot de longueur m aléatoirement avec tiragede.h.
Il permet dans un second temps d'encoder en BCH, par multiplication de g(x), le mot binaire avec encodeBCH.h.
Ensuite, de bruiter ce mots avec ajouterr.h
Enfin, de décoder le mot pour compter le nombre d'erreur, et ensuite trouver les racines (decodeBCH.h).

### Pré-requis
Le décodeur PGZ se divise en 2 parties : Comptage du nombre d'erreurs et localisation des erreurs. Contrairement à BM qui combine les 2 étapes. De plus, l'algorithme n'est pas totalement au point, en particulier sur la recherche de Chien. C'est à peu près au moment de l'implémenter que je me suis rendu compte que cela n'irait pas avec les applications Simpulse.

### Installation


Il suffit simplement de recopier le répertoire en local.

## Démarrage

Lancer le main avec les commandes g++ -Wall main.cpp -o opt et ensuite ./opt comme suit :

![30 11 2022_01 43 52_REC](https://user-images.githubusercontent.com/87069145/204685720-ffd633fb-564a-411d-b948-59e93d9d8cee.gif)

## Auteurs
* **Nouri Huseynov** _alias_ [@chucknouri](https://github.com/chucknouri)
