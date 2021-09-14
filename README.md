# Repliement-prot-MC
Repliement d'un modèle simplifié de protéine par un algorithme de Monte-Carlo

## Environnement

Tout d'abord, il faut installer l'environnement de travail sur votre machine grâce à la commande suivante entrée dans votre invité de commandes : 

 > $ **conda env create --file env-projet.yml**

Puis, pour activer cet environnement de travail, tapez la commande : 

 > $ **conda activate env-projet**

## Lancement du programme

Pour lancer le programme, placez vous dans le bon répertoire où se situe le programme python, puis entrez la commande suivante sur votre invité de commande :

 > $ **python projet.py sequence_proteine temperature**

avec :
projet.py : le nom du programme python du projet
sequence_proteine : la séquence protéique dont vous souhaitez connaître son repliement, en code à une lettre en majuscule
temperature : la température du repliement que vous souhaitez

## Exemple de lancement

 > $ **python projet.py TGLMFWGVYGATT 0.3**

En sortie sur votre invité de commande, vous obtiendrez toutes les conformations de protéine réalisées codée en ASCII. De plus, une fenêtre s'ouvrira avec le graphique de la variation des énergies des conformations de protéine.

## Aide

Si vous avez besoin d'aide, tapez la commande suivante :

 > $ **python projet.py -h**

Vous obtiendez une aide pour renseigner correctement les arguments de séquence et de température.

