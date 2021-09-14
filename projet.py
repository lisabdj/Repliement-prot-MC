#!/usr/bin/env python3

import numpy as np
import argparse
import random
import math
import sys
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Repliement de protéine grâce à un algorithme de Monte Carlo')
parser.add_argument('proteine', metavar="sequence proteique", type=str, nargs="+",
                    help="une séquence protéique à replier")
parser.add_argument('temperature', metavar="T", type=float, nargs="+",
                    help="temperature")

args = parser.parse_args()
print(args)

# Variables globales
SEQ_PROT = sys.argv[1]
TEMP = sys.argv[2]

if len(sys.argv) != 3:
    sys.exit("ERREUR : N'oubliez pas de renseigner la séquence protéique et la température du repliement.")


def traduction(proteine):
    """Traduction de la séquence protéique en acides aminés en séquence HP"""
    aa_hydrophobe = ["A", "I", "L", "M", "F", "W", "V", "Y"]
    proteine_traduite = []
    for i in proteine:
        if i in aa_hydrophobe:
            proteine_traduite.append("H")
        else:
            proteine_traduite.append("P")
    return proteine_traduite


def init_matrice(proteine_traduite):
    """Initialisation de la matrice bidimensionnelle"""
    mat = np.empty((len(proteine_traduite) + 4, len(proteine_traduite) + 4), dtype=str)
    mat.fill("-")

    for i in range(2, len(mat)):
        for j in range(len(proteine_traduite)):
            mat[2 + (len(proteine_traduite) // 2), j + 2] = proteine_traduite[j]

    coordonnees_mat = coordonnees_init(mat)
    return mat, coordonnees_mat


def coordonnees_init(matrice):
    """Calcul des coordonnées des acides aminés sur la matrice bidimensionnelle"""
    coordonnees = []
    for i in range(len(matrice)):
        for j in range(len(matrice)):
            if matrice[i, j] != "-":
                coordonnees.append([i, j])
    return coordonnees


def vshd_moves(matrice, coordonnees, aa_position, seq_HP):
    """VSHD mouvement"""

    # Acide aminé situé aux extrémités
    if aa_position == 0 or aa_position == len(coordonnees) - 1:
        mvt_disponible = vshd_moves_end_moves(matrice, coordonnees, aa_position)
        if mvt_disponible != None:
            if len(mvt_disponible) == 1:
                coordonnees[aa_position] = mvt_disponible[0]
            else:
                i = random.randint(0, 1)
                coordonnees[aa_position] = mvt_disponible[i]

            matrice = matrice_update(coordonnees, seq_HP)
            return matrice, coordonnees

    # Acide aminé non situé aux extrémités
    if aa_position != 0 & aa_position != len(coordonnees) - 1:

        forme = forme_u(coordonnees, aa_position)
        if forme[0] == 1:
            crankshaft = vshd_moves_crankshaft_moves(matrice, coordonnees, aa_position)
            if coordonnees[2] == 0:
                coordonnees[aa_position] = crankshaft[0]
                coordonnees[aa_position + 1] = crankshaft[1]
            elif coordonnees[2] == 1:
                coordonnees[aa_position] = crankshaft[0]
                coordonnees[aa_position - 1] = crankshaft[1]
            matrice = matrice_update(coordonnees, seq_HP)
            return matrice, coordonnees

        elif forme[0] == 0:
            mvt_possible = vshd_moves_corner_moves(matrice, coordonnees, aa_position)
            if mvt_possible != []:
                coordonnees[aa_position] = mvt_possible[0]
                matrice = matrice_update(coordonnees, seq_HP)
                return matrice, coordonnees
    return (0)


def vshd_moves_end_moves(matrice, coordonnees, aa_position):
    """End moves du VSHD moves sur le premier ou le dernier résidu"""
    cases_vides = recherche_cases_vides(matrice, coordonnees, aa_position)
    cases_possibles = []  # stocker les coordonnées de places possibles

    for vide in cases_vides:
        if not (vide in coordonnees):
            cases_possibles.append(vide)

    return cases_possibles


def vshd_moves_corner_moves(matrice, coordonnees, aa_position):
    """Corner moves du VSHD moves, sur tous les résidus sauf les extrémités"""
    cases_vides = recherche_cases_vides(matrice, coordonnees, aa_position)
    cases_possibles = []

    for vide in cases_vides:
        if cases_vides.count(vide) != 1:
            cases_possibles.append(vide)
            return cases_possibles


def vshd_moves_crankshaft_moves(matrice, coordonnees, aa_position):
    """Mouvement du vilebrequin du VSHD, sur tous les res sauf extrémités, à condition d'avoir dans une forme U"""
    cases_possibles = []
    forme, aa_deplacement = forme_u(coordonnees, aa_position)
    cases_vides_aa1 = recherche_cases_vides(matrice, coordonnees, aa_deplacement[0])  # position de l'acide aminé
    cases_vides_aa2 = recherche_cases_vides(matrice, coordonnees, aa_deplacement[1])  # position de l'acide aminé voisin

    if aa_deplacement[2] == "x":
        if aa_deplacement[4] == "+1":
            c1 = [coordonnees[aa_position][0] + 2, coordonnees[aa_position][1]]
            c2 = [coordonnees[aa_position][0] + 2, coordonnees[aa_deplacement[1]][1]]
            if c1 in cases_vides_aa1 and c2 in cases_vides_aa2:
                cases_possibles.append(c1)
                cases_possibles.append(c2)
        elif aa_deplacement[4] == "-1":
            c1 = [coordonnees[aa_position][0] - 2, coordonnees[aa_position][1]]
            c2 = [coordonnees[aa_position][0] - 2, coordonnees[aa_deplacement[1]][1]]
            if c1 in cases_vides_aa1 and c2 in cases_vides_aa2:
                cases_possibles.append(c1)
                cases_possibles.append(c2)
    elif aa_deplacement[2] == "y":
        if aa_deplacement[4] == "+1":
            c1 = [coordonnees[aa_position][0], coordonnees[aa_position][1] + 2]
            c2 = [coordonnees[aa_deplacement[1]][0], coordonnees[aa_position][1] + 2]
            if c1 in cases_vides_aa1 and c2 in cases_vides_aa2:
                cases_possibles.append(c1)
                cases_possibles.append(c2)
        elif aa_deplacement[4] == "-1":
            c1 = [coordonnees[aa_position][0], coordonnees[aa_position][1] - 2]
            c2 = [coordonnees[aa_deplacement[1]][0], coordonnees[aa_position][1] - 2]
            if c1 in cases_vides_aa1 and c2 in cases_vides_aa2:
                cases_possibles.append(c1)
                cases_possibles.append(c2)

    if aa_deplacement[0] - aa_deplacement[1] < 0:
        ordre = 0  # aa_position est le premier acide aminé
    else:
        ordre = 1

    return c1, c2, ordre


def forme_u(coordonnees, aa_position):
    """Déterminer si l'acide aminé est dans une forme U"""

    x = coordonnees[aa_position][0]
    y = coordonnees[aa_position][1]

    forme_u = 0
    aa_u = []
    if coordonnees[aa_position + 1][0] == x and aa_position >= 2 and len(coordonnees) >= aa_position + 4:
        if coordonnees[aa_position + 2][0] == x + 1 or coordonnees[aa_position + 2][0] == x - 1:
            ref = coordonnees[aa_position + 2][0]
            if (coordonnees[aa_position + 3][0] == ref and coordonnees[aa_position - 1][0] == ref and
                    coordonnees[aa_position - 2][0] == ref):
                forme_u = 1
                aa_u.append(aa_position)
                aa_u.append(aa_position + 1)
                aa_u.append("x")
                aa_u.append(0)
                if coordonnees[aa_position + 2][0] == x + 1:
                    aa_u.append("+1")
                elif coordonnees[aa_position + 2][0] == x - 1:
                    aa_u.append("-1")
    elif coordonnees[aa_position + 1][1] == y and aa_position >= 2 and len(coordonnees) >= aa_position + 4:
        if coordonnees[aa_position + 2][1] == y + 1 or coordonnees[aa_position + 2][1] == y - 1:
            ref = coordonnees[aa_position + 2][1]
            if (coordonnees[aa_position + 3][1] == ref and coordonnees[aa_position - 1][1] == ref and
                    coordonnees[aa_position - 2][1] == ref):
                forme_u = 1
                aa_u.append(aa_position)
                aa_u.append(aa_position + 1)
                aa_u.append("y")
                aa_u.append(0)
                if coordonnees[aa_position + 2][1] == y + 1:
                    aa_u.append("+1")
                elif coordonnees[aa_position + 2][1] == y - 1:
                    aa_u.append("-1")
    elif coordonnees[aa_position - 1][0] == x and aa_position >= 3 and len(coordonnees) >= aa_position + 3:
        if coordonnees[aa_position - 2][0] == x + 1 or coordonnees[aa_position - 2][0] == x - 1:
            ref = coordonnees[aa_position - 2][0]
            if (coordonnees[aa_position - 3][0] == ref and coordonnees[aa_position + 1][0] == ref and
                    coordonnees[aa_position + 2][0] == ref):
                forme_u = 1
                aa_u.append(aa_position - 1)
                aa_u.append(aa_position)
                aa_u.append("x")
                aa_u.append(1)
                if coordonnees[aa_position + 2][0] == x + 1:
                    aa_u.append("+1")
                elif coordonnees[aa_position + 2][0] == x - 1:
                    aa_u.append("-1")
    elif coordonnees[aa_position - 1][1] == y and aa_position >= 3 and len(coordonnees) >= aa_position + 3:
        if coordonnees[aa_position - 2][1] == y + 1 or coordonnees[aa_position - 2][1] == y - 1:
            ref = coordonnees[aa_position - 2][1]
            if (coordonnees[aa_position - 3][1] == ref and coordonnees[aa_position + 1][1] == ref and
                    coordonnees[aa_position + 2][1] == ref):
                forme_u = 1
                aa_u.append(aa_position - 1)
                aa_u.append(aa_position)
                aa_u.append("y")
                aa_u.append(1)
                if coordonnees[aa_position + 2][1] == y + 1:
                    aa_u.append("+1")
                elif coordonnees[aa_position + 2][1] == y - 1:
                    aa_u.append("-1")
    else:
        forme_u = 0

    return forme_u, aa_u


def deplacement_possible(matrice, coordonnees, seq_HP):
    acide_amine_deplacable = []
    for i in range(len(coordonnees)):
        retour_fonction = None
        if i == 0 or i == len(coordonnees) - 1:
            retour_fonction = vshd_moves_end_moves(matrice, coordonnees, i)
        else:
            forme = forme_u(coordonnees, i)
            if forme[0] == 1:
                retour_fonction = vshd_moves_crankshaft_moves(matrice, coordonnees, i)
            elif forme[0] == 0:
                retour_fonction = vshd_moves_corner_moves(matrice, coordonnees, i)

        if retour_fonction != None:
            acide_amine_deplacable.append(i)

    return acide_amine_deplacable


def pull_move(matrice, coordonnees, aa_position, seq_HP):  # Fonction non finie
    """Effectue le mouvement de traction"""
    '''diago=recherche_diagonale_vide(matrice, coordonnees, aa_position)
    l=recherche_caseL(matrice, coordonnees, aa_position)
    c=recherche_caseC(matrice, coordonnees, aa_position, diago)
    
    if(c[0]==None):
        coordonnees[aa_position]=c[1]
        matrice=matrice_update(coordonnees, seq_HP)
        return(matrice, coordonnees)
    else:'''


def recherche_diagonale_vide(matrice, coordonnees, aa_position):  # Fonction non utilisée dans le programme
    """Recherche les cases vides en diagonale"""
    diago_vide = []
    x = coordonnees[aa_position][0]
    y = coordonnees[aa_position][1]

    if matrice[x - 1][y - 1] == "-":
        diago_vide.append([x - 1, y - 1])
    if matrice[x - 1][y + 1] == "-":
        diago_vide.append([x - 1, y + 1])
    if matrice[x + 1][y - 1] == "-":
        diago_vide.append([x + 1, y - 1])
    if matrice[x + 1][y + 1] == "-":
        diago_vide.append([x + 1, y + 1])

    return diago_vide


def recherche_caseC(matrice, coordonnees, aa_position, diago):  # Fonction non utilisée dans le programme
    """Trouver la case C adjacente à l'acide aminé et à la diagonale de l'acide aminé"""

    cases_L = recherche_caseL(matrice, coordonnees, aa_position)
    xi = coordonnees[aa_position][0]
    yi = coordonnees[aa_position][1]
    cases_C = []

    for i in range(len(cases_L)):
        xL = cases_L[i][0]
        yL = cases_L[i][1]
        if aa_position == 0 or aa_position == len(coordonnees) - 1:
            if matrice[xi][yL] == "-":
                cases_C.append([xi, yL])
            elif matrice[xL][yi] == "-":
                cases_C.append([xL, yi])
        elif aa_position != 0 and aa_position != len(coordonnees) - 1:
            if coordonnees[aa_position + 1] == [xi, yL] and coordonnees[aa_position - 1] == [xL, yi]:
                cases_C = None
                L = [xL, yL]
                return (cases_C, L)
            elif coordonnees[aa_position + 1] == [xL, yi] and coordonnees[aa_position - 1] == [xi, yL]:
                cases_C = None
                L = [xL, yL]
                return cases_C, L

        if cases_C == None:
            cases_L = L

    return cases_C, cases_L


def recherche_caseL(matrice, coordonnees, aa_position):  # Fonction non utilisée dans le programme
    """Recherche la case L adjacente à la diagonale de l'acide aminé et au voisin de l'acide aminé"""
    cases_L = []
    cases_diagonales = recherche_diagonale_vide(matrice, coordonnees, aa_position)
    cases_vides = recherche_cases_vides(matrice, coordonnees, aa_position)

    for diago in cases_diagonales:
        if diago in cases_vides:
            cases_L.append(diago)

    return cases_L


def recherche_cases_vides(matrice, coordonnees, aa_position):
    """Recherche des cases vides autour du/des voisins de l'acide aminé"""

    cases_vides = []

    # Stocker les coordonnées du/des voisins de l'acide aminé
    if aa_position == 0:
        x = coordonnees[aa_position + 1][0]
        y = coordonnees[aa_position + 1][1]
    elif aa_position == len(coordonnees) - 1:
        x = coordonnees[aa_position - 1][0]
        y = coordonnees[aa_position - 1][1]
    else:
        x = coordonnees[aa_position - 1][0], coordonnees[aa_position + 1][0]
        y = coordonnees[aa_position - 1][1], coordonnees[aa_position + 1][1]

    if type(x) == tuple:
        for xi, yi in zip(x, y):
            if matrice[xi - 1][yi] == "-":
                cases_vides.append([xi - 1, yi])
            if matrice[xi][yi + 1] == "-":
                cases_vides.append([xi, yi + 1])
            if matrice[xi + 1][yi] == "-":
                cases_vides.append([xi + 1, yi])
            if matrice[xi][yi - 1] == "-":
                cases_vides.append([xi, yi - 1])
    if type(x) == int:
        if matrice[x - 1][y] == "-":
            cases_vides.append([x - 1, y])
        if matrice[x][y + 1] == "-":
            cases_vides.append([x, y + 1])
        if matrice[x + 1][y] == "-":
            cases_vides.append([x + 1, y])
        if matrice[x][y - 1] == "-":
            cases_vides.append([x, y - 1])

    return cases_vides


def matrice_update(coordonnees, seq_HP):
    """Mise à jour de la matrice"""
    matrice = np.empty((len(coordonnees) + 4, len(coordonnees) + 4), dtype=str)
    matrice.fill("-")

    for i in range(len(coordonnees)):
        x = coordonnees[i][0]
        y = coordonnees[i][1]
        matrice[x, y] = seq_HP[i]
    return matrice


def premier_mvt(matrice, coordonnees, aa_position, seq_HP):
    """Réalise le premier mouvement sur l'une des deux extrémités"""
    mvt_disponible = recherche_cases_vides(matrice, coordonnees, aa_position)
    i = random.randint(0, len(mvt_disponible) - 1)
    coordonnees[aa_position] = mvt_disponible[i]
    matrice = matrice_update(coordonnees, seq_HP)
    return matrice, coordonnees


def energie_hydrophobe(coordonnees, seq_HP):
    """Calcul de l'énergie"""
    pos_H = []
    energie = 0
    for i in range(len(seq_HP)):
        if seq_HP[i] == "H":
            pos_H.append(i)
    for h in range(len(pos_H)):
        for h1 in range(h, len(pos_H)):
            if not (h1 == h + 1 or h1 == h):
                x1 = coordonnees[h][0]
                y1 = coordonnees[h][1]

                x2 = coordonnees[h1][0]
                y2 = coordonnees[h1][1]

                if x1 == x2 or y1 == y2:
                    energie += 1
    return -energie


def procedure_MC(n, matrice, coordonnees, seq_HP):
    """Procédure Monte Carlo : Mouvement VSHD et calcul d'énergie"""
    matrice_new = matrice
    coordonnees_new = coordonnees
    energie_procedure = []
    liste_matrice = []
    liste_coordonnees = []
    for i in range(n):
        aa = random.choice(deplacement_possible(matrice_new, coordonnees_new, prot_HP))
        matrice_new, coordonnees_new = vshd_moves(matrice_new, coordonnees_new, aa, prot_HP)
        diffE = energie_hydrophobe(coordonnees_new, seq_HP) - energie_hydrophobe(coordonnees, seq_HP)
        if diffE <= 0:
            matrice = matrice_new
            coordonnees = coordonnees
            print(matrice)
        else:
            q = random.randint(0, 1)
            if q > math.exp(-diffE / TEMP):
                matrice = matrice_new
                coordonnees = coordonnees_new
                print(matrice)

        liste_matrice.append(matrice)
        liste_coordonnees.append(coordonnees)
        energie_procedure.append(energie_hydrophobe(coordonnees, seq_HP))

    return liste_matrice, liste_coordonnees, energie_procedure


# Main

## Traduction de la protéine selon le modèle HP
prot_HP = (traduction(SEQ_PROT))

## Initialisation : conformation C0
matrice, coordonnees = init_matrice(prot_HP)
print(matrice)

## Conformation C1 sur l'une des deux extrémités
N = len(coordonnees) - 1
aa = random.randint(0, 1)
if aa == 0:
    matrice, coordonnees = premier_mvt(matrice, coordonnees, aa, prot_HP)
elif aa == 1:
    matrice, coordonnees = premier_mvt(matrice, coordonnees, N, prot_HP)
print(matrice)

## Conformation C2 à Cn
matrice, coordonnees, energie = procedure_MC(100, matrice, coordonnees, prot_HP)

## Graphique énergie
# print(energie)
plt.title("Energie des conformations")
plt.plot([i for i in range(len(energie))], energie)
plt.xlabel("Temps")
plt.ylabel("Energie")
plt.show()
