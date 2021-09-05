#!/usr/bin/env python3

import numpy as np

def traduction(proteine):
	aa_hydrophobe=["A", "I", "L", "M", "F", "W", "V", "Y"]
	proteine_traduite=[]
	for i in proteine:
		if (i in aa_hydrophobe):
			proteine_traduite.append("H")
		else:			
			proteine_traduite.append("P")
	return(proteine_traduite)


def init_matrice(proteine_traduite):
	mat=np.empty((len(proteine_traduite)+4, len(proteine_traduite)+4), dtype= str)
	mat.fill("-")

	for i in range(2, len(mat)):
		for j in range(len(proteine_traduite)):
			mat[2+(len(proteine_traduite)//2), j+2]=proteine_traduite[j]
			
	return(mat)


prot_HP=(traduction("AILMFWVYGA"))
print(prot_HP)
print(init_matrice(prot_HP))


		