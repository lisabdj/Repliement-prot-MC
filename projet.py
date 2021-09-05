#!/usr/bin/env python3

import numpy as np

def traduction(proteine):
	aa_hydrophobe=["A", "I", "L", "M", "F", "W", "V", "Y"]
	proteine_traduite=[]

	for i in proteine:
		#&print(i)
		if (i in aa_hydrophobe):
				
			proteine_traduite.append("H")
				#h=0
				#print(i, "est polaire")
		else:			
			proteine_traduite.append("P")
				#print(i, "est hydrophobe")
	#chaine_traduite="".join(proteine_traduite)
	return(proteine_traduite)


def matrice_bidimensionnelle(proteine_traduite):
	mat=np.zeros((len(proteine_traduite)+4,len(proteine_traduite)+4))
	for i in range(2, len(mat)):
		for aa in proteine_traduite:
			mat[i, 2+len(proteine_traduite)/2]=aa
			
	return(mat)



prot_HP=(traduction("AILMFWVYGA"))
print(prot_HP)
print(matrice_bidimensionnelle(prot_HP))


		