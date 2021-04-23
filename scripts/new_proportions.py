prop_surface = [i[4] / i[0] for i in list_taille_nombre_residus]
list_taille_nombre_residus
sum([prop_surface[0], prop_interface1[0], prop_interface2[0]])
list_taille_nombre_residus
from Bio.SeqIO.FastaIO import FastaIterator
help(FastaIterator)
tetramer_names
dSASA
dSASA.keys()
dSASA["notInterfaceRes"]
a = set(dSASA["notInterfaceRes"])
b = set(dSASA["surfaceRes"])
a
b
a.intersection(b)
a == b
len(a), len(b), len(a.intersection(b))
not_interface = a = set(dSASA["notInterfaceRes"])
surfacee = b = set(dSASA["surfaceRes"])
len(not_interface), len(surface), len(not_interface.intersection(surface))
surface = b = set(dSASA["surfaceRes"])
len(not_interface), len(surface), len(not_interface.intersection(surface))
list_taille_nombre_residus[0]
sum(list_taille_nombre_residus[0])
sum(list_taille_nombre_residus[0][1:])
tetramer_names[0]
aber = list_taille_nombre_residus[0][1:].copy()
aber
sum(aber)
aber[0] -= aber[-1]
sum(aber)
aber
list_taille_nombre_residus[0]
tetramer_names
tetramer_names[0]
list_taille_nombre_residus[0]
sum(aber)
sum(list_taille_nombre_residus[0][1:])
aber
list_taille_nombre_residus[0]
list_taille_nombre_residus[0] %/%
list_taille_nombre_residus[0] %/% 4
norm = [ residu / list_taille_nombre_residus[0][0] for residu in list_taille_nombre_residus[0] ]
norm
sum(norm[1:])
prop_interface1
list_taille_nombre_residus
%history -f new_proportions.py
