""" 
    This file contains code not authored by the repository owners.

    File taken from a course example by Anne Lopes.
    https://www.researchgate.net/profile/Anne-Lopes
    Course : Python applied to structural bioinformatics.
"""
import math, string, sys
from ..utils.customobjs import ObjDict


def PDB_parser(infile):
    """version adaptee du parser PDB qui est retourne
    - un dico PDB (meme format que d'habitude
    input: nom du fichier pdb
    output : dPDB (dico PDB)

    """
    # lecture du fichier PDB
    # -----------------------
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    # var ini
    # ---------
    dPDB = ObjDict({})
    atomlist = []
    dPDB["chains"] = []

    # parcoure le PDB
    # -----------------
    for line in lines:
        if line[0:4] == "ATOM":

            # on recupere l'info de la chaine
            chain = line[21]

            # si la chaine n'existe pas, on cree la cle correspondante et on ajoute la chaine a la liste des chaines
            if not chain in dPDB["chains"]:
                dPDB["chains"].append(chain)  # ajout de "chain" a la liste des chaines
                dPDB[chain] = ObjDict({})  # creation du sous-dico pour la chaine
                # on prepare la structure de donnees pour cette chaine
                dPDB[chain]["reslist"] = []

            # on recupere l'info du residu
            curres = "%s" % (line[22:26]).strip()

            # si le residu pour cette chaine "chain" n'existe pas, on cree la cle correspondante et on ajoute le res a la liste des res
            if not curres in dPDB[chain]["reslist"]:
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = ObjDict({})
                # on prepare la structure de donnees pour ce residu
                dPDB[chain][curres]["atomlist"] = []
                # on recupere l'info du residu
                dPDB[chain][curres]["resname"] = line[17:20].strip()

            # on recupere les info pour l'atome de ce res de cette chaine (type atomique + coords x, y, z)
            atomtype = line[12:16].strip()
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = ObjDict({})
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
            dPDB[chain][curres][atomtype]["bfactor"] = line[60:67].strip()

    return dPDB


def writePDB(dPDB, filout="out.pdb", bfactor=False):
    """purpose: according to the coordinates in dPDB, writes the corresponding PDB file.
    If bfactor = True, writes also the information corresponding to the key bfactor
    of each residue (one key per residue) in dPDB.
    input: a dico with the dPDB format
    output: PDB file.
    """

    fout = open(filout, "w")

    for chain in dPDB["chains"]:
        for res in dPDB[chain]["reslist"]:
            for atom in dPDB[chain][res]["atomlist"]:
                if bfactor:
                    fout.write(
                        "ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00%7.3f X X\n"
                        % (
                            dPDB[chain][res][atom]["id"],
                            atom,
                            dPDB[chain][res]["resname"],
                            chain,
                            res,
                            dPDB[chain][res][atom]["x"],
                            dPDB[chain][res][atom]["y"],
                            dPDB[chain][res][atom]["z"],
                            dPDB[chain][res]["bfactor"],
                        )
                    )
                else:
                    fout.write(
                        "ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"
                        % (
                            dPDB[chain][res][atom]["id"],
                            atom,
                            dPDB[chain][res]["resname"],
                            chain,
                            res,
                            dPDB[chain][res][atom]["x"],
                            dPDB[chain][res][atom]["y"],
                            dPDB[chain][res][atom]["z"],
                        )
                    )

    fout.close()
