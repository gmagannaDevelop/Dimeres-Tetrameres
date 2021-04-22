"""
    Utility functions implemented by Theo Roncalli et al.
"""

import math
import numpy as np
import freesasa


def PDB_parser(infile):
    """Adapted draft of PDB parser
    Input: name of PDB file
    Output : dPDB (dico PDB)
    """
    # lecture du fichier PDB
    # -----------------------
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    # var ini
    # ---------
    dPDB = {}
    dPDB["reslist"] = []

    # parcoure le PDB
    # -----------------
    for line in lines:
        if line[0:4] == "ATOM":

            # on recupere l'info du residu
            curres = "{:s}".format(line[22:26]).strip()

            # si le residu pour cette chaine "chain" n'existe pas, on cree la cle correspondante et on ajoute le res a la liste des res
            if not curres in dPDB["reslist"]:
                dPDB["reslist"].append(curres)
                dPDB[curres] = {}
                # on prepare la structure de donnees pour ce residu
                dPDB[curres]["atomlist"] = []
                # on recupere l'info du residu
                dPDB[curres]["resname"] = line[17:20].strip()
                dPDB[curres]["interface"] = int(line[60:67].strip()[0])

            # on recupere les info pour l'atome de ce res de cette chaine (type atomique + coords x, y, z)
            atomtype = line[12:16].strip()
            dPDB[curres]["atomlist"].append(atomtype)
            dPDB[curres][atomtype] = {}
            dPDB[curres][atomtype]["x"] = float(line[30:38])
            dPDB[curres][atomtype]["y"] = float(line[38:46])
            dPDB[curres][atomtype]["z"] = float(line[46:54])
            dPDB[curres][atomtype]["id"] = line[6:11].strip()
    #            dPDB[curres][atomtype]["bfactor"] = line[60:67].strip()

    return dPDB


def SASA_parser(infile):
    """Enhanced draft of PDB parser
    Input: name of PDB file
    Output : dSASA (dico PDB + SASA)
    """

    structure = freesasa.Structure(infile)
    result = freesasa.calc(structure)
    SASA = result.residueAreas()

    dSASA = PDB_parser(infile)

    tSASA = 0

    for resid in dSASA["reslist"]:
        dSASA[resid]["tSASA"] = SASA["A"][resid].total
        dSASA[resid]["rSASA"] = SASA["A"][resid].relativeTotal
        tSASA += SASA["A"][resid].total

    dSASA["tSASA"] = tSASA

    dSASA["notInterfaceRes"] = []
    dSASA["surfaceRes"] = []
    dSASA["interfaceRes1"] = []
    dSASA["interfaceRes2"] = []

    for resid in dSASA["reslist"]:
        if dSASA[resid]["interface"] == 0:
            dSASA["notInterfaceRes"].append(resid)
            if dSASA[resid]["rSASA"] >= 0.25:
                dSASA["surfaceRes"].append(resid)
        elif dSASA[resid]["interface"] == 1:
            dSASA["interfaceRes1"].append(resid)
        elif dSASA[resid]["interface"] == 2:
            dSASA["interfaceRes2"].append(resid)

    dSASA["tailleInterface1"] = 0.0
    dSASA["tailleInterface2"] = 0.0

    for resid in dSASA["interfaceRes1"]:
        dSASA["tailleInterface1"] += dSASA[resid]["tSASA"]
    for resid in dSASA["interfaceRes2"]:
        dSASA["tailleInterface2"] += dSASA[resid]["tSASA"]

    dSASA["seqSurface"] = []
    dSASA["seqInterface1"] = []
    dSASA["seqInterface2"] = []

    for resid in dSASA["surfaceRes"]:
        dSASA["seqSurface"].append(dSASA[resid]["resname"])
    for resid in dSASA["interfaceRes1"]:
        dSASA["seqInterface1"].append(dSASA[resid]["resname"])
    for resid in dSASA["interfaceRes2"]:
        dSASA["seqInterface2"].append(dSASA[resid]["resname"])

    AA = [
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "GLY",
        "HIS",
        "ILE",
        "LYS",
        "LEU",
        "MET",
        "ASN",
        "PRO",
        "GLN",
        "ARG",
        "SER",
        "THR",
        "VAL",
        "TRP",
        "TYR",
    ]

    dSASA["freqSurface"] = []
    dSASA["freqInterface1"] = []
    dSASA["freqInterface2"] = []

    for resid in AA:
        dSASA["freqSurface"].append(dSASA["seqSurface"].count(resid))
        dSASA["freqInterface1"].append(dSASA["seqInterface1"].count(resid))
        dSASA["freqInterface2"].append(dSASA["seqInterface2"].count(resid))

    return dSASA


def countAA(dPDB, resname):
    """Count the number of an amino acid in each chain.
    Input: name of PDB file, name of amino acid in 3 letter code
    Output : a dico containing the number of the amino acid for each chain.
    """

    dAAPerChain = {}

    for chain in dPDB["chains"]:
        n = 0

        for aa in dPDB["reslist"]:
            if dPDB[aa]["resname"] == resname:
                n += 1

        dAAPerChain = n

    return dAAPerChain


def writePDB(dPDB, filout="out.pdb", bfactor=False):
    """According to the coordinates in dPDB, writes the corresponding PDB file.
    If bfactor = True, writes also the information corresponding to the key bfactor
    of each residue (one key per residue) in dPDB.
    input: a dico with the dPDB format
    output: PDB file.
    """

    fout = open(filout, "w")

    for chain in dPDB["chains"]:
        for res in dPDB["reslist"]:
            for atom in dPDB[res]["atomlist"]:
                if bfactor:
                    fout.write(
                        "ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00%7.3f X X\n"
                        % (
                            dPDB[res][atom]["id"],
                            atom,
                            dPDB[res]["resname"],
                            chain,
                            res,
                            dPDB[res][atom]["x"],
                            dPDB[res][atom]["y"],
                            dPDB[res][atom]["z"],
                            dPDB[res]["bfactor"],
                        )
                    )
                else:
                    fout.write(
                        "ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"
                        % (
                            dPDB[res][atom]["id"],
                            atom,
                            dPDB[res]["resname"],
                            chain,
                            res,
                            dPDB[res][atom]["x"],
                            dPDB[res][atom]["y"],
                            dPDB[res][atom]["z"],
                        )
                    )

    fout.close()


def ComputeDistAtom(dAtom1, dAtom2):
    """Compute the distance between two atoms.
    Input: dAtom1, dAtom2 which are dico corresponding to atom 1 and atom 2 respectively.
           For instance dAtom = dPDB[resnumber][atom].
    Output: distance (float)
    """

    coord1 = [dAtom1["x"], dAtom1["y"], dAtom1["z"]]
    coord2 = [dAtom2["x"], dAtom2["y"], dAtom2["z"]]

    dist = math.sqrt(
        (coord1[0] - coord2[0]) ** 2
        + (coord1[1] - coord2[1]) ** 2
        + (coord1[2] - coord2[2]) ** 2
    )

    return float(dist)


def ComputeBarycenter(dRes):
    """Compute the barycenter of an amino acid.
    Input: a dico containing information about a residual.
           For instance dBPDB[resnumber].
    Output: coordinates (x,y,z)
    """

    n = len(dRes["atomlist"])
    x = []
    y = []
    z = []

    for atom in dRes["atomlist"]:
        x.append(float(dRes[atom]["x"]))
        y.append(float(dRes[atom]["y"]))
        z.append(float(dRes[atom]["z"]))

    xbar = sum(x) / n
    ybar = sum(y) / n
    zbar = sum(z) / n

    return (xbar, ybar, zbar)


def ComputeDistance(dRes1, dRes2, method="center"):
    """Compute the distance between two amino acids.
    Input: two dico containing information about a residual and method used for computing distance.
           For instance : dRes = dPDB[resnumber].
    method 'center': computes the distance between the two centers of mass
                     of the two residues and returns it (default = 'center')
    method 'atom': computes the distance between all the atoms of the two
                   residues and returns the smallest distance.
    Output : distance between two residues.
    """

    dist = float("Inf")

    if method == "center":

        coord1 = ComputeBarycenter(dRes1)
        coord2 = ComputeBarycenter(dRes2)

        dist = math.sqrt(
            (coord1[0] - coord2[0]) ** 2
            + (coord1[1] - coord2[1]) ** 2
            + (coord1[2] - coord2[2]) ** 2
        )

    elif method == "atom":

        distances = []

        for atom1 in dRes1["atomlist"]:

            coord1 = [dRes1[atom1]["x"], dRes1[atom1]["y"], dRes1[atom1]["z"]]

            for atom2 in dRes2["atomlist"]:

                coord2 = [dRes2[atom2]["x"], dRes2[atom2]["y"], dRes2[atom2]["z"]]

                distances.append(
                    math.sqrt(
                        (coord1[0] - coord2[0]) ** 2
                        + (coord1[1] - coord2[1]) ** 2
                        + (coord1[2] - coord2[2]) ** 2
                    )
                )

        dist = min(distances)

    return dist


def ComputeDistanceMatrix(dPDB, chain1, chain2, method="center"):
    """Compute the distance matrix, i.e. the distance between each residu of two chains (or of a unique chain).
    Input: name of PDB file, name of chain1 and chain2 and method used for computing distance.
    method 'center': computes the distance between the two centers of mass
                     of the two residues and returns it (default = 'center')
    method 'atom': computes the distance between all the atoms of the two
                   residues and returns the smallest distance.
    Output : distance-based matrix
    """

    max_i = max(list(map(int, dPDB[chain1]["reslist"])))
    max_j = max(list(map(int, dPDB[chain2]["reslist"])))

    distMat = np.empty((max_i + 1, max_j + 1))
    distMat[:] = np.nan

    for i in dPDB[chain1]["reslist"]:
        for j in dPDB[chain2]["reslist"]:
            distMat[int(i), int(j)] = ComputeDistance(
                dPDB[chain1][i], dPDB[chain2][j], method
            )

    return distMat


def ExtractContactResidues(distMat, threshold=5):
    """returns pairs of residues in contacts (threshold).
    Input: distance matrix and threshold (default = 5)
    Output : list of lists which contain the residues in contacts.
    """

    contacts = []

    for i in range(distMat.shape[0]):
        for j in range(i + 1, distMat.shape[1]):
            if distMat[i, j] <= threshold:
                contacts.append([i, j])

    return contacts