from dimertetramer.parsing.pdb import PDB_parser, pdb_to_csv
from dimertetramer.utils.customobjs import Path as path, ObjDict as odict

# from dimertetramer.calculations.theo import *
from dimertetramer.calculations import misc as PDBTools
import matplotlib.pyplot as plt

list_tetramere = [
    "1non",
    "1t2a",
    "3mds",
    "1b9b",
    "3pgm",
    "1ub3",
    "1lk5",
    "1m3k",
    "1inl",
    "1f8f",
    "1f8w",
    "1a4s",
    "1a16",
    "1eyi",
    "1rli",
    "1fo6",
    "1m41",
    "1qsm",
    "1bfd",
    "1j2r",
    "1dq8",
]

list_taille_nombre_residus = []
list_taille_SASA = []
list_frequence_types = []

HP = ["ALA", "CYS", "PHE", "ILE", "LEU", "MET", "VAL", "TRP", "TYR"]  # hydrophobe
POL = ["ASN", "GLN", "SER", "THR"]  # polaire
CHG = ["ASP", "GLU", "HIS", "LYS", "ARG"]  # chargé

# Note : la proline et la glycine sont considerees comme des cas a part.

mean_consurf_interface1 = []
mean_consurf_interface2 = []

proj_root = path(".")
assets = proj_root.here("assets").dglob("*")
dimers = assets["Dimers"].dglob("*.pdb")
tetramers = assets["Tetramers"].dglob("*.pdb")
interfaces = assets["Dim_Tet_interfaces"].dglob("*.pdb")

for tetramere in list_tetramere:

    infile = "assets/Dim_Tet_interfaces/{:s}_inter.pdb".format(tetramere)
    print(tetramere)
    dSASA = PDBTools.SASA_parser(infile)

    infile = "assets/Consurf/{:s}_Tet_pdb_With_Conservation_Scores.pdb".format(
        tetramere
    )

    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    for line in lines:
        if line[0:4] == "ATOM" and line[60:67].strip() != "":

            curres = "{:s}".format(line[22:26]).strip()
            dSASA[curres]["consurf"] = float(line[60:67].strip())

    list_taille_nombre_residus.append(
        [
            len(dSASA["reslist"]),
            len(dSASA["notInterfaceRes"]),
            len(dSASA["interfaceRes1"]),
            len(dSASA["interfaceRes2"]),
        ]
    )
    list_taille_SASA.append(
        [dSASA["tSASA"], dSASA["tailleInterface1"], dSASA["tailleInterface2"]]
    )

    HP_surface = 0
    HP_interface1 = 0
    HP_interface2 = 0
    POL_surface = 0
    POL_interface1 = 0
    POL_interface2 = 0
    CHG_surface = 0
    CHG_interface1 = 0
    CHG_interface2 = 0

    for resid in dSASA["surfaceRes"]:
        if dSASA[resid]["resname"] in HP:
            HP_surface += 1
        elif dSASA[resid]["resname"] in POL:
            POL_surface += 1
        elif dSASA[resid]["resname"] in CHG:
            CHG_surface += 1

    for resid in dSASA["interfaceRes1"]:
        if dSASA[resid]["resname"] in HP:
            HP_interface1 += 1
        elif dSASA[resid]["resname"] in POL:
            POL_interface1 += 1
        elif dSASA[resid]["resname"] in CHG:
            CHG_interface1 += 1

    for resid in dSASA["interfaceRes2"]:
        if dSASA[resid]["resname"] in HP:
            HP_interface2 += 1
        elif dSASA[resid]["resname"] in POL:
            POL_interface2 += 1
        elif dSASA[resid]["resname"] in CHG:
            CHG_interface2 += 1

    list_frequence_types.append(
        [
            HP_surface,
            HP_interface1,
            HP_interface2,
            POL_surface,
            POL_interface1,
            POL_interface2,
            CHG_surface,
            CHG_interface1,
            CHG_interface2,
        ]
    )

    list_consurf_interface1 = []
    list_consurf_interface2 = []

    for resid in dSASA["interfaceRes1"]:
        list_consurf_interface1.append(dSASA[resid]["consurf"])

    for resid in dSASA["interfaceRes2"]:
        list_consurf_interface2.append(dSASA[resid]["consurf"])

    mean_consurf_interface1.append(
        sum(list_consurf_interface1) / len(list_consurf_interface1)
    )
    mean_consurf_interface2.append(
        sum(list_consurf_interface2) / len(list_consurf_interface2)
    )
