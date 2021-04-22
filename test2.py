import matplotlib.pyplot as plt

from dimertetramer.parsing.pdb import PDB_parser, pdb_to_csv
from dimertetramer.utils.customobjs import Path as path, ObjDict as odict
from dimertetramer.calculations import misc as PDBTools
from dimertetramer.calculations.constants import AminoAcids as aa

tetramer_names = [
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

# Note : la proline et la glycine sont considerees comme des cas a part.

mean_consurf_interface1 = []
mean_consurf_interface2 = []

proj_root = path(".")
assets = proj_root.here("assets").dglob("*")
dimers = assets["Dimers"].dglob("*.pdb")
tetramers = assets["Tetramers"].dglob("*.pdb")
interfaces = assets["Dim_Tet_interfaces"].dglob("*.pdb")

for tetramere in tetramer_names:

    tetramer_file = f"assets/Dim_Tet_interfaces/{tetramere}_inter.pdb"
    print(tetramere)
    dSASA = PDBTools.SASA_parser(tetramer_file)

    consurf_file = f"assets/Consurf/{tetramere}_Tet_pdb_With_Conservation_Scores.pdb"

    with open(consurf_file, "r") as f:
        lines = f.readlines()

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

    for residue in dSASA["surfaceRes"]:
        if dSASA[residue]["resname"] in aa.hydrophobic:
            HP_surface += 1
        elif dSASA[residue]["resname"] in aa.polar:
            POL_surface += 1
        elif dSASA[residue]["resname"] in aa.charged:
            CHG_surface += 1

    for residue in dSASA["interfaceRes1"]:
        if dSASA[residue]["resname"] in aa.hydrophobic:
            HP_interface1 += 1
        elif dSASA[residue]["resname"] in aa.polar:
            POL_interface1 += 1
        elif dSASA[residue]["resname"] in aa.charged:
            CHG_interface1 += 1

    for residue in dSASA["interfaceRes2"]:
        if dSASA[residue]["resname"] in aa.hydrophobic:
            HP_interface2 += 1
        elif dSASA[residue]["resname"] in aa.polar:
            POL_interface2 += 1
        elif dSASA[residue]["resname"] in aa.charged:
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

    for residue in dSASA["interfaceRes1"]:
        list_consurf_interface1.append(dSASA[residue]["consurf"])

    for residue in dSASA["interfaceRes2"]:
        list_consurf_interface2.append(dSASA[residue]["consurf"])

    mean_consurf_interface1.append(
        sum(list_consurf_interface1) / len(list_consurf_interface1)
    )
    mean_consurf_interface2.append(
        sum(list_consurf_interface2) / len(list_consurf_interface2)
    )
