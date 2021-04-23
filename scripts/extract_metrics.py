import pandas as pd

from functools import reduce
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

proj_root = path("..")
assets = proj_root.here("assets").dglob("*")
dimers = assets["Dimers"].dglob("*.pdb")
tetramers = assets["Tetramers"].dglob("*.pdb")
interfaces = assets["Dim_Tet_interfaces"].dglob("*.pdb")

aa_counts_interface1 = []
aa_counts_interface2 = []
aa_counts_surface = []

for iii, tetramere in enumerate(tetramer_names):

    tetramer_file = f"../assets/Dim_Tet_interfaces/{tetramere}_inter.pdb"
    print(iii, tetramere)
    dSASA = PDBTools.SASA_parser(tetramer_file)

    consurf_file = f"../assets/Consurf/{tetramere}_Tet_pdb_With_Conservation_Scores.pdb"

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

    _aa_surface = [0] * len(aa.all)
    for residue in dSASA["surfaceRes"]:
        _resname = dSASA[residue]["resname"]
        _j = aa.all.index(_resname)
        _aa_surface[_j] += 1

        if _resname in aa.hydrophobic:
            HP_surface += 1
        elif _resname in aa.polar:
            POL_surface += 1
        elif _resname in aa.charged:
            CHG_surface += 1
    aa_counts_surface.append(_aa_surface)

    _aa_interface1 = [0] * len(aa.all)
    for residue in dSASA["interfaceRes1"]:
        _resname = dSASA[residue]["resname"]
        _j = aa.all.index(_resname)
        _aa_interface1[_j] += 1

        if _resname in aa.hydrophobic:
            HP_interface1 += 1
        elif _resname in aa.polar:
            POL_interface1 += 1
        elif _resname in aa.charged:
            CHG_interface1 += 1
    aa_counts_interface1.append(_aa_interface1)

    _aa_interface2 = [0] * len(aa.all)
    for residue in dSASA["interfaceRes2"]:
        _resname = dSASA[residue]["resname"]
        _j = aa.all.index(_resname)
        _aa_interface2[_j] += 1

        if _resname in aa.hydrophobic:
            HP_interface2 += 1
        elif _resname in aa.polar:
            POL_interface2 += 1
        elif _resname in aa.charged:
            CHG_interface2 += 1
    aa_counts_interface2.append(_aa_interface2)

    list_frequence_types.append(
        [
            HP_surface / len(dSASA["surfaceRes"]),
            HP_interface1 / len(dSASA["interfaceRes1"]),
            HP_interface2 / len(dSASA["interfaceRes2"]),
            POL_surface / len(dSASA["surfaceRes"]),
            POL_interface1 / len(dSASA["interfaceRes1"]),
            POL_interface2 / len(dSASA["interfaceRes2"]),
            CHG_surface / len(dSASA["surfaceRes"]),
            CHG_interface1 / len(dSASA["interfaceRes1"]),
            CHG_interface2 / len(dSASA["interfaceRes2"]),
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

surface_HP = HP_surface = [i[0] for i in list_frequence_types]
interface1_HP = HP_interface1 = [i[1] for i in list_frequence_types]
interface2_HP = HP_interface2 = [i[2] for i in list_frequence_types]
surface_POL = POL_surface = [i[3] for i in list_frequence_types]
interface1_POL = POL_interface1 = [i[4] for i in list_frequence_types]
interface2_POL = POL_interface2 = [i[5] for i in list_frequence_types]
surface_CHG = CHG_surface = [i[6] for i in list_frequence_types]
interface1_CHG = CHG_interface1 = [i[7] for i in list_frequence_types]
interface2_CHG = CHG_interface2 = [i[8] for i in list_frequence_types]

prop_not_interface = [i[1] / i[0] for i in list_taille_nombre_residus]
prop_interface1 = [i[2] / i[0] for i in list_taille_nombre_residus]
prop_interface2 = [i[3] / i[0] for i in list_taille_nombre_residus]

SASA_total = [sasas[0] for sasas in list_taille_SASA]
SASA_interface1 = [sasas[1] for sasas in list_taille_SASA]
SASA_interface2 = [sasas[2] for sasas in list_taille_SASA]

lists = {
    i: eval(i) for i in dir() if isinstance(eval(i), list) and not i.startswith("_")
}

flatten_nested_aa_lists = lambda aa_list, name: {
    f"{name} count {aa.all[i]}": [aa_list[j][i] for j in range(len(aa_list))]
    for i in range(len(aa.all))
}

aa_surface = flatten_nested_aa_lists(aa_counts_surface, "surface")
aa_interface1 = flatten_nested_aa_lists(aa_counts_interface1, "interface1")
aa_interface2 = flatten_nested_aa_lists(aa_counts_interface2, "interface2")

aa_occurrences = {}
aa_occurrences.update(aa_surface)
aa_occurrences.update(aa_interface1)
aa_occurrences.update(aa_interface2)

metrics_dict = {
    "SASA total": SASA_total,
    "SASA interface1": SASA_interface1,
    "SASA interface2": SASA_interface2,
    "consurf interface1": mean_consurf_interface1,
    "consurf interface2": mean_consurf_interface2,
    "relative size not interface": prop_not_interface,
    "relative size interface1": prop_interface1,
    "relative size interface2": prop_interface2,
    "surface hydrophobic": surface_HP,
    "surface polar": surface_POL,
    "surface charged": surface_CHG,
    "interface1 hydrophobic": interface1_HP,
    "interface1 polar": interface1_POL,
    "interface1 charged": interface1_CHG,
    "interface2 hydrophobic": interface2_HP,
    "interface2 polar": interface2_POL,
    "interface2 charged": interface2_CHG,
}

metrics_dict.update(aa_occurrences)

metrics_df: pd.DataFrame = pd.DataFrame(metrics_dict, index=tetramer_names)
metrics_df.index.name = "protein"
