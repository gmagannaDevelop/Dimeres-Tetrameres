import pandas as pd

from dimertetramer.parsing.pdb import PDB_parser, pdb_to_csv
from dimertetramer.utils.customobjs import Path as path, ObjDict as odict
from dimertetramer.calculations.misc import *
from dimertetramer.calculations.sasa import sasa_from_file

proj_root = path(".")
assets = proj_root.here("assets").dglob("*")

dimers = assets["Dimers"].dglob("*.pdb")
tetramers = assets["Tetramers"].dglob("*.pdb")
interfaces = assets["Dim_Tet_interfaces"].dglob("*.pdb")


sasas = odict(
    {
        "dimers": odict(
            {
                molecule.replace(".pdb", ""): sasa_from_file(file)
                for molecule, file in dimers.items()
            }
        ),
        "tetramers": odict(
            {
                molecule.replace(".pdb", ""): sasa_from_file(file)
                for molecule, file in tetramers.items()
            }
        ),
        "interfaces": odict(
            {
                molecule.replace(".pdb", ""): sasa_from_file(file)
                for molecule, file in interfaces.items()
            }
        ),
    }
)
