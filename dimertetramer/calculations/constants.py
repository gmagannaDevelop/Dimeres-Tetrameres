from ..utils.customobjs import ObjDict

AminoAcids = ObjDict(
    {
        "all": [
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
        ],
        "hydrophobic": ["ALA", "CYS", "PHE", "ILE", "LEU", "MET", "VAL", "TRP", "TYR"],
        "polar": ["ASN", "GLN", "SER", "THR"],
        "charged": ["ASP", "GLU", "HIS", "LYS", "ARG"],
    }
)
