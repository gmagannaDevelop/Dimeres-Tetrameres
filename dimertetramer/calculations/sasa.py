import pathlib

import Bio.PDB as pdb
import freesasa


from typing import (
    Dict,
    List,
    Tuple,
    NoReturn,
    Union,
    Optional,
    Any,
    Callable,
    Iterable,
    Type,
)

# Add custom type hints
Sasa: Type = Dict[str, Dict[str, freesasa.ResidueArea]]


def sasa_from_file(file: Union[str, pathlib.Path]) -> Sasa:
    """Get the freesasa.Result.residueAreas() dictionary
    obtained after parsing a PDB file to a freesasa.Structure
    and calling fresasa.calc() on it.
    """
    if isinstance(file, str):
        file = pathlib.Path(file)
    elif isinstance(file, pathlib.Path):
        pass
    else:
        raise TypeError("Invalid argument type. File should be 'str' or pathlib.Path")

    if not file.exists():
        raise FileNotFoundError(f"File {file.absolute().as_posix()} does not exist.")

    _struct = freesasa.Structure(file.absolute().as_posix())
    _sasa = freesasa.calc(_struct)

    return _sasa.residueAreas()


def delta_r_sasa(monomer: Sasa, complex: Sasa, chain: str, residue: str):
    """Return the difference calculated as:
    monomer.chain.residue.rSASA - complex.chain.residue.rSASA"""
    return monomer[chain][residue].relativeTotal - complex[chain][residue].relativeTotal


def residue_surface(sasa: Sasa, threshold: Optional[float] = 0.25):
    """ """
    liste_surface = []
    for chain in sasa:
        for residue in sasa[chain]:
            if sasa[chain][residue].relativeTotal > 0.25:
                liste_surface.append(residue)
    return liste_surface


def residue_interface(monomer, complex):
    """ """
    liste_interface = []
    for chain in monomer:
        for residue in monomer[chain]:
            if delta_r_sasa(monomer, complex, chain, residue) > 0:
                liste_interface.append(residue)
    return liste_interface