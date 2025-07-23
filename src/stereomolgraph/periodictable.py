from __future__ import annotations

from dataclasses import dataclass
from functools import total_ordering
from types import MappingProxyType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any
    


@total_ordering
@dataclass
class Element:
    """
    represents chemical element

    :param symbol: atomic symbol
    :param atomic_nr: atomic number of element
    """

    symbol: str
    atomic_nr: int

    def __hash__(self) -> int:
        return self.atomic_nr

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, Element):
            return self.atomic_nr == other.atomic_nr
        else:
            return NotImplemented

    def __lt__(self, other: Any) -> bool:
        if isinstance(other, Element):
            return self.atomic_nr < other.atomic_nr
        else:
            return NotImplemented

    def __repr__(self) -> str:
        return self.symbol

    def __str__(self) -> str:
        return self.symbol

    def __format__(self, __format_spec: str) -> str:
        if __format_spec == 's':
            return self.symbol
        elif __format_spec == 'd':
            return str(self.atomic_nr)
        else:
            raise TypeError('unsupported format string passed to Element.'
                            '__format__')

# added underscore so that these are not imported with from constants import *
_SYMBOLS = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
            9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
            16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
            23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni',
            29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se',
            35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
            41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd',
            47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te',
            53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce',
            59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
            65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
            71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os',
            77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb',
            83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra',
            89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu',
            95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
            101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db',
            106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
            111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc',
            116: 'Lv', 117: 'Ts', 118: 'Og'}

_ATOMIC_NRS = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
               'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
               'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
               'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
               'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
               'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
               'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
               'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
               'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
               'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
               'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
               'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
               'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
               'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
               'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
               'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
               'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102,
               'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107,
               'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
               'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117,
               'Og': 118}

#    covalent radius in (Pekka Pyykkö and Michiko
#       Atsumi. Molecular Double-Bond Covalent Radii for Elements Li-E112.
#       Chemistry - A European Journal, 15(46):12770–12779, nov 2009.
#       URL: http://doi.wiley.com/10.1002/chem.200901472,
#       doi:10.1002/chem.200901472.)

_ELEMENT_COVALENT_RADII = {'H': 0.32, 'He': 0.46, 'Li': 1.33, 'Be': 1.02,
                           'B': 0.85, 'C': 0.75, 'N': 0.71, 'O': 0.63,
                           'F': 0.64, 'Ne': 0.67, 'Na': 1.55, 'Mg': 1.39,
                           'Al': 1.26, 'Si': 1.16, 'P': 1.11, 'S': 1.03,
                           'Cl': 0.99, 'Ar': 0.96, 'K': 1.96, 'Ca': 1.71,
                           'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22,
                           'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11, 'Ni': 1.1,
                           'Cu': 1.12, 'Zn': 1.18, 'Ga': 1.24, 'Ge': 1.21,
                           'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17,
                           'Rb': 2.1, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54,
                           'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25,
                           'Rh': 1.25, 'Pd': 1.2, 'Ag': 1.28, 'Cd': 1.36,
                           'In': 1.42, 'Sn': 1.4, 'Sb': 1.4, 'Te': 1.36,
                           'I': 1.33, 'Xe': 1.31, 'Cs': 2.32, 'Ba': 1.96,
                           'La': 1.8, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74,
                           'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68, 'Gd': 1.69,
                           'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65,
                           'Tm': 1.64, 'Yb': 1.7, 'Lu': 1.62, 'Hf': 1.52,
                           'Ta': 1.46, 'W': 1.37, 'Re': 1.31, 'Os': 1.29,
                           'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33,
                           'Tl': 1.44, 'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45,
                           'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01,
                           'Ac': 1.86, 'Th': 1.75, 'Pa': 1.69, 'U': 1.7,
                           'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66,
                           'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67,
                           'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57,
                           'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34,
                           'Mt': 1.29, 'Ds': 1.28, 'Rg': 1.21, 'Cn': 1.22,
                           'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75,
                           'Ts': 1.65, 'Og': 1.57}


# put all elements in a table searchable by atomic number and symbol
_PERIODIC_TABLE:dict[int|str|Element, Element] = {
    sym : Element(sym,
                  _ATOMIC_NRS[sym],)
    for sym in _ATOMIC_NRS}

_PERIODIC_TABLE.update({
    nr : _PERIODIC_TABLE[_SYMBOLS[nr]]
    for nr in _SYMBOLS})

_PERIODIC_TABLE.update({
    _PERIODIC_TABLE[_SYMBOLS[nr]] : _PERIODIC_TABLE[_SYMBOLS[nr]]
    for nr in _SYMBOLS})

_PERIODIC_TABLE.update({
    sym.upper() : _PERIODIC_TABLE[sym]
    for sym in _ATOMIC_NRS})

_PERIODIC_TABLE.update({
    sym.lower() : _PERIODIC_TABLE[sym]
    for sym in _ATOMIC_NRS})


# use MappingProxyType to make the public table immutable
PERIODIC_TABLE: Mapping[str|int|Element, Element] = MappingProxyType(
                                                            _PERIODIC_TABLE)

COVALENT_RADII: Mapping[Element, float] = MappingProxyType(
    {Element(sym, _ATOMIC_NRS[sym]): _ELEMENT_COVALENT_RADII[sym]
     for sym in _ELEMENT_COVALENT_RADII})