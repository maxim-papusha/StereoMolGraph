from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping

    
Element = int
ElementLike = Element | str | int


SYMBOLS: dict[Element, str] = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
    10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S',
    17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V',
    24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb',
    38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru',
    45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb',
    52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce',
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
    116: 'Lv', 117: 'Ts', 118: 'Og',
    }

_COVALENT_RADII: dict[Element, float] = {
    1: 0.32, 2: 0.46, 3: 1.33, 4: 1.02, 5: 0.85, 6: 0.75, 7: 0.71, 8: 0.63,
    9: 0.64, 10: 0.67, 11: 1.55, 12: 1.39, 13: 1.26, 14: 1.16, 15: 1.11,
    16: 1.03, 17: 0.99, 18: 0.96, 19: 1.96, 20: 1.71, 21: 1.48, 22: 1.36,
    23: 1.34, 24: 1.22, 25: 1.19, 26: 1.16, 27: 1.11, 28: 1.1, 29: 1.12,
    30: 1.18, 31: 1.24, 32: 1.21, 33: 1.21, 34: 1.16, 35: 1.14, 36: 1.17,
    37: 2.1, 38: 1.85, 39: 1.63, 40: 1.54, 41: 1.47, 42: 1.38, 43: 1.28,
    44: 1.25, 45: 1.25, 46: 1.2, 47: 1.28, 48: 1.36, 49: 1.42, 50: 1.4,
    51: 1.4, 52: 1.36, 53: 1.33, 54: 1.31, 55: 2.32, 56: 1.96, 57: 1.8,
    58: 1.63, 59: 1.76, 60: 1.74, 61: 1.73, 62: 1.72, 63: 1.68, 64: 1.69,
    65: 1.68, 66: 1.67, 67: 1.66, 68: 1.65, 69: 1.64, 70: 1.7, 71: 1.62,
    72: 1.52, 73: 1.46, 74: 1.37, 75: 1.31, 76: 1.29, 77: 1.22, 78: 1.23,
    79: 1.24, 80: 1.33, 81: 1.44, 82: 1.44, 83: 1.51, 84: 1.45, 85: 1.47,
    86: 1.42, 87: 2.23, 88: 2.01, 89: 1.86, 90: 1.75, 91: 1.69, 92: 1.7,
    93: 1.71, 94: 1.72, 95: 1.66, 96: 1.66, 97: 1.68, 98: 1.68, 99: 1.65,
    100: 1.67, 101: 1.73, 102: 1.76, 103: 1.61, 104: 1.57, 105: 1.49,
    106: 1.43, 107: 1.41, 108: 1.34, 109: 1.29, 110: 1.28, 111: 1.21,
    112: 1.22, 113: 1.36, 114: 1.43, 115: 1.62, 116: 1.75, 117: 1.65,
    118: 1.57
}


# put all elements in a table searchable by atomic number and symbol

_PERIODIC_TABLE:dict[int|str, Element] = {}
    
for atomic_nr, symbol in SYMBOLS.items():
    _PERIODIC_TABLE[symbol] = atomic_nr
    _PERIODIC_TABLE[atomic_nr] = atomic_nr
    _PERIODIC_TABLE[symbol.upper()] = atomic_nr
    _PERIODIC_TABLE[symbol.lower()] = atomic_nr


PERIODIC_TABLE: Mapping[str|int|Element, Element] = MappingProxyType(
                                                            _PERIODIC_TABLE)
"Mapping of atomic numbers and symbols to Element objects."


COVALENT_RADII: Mapping[Element, float] = MappingProxyType(_COVALENT_RADII)
"""Covalent radii of elements in Angstrom.
(Pekka Pyykkö and Michiko Atsumi.
Molecular Double-Bond Covalent Radii for Elements Li-E112.
Chemistry - A European Journal, 15(46):12770–12779, nov 2009.
doi:10.1002/chem.200901472.)"""

assert PERIODIC_TABLE["C"] in PERIODIC_TABLE.keys()