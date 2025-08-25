
from collections import Counter
from collections.abc import Iterable, Hashable

import numpy as np

import math

def jaccard_distance_np(a: np.ndarray, b: np.ndarray) -> float:
    a = np.unique(a)
    b = np.unique(b)
    intersection = np.intersect1d(a, b).size
    union = np.union1d(a, b).size
    if union == 0:
        return 1.0
    return 1 - (intersection / union)

def jaccard_distance(a: Iterable[Hashable]|Counter,
                     b: Iterable[Hashable]|Counter) -> float:
    """Compute the Jaccard distance between two sets.
    Jaccard distance = 1 - (|A ∩ B| / |A ∪ B|)
    """
    a_counter = Counter(a)
    b_counter = Counter(b)

    union = sum((a_counter | b_counter).values())
    if union == 0:
        return 1.0
    intersection = sum((a_counter & b_counter).values())
    return 1 - (intersection / union)


def cosine_distance(a: Iterable[Hashable]|Counter,
                    b: Iterable[Hashable]|Counter) -> float:
    """Compute the cosine distance between two containers.
    Treated as frequency vectors
    cosine distance = 1 - (sum(a_i * b_i) / (sqrt(sum(a_i²)) * sqrt(sum(b_i²))))
    """
    a = Counter(a)
    b = Counter(b)
    keys = set(a.keys()) | set(b.keys())
    dot_product = sum(a.get(k, 0) * b.get(k, 0) for k in keys)
    mag_a = math.sqrt(sum(v**2 for v in a.values()))
    mag_b = math.sqrt(sum(v**2 for v in b.values()))
    if mag_a == 0 or mag_b == 0:
        return 1.0
    return 1 - (dot_product / (mag_a * mag_b))

class FingerprintGenerator:
    ...

class ECFPGenerator:
    ...

class FCFPxCGenerator:
    ...

class MAPGenerator:
    ...

class MAPxCGenerator:
    ...

class MACCSGenerator:
    ...

