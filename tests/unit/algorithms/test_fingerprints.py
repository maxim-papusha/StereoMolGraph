from collections import Counter

import pytest

from stereomolgraph.algorithms.fingerprints import (
    cosine_distance,
    jaccard_distance,
)


class TestJaccardDistance:
    def test_set_jaccard(self):
        assert jaccard_distance({1, 2}, {2, 3}) == pytest.approx(2 / 3)
        assert jaccard_distance(set(), set()) == 1.0  # Edge case

    def test_counter_jaccard(self):
        c1 = Counter("hello")  # h:1, e:1, l:2, o:1
        c2 = Counter("world")  # w:1, o:1, r:1, l:1, d:1
        # Intersection: o:1 + l:1 = 2
        # Union: h:1 + e:1 + l:2 + o:1 + w:1 + r:1 + d:1 = 8
        assert jaccard_distance(c1, c2) == pytest.approx(1 - (2 / 8))


class TestCosineDistance:
    def test_sets(self):
        assert cosine_distance({1, 2}, {2, 3}) == pytest.approx(0.5)
        assert cosine_distance(set(), {1}) == 1.0  # Empty set

    def test_counters(self):
        c1 = Counter("hello")  # h:1, e:1, l:2, o:1
        c2 = Counter("world")  # w:1, o:1, r:1, l:1, d:1
        assert cosine_distance(c1, c2) == pytest.approx(0.49290744716289014)

    def test_edge_cases(self):
        # Zero vector
        assert cosine_distance(
            Counter([0, 0]), Counter([1, 1])
        ) == pytest.approx(1.0)
        # Identical vectors
        assert cosine_distance(
            Counter([1, 2, 3]), Counter([1, 2, 3])
        ) == pytest.approx(0.0)
