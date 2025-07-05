import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pytest

from stereomolgraph.cartesian import Geometry


class TestGeometry:
    @pytest.mark.parametrize(
        "xyz_content, comment, coords",
        [
            (  # with comment
                textwrap.dedent("""\
                    2
                        this is a comment comment with whitespace
                    C                 -3.7    0.02    0.2
                    C                 -3.1   -1.1   -0.2
                """),
                "this is a comment comment with whitespace",
                [[-3.7, 0.02, 0.2], [-3.1, -1.1, -0.2]],
            ),
            (  # with empty lines at end
                textwrap.dedent("""\
                    2

                    C                 -3.7    0.02    0.2
                    C                  0.0    0.0     0.0

                """),
                "",
                [[-3.7, 0.02, 0.2], [0.0, 0.0, 0.0]],
            ),
        ],
        ids=["with comment", "empty lines at end"],
    )
    def test_from_xyz_file(self, xyz_content, comment, coords):
        with tempfile.NamedTemporaryFile("w+", delete=False) as tmp:
            tmp.write(xyz_content)
            tmp.flush()
            fake_path = Path(tmp.name)

            geo = Geometry.from_xyz_file(fake_path)

        np.testing.assert_equal(geo.coords, coords)
