import sys
from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient

EXAMPLES_FOLDER = Path(__file__).parent.parent / "examples"

NOTEBOOKS = [
    "compare_molecules.ipynb",
    "phosgenation.ipynb",
    #"werner.ipynb",
    "visualization.ipynb",
]


@pytest.mark.example
@pytest.mark.parametrize("notebook", NOTEBOOKS)
def test_notebook(notebook):
    notebook_path = EXAMPLES_FOLDER / notebook

    with Path(notebook_path).open() as f:
        executed_nb = nbformat.read(f, as_version=4)

    client = NotebookClient(
        executed_nb,
        timeout=600,
        kernel_name="python",
        kernel_path=Path(sys.executable).resolve(),
    )
    client.execute()

    with Path(notebook_path).open() as f:
        nb = nbformat.read(f, as_version=4)

    for ex_cell, cell in zip(executed_nb.cells, nb.cells):
        assert ex_cell.cell_type == cell.cell_type

        if ex_cell.cell_type == "code":
            assert ex_cell.source == cell.source

            for ex_output, output in zip(ex_cell.outputs, cell.outputs):
                assert ex_output.output_type == output.output_type

                if ex_output.output_type == "stream":
                    assert ex_output.name == output.name
                    assert ex_output.text == output.text

                elif ex_output.output_type == "execute_result":
                    assert ex_output.data == output.data
                    assert ex_output.metadata == output.metadata

                elif ex_output.output_type == "display_data":
                    assert ex_output.data == output.data
                    assert ex_output.metadata == output.metadata

                elif ex_output.output_type == "error":
                    assert ex_output.ename == output.ename
                    assert ex_output.evalue == output.evalue
