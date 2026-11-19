"""Export Sphinx-generated Markdown docs into a versioned repository mirror.

The original Sphinx sources in docs/ remain authoritative.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REPO_LLM_INDEX = """# StereoMolGraph

> Stereochemistry-aware molecular and reaction graphs for Python.

Agent-readable Markdown generated from Sphinx on commits to the `doc_md` branch.

## Start here

- [Documentation overview](docs_md/index.md)
- [Agent usage guide](docs_md/agent-guide.md)
- [Installation](docs_md/installation.md)

## API reference

- [MolGraph](docs_md/reference/graphs/molgraph.md)
- [StereoMolGraph](docs_md/reference/graphs/stereomolgraph.md)
- [CondensedReactionGraph](docs_md/reference/graphs/condensedreactiongraph.md)
- [StereoCondensedReactionGraph](docs_md/reference/graphs/stereocondensedreactiongraph.md)
- [Algorithms](docs_md/reference/algorithms/index.md)

## Complete documentation

- [Page index](docs_md/llms.txt)
- [Full documentation](docs_md/llms-full.txt)
- [Browsable documentation](https://stereomolgraph.readthedocs.io/)
"""

MIRROR_README = """# Generated Markdown documentation

**Do not edit files in this directory manually.** The documentation sources
are in [`docs/`](../docs/) and Python docstrings. GitHub Actions regenerates this
mirror on every push to `doc_md`.

- [Documentation overview](index.md)
- [Agent guide](agent-guide.md)
- [LLM index](llms.txt)
- [Full documentation](llms-full.txt)

Regenerate locally from the repository root:

```bash
python -m pip install -e '.[docs]'
python -m sphinx -b html docs docs/_build/html
python scripts/export_markdown_docs.py
```
"""


def export(build_dir: Path, output_dir: Path, repo_index: Path) -> int:
    """Copy Markdown and linked images, removing obsolete generated output."""
    build_dir = build_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir == build_dir or build_dir in output_dir.parents:
        raise ValueError("Output directory must be outside the Sphinx build")

    required = ("index.md", "llms.txt", "llms-full.txt",
                "reference/graphs/stereomolgraph.md")
    missing = [name for name in required if not (build_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(
            f"Incomplete Markdown build at {build_dir}; missing: {', '.join(missing)}"
        )

    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    pages = sorted(
        p for p in build_dir.rglob("*.md")
        if "_build" not in p.relative_to(build_dir).parts
    )
    for page in pages:
        destination = output_dir / page.relative_to(build_dir)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(page, destination)

    for filename in ("llms.txt", "llms-full.txt"):
        shutil.copy2(build_dir / filename, output_dir / filename)

    for image_dir in ("_images", "_static/img"):
        source = build_dir / image_dir
        if source.is_dir():
            shutil.copytree(source, output_dir / image_dir, dirs_exist_ok=True)

    (output_dir / "README.md").write_text(MIRROR_README, encoding="utf-8")
    repo_index.write_text(REPO_LLM_INDEX, encoding="utf-8")
    return len(pages)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, default=ROOT / "build/docs/html")
    parser.add_argument("--output", type=Path, default=ROOT / "docs_md")
    args = parser.parse_args()
    count = export(args.build, args.output, ROOT / "llms.txt")
    print(f"Exported {count} Markdown pages to {args.output}")


if __name__ == "__main__":
    main()
