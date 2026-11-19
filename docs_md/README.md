# Generated Markdown documentation

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
