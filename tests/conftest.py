# tests/conftest.py
import pytest
from pathlib import Path

def pytest_collection_modifyitems(items, config):
    """Auto-mark tests based on their folder using pathlib"""
    for item in items:
        test_path = Path(item.fspath)
        
        # Get the parts of the path relative to tests directory
        try:
            # Try to get path relative to tests directory
            rel_to_tests = test_path.relative_to(Path("tests"))
            folder_name = rel_to_tests.parts[0]  # First directory under tests/
            
            if folder_name == "unit":
                item.add_marker(pytest.mark.unit)
            elif folder_name == "hypothesis":
                item.add_marker(pytest.mark.hypothesis)
            #elif folder_name == "examples":
            #    item.add_marker(pytest.mark.examples)
                
        except ValueError:
            # Test not under tests/ directory, skip marking
            pass