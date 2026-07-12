from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def test_preprocessor(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    input_file = root / "examples" / "single_element_template.inp"
    config = json.loads((root / "examples" / "config_template.json").read_text())
    config["node_centers"] = [2.5]
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    output = tmp_path / "processed.inp"
    subprocess.run(
        [sys.executable, str(root / "preprocessing" / "frp_grid_preprocess.py"),
         str(input_file), str(config_file), "-o", str(output)],
        check=True,
    )
    text = output.read_text()
    assert "*USER ELEMENT, TYPE=U1" in text
    assert "*UEL PROPERTY, ELSET=UEL_E_1" in text
    assert "*ELEMENT, TYPE=U1, ELSET=UEL_E_1" in text
    assert "TYPE=U1P" not in text
