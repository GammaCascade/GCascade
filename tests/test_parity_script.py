from pathlib import Path
import subprocess
import sys


def test_parity_script_runs_without_fixtures():
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, str(root / "scripts" / "run_parity.py")],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "No parity fixtures found" in result.stdout
