from pathlib import Path
import subprocess
import sys


def test_parity_script_runs_without_fixtures(tmp_path):
    root = Path(__file__).resolve().parents[1]
    empty_fixtures = tmp_path / "empty_fixtures"
    empty_fixtures.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "run_parity.py"),
            "--fixtures-root",
            str(empty_fixtures),
        ],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "No parity fixtures found" in result.stdout
