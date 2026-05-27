"""
run_g_factor.py
Run extract_g_factor.R via Rscript.
"""

from pathlib import Path
import subprocess
import shutil
import sys

BASE_DIR = Path(__file__).resolve().parent
R_SCRIPT = BASE_DIR / "extract_g_factor.R"

RSCRIPT_CANDIDATES = [
    r"R\R-4.5.1\bin\Rscript.exe",
    r"R\R-4.4.1\bin\Rscript.exe",
    "Rscript",
]

def find_rscript():
    for r in RSCRIPT_CANDIDATES:
        if shutil.which(r) or Path(r).exists():
            return r
    return None

def main():
    rscript = find_rscript()
    if rscript is None:
        sys.exit("ERROR: Rscript not found.")

    result = subprocess.run(
        [rscript, "--vanilla", str(R_SCRIPT)],
        cwd=BASE_DIR,
    )

    if result.returncode != 0:
        sys.exit(f"Rscript failed ({result.returncode})")

    print("Done.")

if __name__ == "__main__":
    main()