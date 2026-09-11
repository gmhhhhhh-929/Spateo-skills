"""Exercise the packaged modules, not another checkout or installed core."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'skills/spatial-slice-quality-qc/scripts'))
from activate_candidate import activate
activate()
