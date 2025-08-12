from pathlib import Path
def ensure_dir(p: Path):
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
