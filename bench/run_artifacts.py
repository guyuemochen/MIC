"""为每次基准运行分配独立的输出目录。"""

from pathlib import Path
from tempfile import mkdtemp


def new_run_dir(name: str) -> Path:
    root = Path(__file__).resolve().parent / "runs"
    root.mkdir(exist_ok=True)
    return Path(mkdtemp(prefix=f"{name}-", dir=root))
