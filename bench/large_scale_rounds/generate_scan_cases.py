#!/usr/bin/env python3
"""生成固定种子的扫描留出集。

在仓库根目录运行：
    python3 bench/large_scale_rounds/generate_scan_cases.py
算例与原始清单位于 /tmp/mic-rounds-scan-holdouts，清单副本位于
bench/large_scale_rounds/results/scan-holdouts.json。
"""
import importlib.util
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path('/tmp/mic-rounds-scan-holdouts')
MANIFEST = ROOT / 'bench/large_scale_rounds/results/scan-holdouts.json'
SOURCE = ROOT / 'bench/experiments/generate_scan_holdouts.py'

spec = importlib.util.spec_from_file_location('generate_scan_holdouts', SOURCE)
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
generator.CASES = [
    ('scan2-sparse-wide', 2400, 500, 0.05, 2026111201),
    ('scan2-sparse-tall', 1000, 1500, 0.08, 2026111202),
    ('scan2-medium', 1600, 800, 0.30, 2026111203),
    ('scan2-dense-wide', 2000, 500, 1.0, 2026111204),
]
generator.generate(OUTPUT)
MANIFEST.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(OUTPUT / 'manifest.json', MANIFEST)
