"""从仓库源码构建，逐项核对 Rust、Python 和 Java 独立建立的四个模型。"""

import difflib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "planner-java"


def run(args, env=None, capture=False):
    return subprocess.run([str(arg) for arg in args], cwd=str(ROOT), env=env,
                          check=True, stdout=subprocess.PIPE if capture else None,
                          encoding="utf8").stdout


class PlannerBindingsParity(unittest.TestCase):
    def test_four_business_models(self):
        run([sys.executable, JAVA / "build.py"])
        expected = run(["cargo", "run", "--quiet", "--release", "--manifest-path",
                        ROOT / "planner/Cargo.toml", "--example", "business_models"], capture=True)
        self.assertEqual([s.split("\t")[0] for s in expected.splitlines()],
                         ["roster", "timetable", "delivery", "factory"])
        target = JAVA / "target"
        jar = target / "micplanner.jar"
        suffix = "dll" if os.name == "nt" else "dylib" if sys.platform == "darwin" else "so"
        prefix = "" if os.name == "nt" else "lib"
        env = dict(os.environ,
                   MIC_PLANNER_LIB=str(target / (prefix + "planner_ffi." + suffix)),
                   MIC_PLANNER_JNI=str(target / (prefix + "micplanner_jni." + suffix)))
        python_output = run([sys.executable, ROOT / "micpython/examples/planner_cases.py"],
                            env=env, capture=True)
        home = os.environ.get("JAVA_HOME")
        bin_dir = Path(home) / "bin" if home else Path("")
        with tempfile.TemporaryDirectory() as classes:
            run([bin_dir / "javac", "-encoding", "UTF-8", "-cp", jar, "-d", classes,
                 JAVA / "examples/PlannerCases.java", JAVA / "tests/BindingTest.java"])
            classpath = os.pathsep.join([classes, str(jar)])
            run([bin_dir / "java", "-cp", classpath, "BindingTest"], env=env)
            java_output = run([bin_dir / "java", "-cp", classpath, "PlannerCases"],
                              env=env, capture=True)
        for language, actual in (("Python", python_output), ("Java", java_output)):
            with self.subTest(language=language):
                self.assertEqual(expected, actual, "".join(difflib.unified_diff(
                    expected.splitlines(True), actual.splitlines(True),
                    fromfile="Rust", tofile=language)))


if __name__ == "__main__":
    unittest.main()
