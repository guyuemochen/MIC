"""构建零依赖 Java 8 JAR、Rust 动态库与 JNI 薄桥，产物统一放入 target。"""
import os
import argparse
from pathlib import Path
import shutil
import subprocess
import sys


HERE = Path(__file__).resolve().parent
FFI = HERE.parent / "planner-ffi"
TARGET = HERE / "target"


def run(args):
    subprocess.run([str(arg) for arg in args], check=True)


def jdk():
    if os.getenv("JAVA_HOME"):
        home = Path(os.environ["JAVA_HOME"])
    else:
        java = shutil.which("javac")
        if not java:
            raise SystemExit("需要安装 JDK 8 或更高版本，或设置 JAVA_HOME")
        home = Path(java).resolve().parent.parent
        if not (home / "include" / "jni.h").is_file() and sys.platform == "darwin":
            home = Path(subprocess.check_output(["/usr/libexec/java_home"], text=True).strip())
    if not (home / "include" / "jni.h").is_file():
        raise SystemExit("JAVA_HOME 需要指向包含 include/jni.h 的完整 JDK")
    return home


def build():
    home = jdk()
    classes = TARGET / "classes"
    classes.mkdir(parents=True, exist_ok=True)
    libdir = FFI / "target" / "release"
    run(["cargo", "build", "--release", "--manifest-path", FFI / "Cargo.toml",
         "--target-dir", FFI / "target"])
    suffix = ".exe" if os.name == "nt" else ""
    javac = home / "bin" / ("javac" + suffix)
    version = subprocess.check_output([str(javac), "-version"], stderr=subprocess.STDOUT, text=True)
    flags = ["-source", "8", "-target", "8"] if "1.8." in version else ["--release", "8"]
    run([javac] + flags + ["-encoding", "UTF-8", "-d", classes] + sorted((HERE / "src").rglob("*.java")))
    run([home / "bin" / ("jar" + suffix), "cf", TARGET / "micplanner.jar", "-C", classes, "org/micplanner"])
    platform = {"darwin": "darwin", "win32": "win32", "linux": "linux"}.get(sys.platform)
    if platform is None:
        raise SystemExit("目前支持 macOS、Linux 和 Windows")
    libname = {"darwin": "libplanner_ffi.dylib", "win32": "planner_ffi.dll",
               "linux": "libplanner_ffi.so"}[platform]
    shutil.copy2(str(libdir / libname), str(TARGET / libname))
    headers = [home / "include", home / "include" / platform, FFI / "include"]
    source = HERE / "native" / "bridge.c"
    if platform == "win32":
        # MSVC 从开发者命令行运行；GNU 工具链使用导入库自动链接。
        if (libdir / "planner_ffi.dll.lib").is_file():
            run([os.getenv("CC", "cl"), "/nologo", "/LD", "/O2"] +
                ["/I" + str(p) for p in headers] + [source, libdir / "planner_ffi.dll.lib",
                 "/Fo" + str(TARGET / "bridge.obj"), "/link", "/OUT:" + str(TARGET / "micplanner_jni.dll")])
            return
    output = {"win32": "micplanner_jni.dll", "darwin": "libmicplanner_jni.dylib",
              "linux": "libmicplanner_jni.so"}[platform]
    cmd = [os.getenv("CC", "cc"), "-std=c11", "-Wall", "-Wextra", "-O2"]
    if platform != "win32":
        cmd += ["-fPIC"]
    cmd += ["-dynamiclib" if platform == "darwin" else "-shared"]
    cmd += ["-I" + str(p) for p in headers] + [source, "-L" + str(libdir), "-lplanner_ffi"]
    if platform == "linux":
        cmd += ["-Wl,-rpath,$ORIGIN"]
    if platform == "darwin":
        cmd += ["-Wl,-rpath,@loader_path"]
    run(cmd + ["-o", TARGET / output])
    if platform == "darwin":
        # Rust 默认的安装名可能指向构建机绝对路径，使复制后的 JAR 无法加载。
        install_id = subprocess.check_output(["otool", "-D", str(TARGET / libname)], text=True).splitlines()[1]
        run(["install_name_tool", "-change", install_id, "@loader_path/" + libname, TARGET / output])
        run(["codesign", "--force", "--sign", "-", TARGET / output])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", action="store_true", help="构建后在所选 JDK 上运行绑定验收")
    options = parser.parse_args()
    build()
    if options.test:
        home = jdk()
        suffix = ".exe" if os.name == "nt" else ""
        classes = TARGET / "test-classes"
        classes.mkdir(exist_ok=True)
        jar = TARGET / "micplanner.jar"
        run([home / "bin" / ("javac" + suffix), "-encoding", "UTF-8", "-cp", jar,
             "-d", classes, HERE / "tests" / "BindingTest.java"])
        run([home / "bin" / ("java" + suffix), "-cp", os.pathsep.join(map(str, [jar, classes])), "BindingTest"])
