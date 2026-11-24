"""独立加载规划器动态库。首次求解时才加载。"""
import ctypes as C
import os
import sys
from functools import lru_cache
from pathlib import Path


@lru_cache(None)
def library():
    filename = ("planner_ffi.dll" if os.name == "nt" else
                "libplanner_ffi.dylib" if sys.platform == "darwin" else "libplanner_ffi.so")
    override = os.getenv("MIC_PLANNER_LIB")
    path = Path(override) if override else Path(__file__).resolve().parents[2] / "planner-ffi" / "target" / "release" / filename
    try:
        lib = C.CDLL(str(path))
        signatures = {
            "planner_abi_version": (C.c_uint32, []),
            "planner_create": (C.c_uint64, []),
            "planner_destroy": (None, [C.c_uint64]),
            "planner_run": (C.c_int32, [C.c_uint64, C.c_void_p, C.c_size_t]),
            "planner_size": (C.c_size_t, [C.c_uint64]),
            "planner_read": (C.c_int32, [C.c_uint64, C.c_void_p, C.c_size_t]),
        }
        for name, (ret, args) in signatures.items():
            func = getattr(lib, name)
            func.restype, func.argtypes = ret, args
        if lib.planner_abi_version() != 1:
            raise RuntimeError("规划器需要 ABI 版本 1")
    except (OSError, AttributeError) as error:
        raise RuntimeError("不能加载规划器 %s；请构建 planner-ffi 或设置 MIC_PLANNER_LIB" % path) from error
    return lib


def run(lib, handle, data):
    code = lib.planner_run(handle, data, len(data))
    size = lib.planner_size(handle)
    if size == C.c_size_t(-1).value:
        raise RuntimeError("规划器句柄无效")
    buffer = C.create_string_buffer(size)
    if lib.planner_read(handle, buffer, size) != 0:
        raise RuntimeError("读取规划器结果失败")
    value = buffer.raw.decode("utf-8")
    if code:
        raise ValueError(value)
    return value
