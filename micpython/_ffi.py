# -*- coding: utf-8 -*-
"""micpython 的 ctypes FFI 绑定层（零第三方依赖，仅标准库，Py3.8+）。

按 milp-ffi（WP1，``milp-ffi/src/api.rs`` 的 22 个 ``#[no_mangle] pub
extern "C"`` 符号）的 C ABI 契约工作：

- **符号绑定**：``EXPECTED_SYMBOLS`` 列出全部 22 个符号名；装载后逐一
  ``hasattr`` 校验并按 ``_SIGNATURES`` 显式声明 argtypes/restype（禁默认
  转换；``milp_ffi_model_destroy`` 返回 void，restype=None）；
- **错误处理**：``check(rc)`` 对 rc!=0 经 ``milp_ffi_last_error`` 取 UTF-8
  粘滞消息（errors='replace'），按 ``_EXC_BY_CODE``（覆盖 1..=11 全码，
  1/8→ArgumentError）抛对应异常；
- **消息缓冲契约**：``milp_ffi_last_error(buf, cap)`` 返回所需字节数（不含
  NUL）；buf=NULL/cap=0 只报长度；否则写 min(cap-1,len) 字节并补 NUL。
  所有缓冲归调用方（本模块）所有，FFI 侧绝不跨界 malloc/free——本模块用
  ``ctypes.create_string_buffer`` 分配，4096 起步，返回长度超过缓冲则二次
  调用重取；
- **DLL 装载**（惰性单例，``import micpython`` 无副作用）：环境变量
  ``MIC_FFI_LIB`` 非空→视为库文件路径；否则以本包目录上一级为工作区根，
  依次探测 ``root/milp-ffi/target/release/`` 下 ``milp_ffi.dll`` →
  ``libmilp_ffi.so`` → ``libmilp_ffi.dylib``；全失败抛 ``FfiLoadError``
  （消息列全部候选路径 + 构建提示）。装载成功后校验
  ``milp_ffi_abi_version() == 1``。``_reset_cache()`` 仅供测试；
- **GIL**：ctypes 在 CPython 下自动释放 GIL，长求解不阻塞解释器。

异常层次定义在本模块（而非 core.py）：``check``/``_EXC_BY_CODE`` 在导入时
即需引用这些类，若放在 core.py 会造成循环导入；core.py 再导出，用户侧入口
``from micpython import ModelError`` 不变。
"""

import ctypes
import os

ABI_VERSION = 1

# last_error 初始缓冲（字节）；返回长度超过缓冲时按需扩容重取。
_ERROR_BUFFER_SIZE = 4096
# format_number 输出上限（"{:e}" 最长约 24 字节，64 富余）。
_FORMAT_BUFFER_SIZE = 64

_LIB_FILE_NAMES = ('milp_ffi.dll', 'libmilp_ffi.so', 'libmilp_ffi.dylib')

_BUILD_HINT = ('提示：先在工作区根执行 '
               'cargo build --manifest-path milp-ffi/Cargo.toml --release，'
               '或设 MIC_FFI_LIB 指向已构建库')


# ===========================================================================
# 异常层次（core.py 再导出；.code 为 FFI 错误码，纯 Python 侧错误 code=0）
# ===========================================================================

class MicError(Exception):
    """micpython 统一异常基类（.message / .code）。"""

    def __init__(self, message='', code=0):
        Exception.__init__(self, message)
        self.message = message
        self.code = code


class FfiLoadError(MicError):
    """动态库装载失败（路径不存在 / 无法加载 / 符号缺失 / ABI 版本不符）。"""


class ArgumentError(MicError, ValueError):
    """FFI 码 1/8/9：空指针（公共 API 不可达）、UTF-8 编码失败（公共 API
    不可达：str.encode 已先行抛 UnicodeEncodeError）、参数值域/类型非法。

    亦为 ``ValueError`` 子类（B3）：『参数值非法』按 Python 惯例可被
    ``except ValueError`` 捕获——含内嵌 NUL 的名字/字符串值即经此抛出
    （core.py ``_check_no_nul``）。多重继承只放宽可捕获面，既有
    ``except ArgumentError`` / ``except MicError`` 语义不变。"""


class ModelError(MicError):
    """FFI 码 2：模型非法（变量重名、撞保留字、lb>ub、Binary 界越界等）。"""


class NumericalError(MicError):
    """FFI 码 3：数值异常（终检失败、NaN/Inf 等）。"""


class ParseError(MicError):
    """FFI 码 4：解析错误（FFI 层保留映射完整性）。"""


class IoError(MicError):
    """FFI 码 5：IO 错误（FFI 层保留映射完整性）。"""


class InternalError(MicError):
    """FFI 码 6：内部不变量被破坏。"""


class PanicError(MicError):
    """FFI 码 7：Rust panic 被 catch_unwind 捕获；句柄自此中毒（仅 destroy 可用）。"""


class UnknownParameterError(MicError):
    """FFI 码 10：参数名不存在。"""


class NoSolutionError(MicError):
    """FFI 码 11：未求解或本次求解无解时取解/统计/界。"""


# 错误码 → 异常类（全码覆盖 1..=11；1/8 均映射 ArgumentError）。
_EXC_BY_CODE = {
    1: ArgumentError,   # ERR_NULL（公共 API 不可达：ctypes 传参不为 NULL）
    2: ModelError,      # ERR_INVALID_MODEL
    3: NumericalError,  # ERR_NUMERICAL
    4: ParseError,      # ERR_PARSE
    5: IoError,         # ERR_IO
    6: InternalError,   # ERR_INTERNAL
    7: PanicError,      # ERR_PANIC
    8: ArgumentError,   # ERR_UTF8（公共 API 不可达：str.encode('utf-8') 先行抛出）
    9: ArgumentError,   # ERR_ARG（枚举越界/值域非法/lb 非有限等）
    10: UnknownParameterError,  # ERR_UNKNOWN_PARAM
    11: NoSolutionError,        # ERR_NO_SOLUTION
}


# ===========================================================================
# 22 个符号与 C 签名（milp-ffi/src/api.rs 与其 README 函数总表为准）
# ===========================================================================

EXPECTED_SYMBOLS = (
    'milp_ffi_abi_version',
    'milp_ffi_last_error',
    'milp_ffi_format_number',
    'milp_ffi_model_create',
    'milp_ffi_model_destroy',
    'milp_ffi_model_num_vars',
    'milp_ffi_model_num_constraints',
    'milp_ffi_add_var',
    'milp_ffi_set_var_bounds',
    'milp_ffi_add_constraint',
    'milp_ffi_set_objective',
    'milp_ffi_set_param_real',
    'milp_ffi_set_param_int',
    'milp_ffi_set_param_str',
    'milp_ffi_solve',
    'milp_ffi_get_status',
    'milp_ffi_get_objective',
    'milp_ffi_get_bound',
    'milp_ffi_get_gap_percent',
    'milp_ffi_get_solution',
    'milp_ffi_get_var_value',
    'milp_ffi_get_stat',
)

_c_char_p = ctypes.c_char_p
_c_void_p = ctypes.c_void_p
_c_int32 = ctypes.c_int32
_c_uint32 = ctypes.c_uint32
_c_int64 = ctypes.c_int64
_c_uint64 = ctypes.c_uint64
_c_double = ctypes.c_double
_p_char = ctypes.POINTER(ctypes.c_char)      # char* 出参缓冲（create_string_buffer）
_p_int32 = ctypes.POINTER(ctypes.c_int32)    # int32_t* 出参
_p_uint32 = ctypes.POINTER(ctypes.c_uint32)  # uint32_t* 出参/入参数组
_p_int64 = ctypes.POINTER(ctypes.c_int64)    # int64_t* 出参（get_stat）
_p_double = ctypes.POINTER(ctypes.c_double)  # double* 出参/入参数组

# name -> (argtypes, restype)，与 milp-ffi README「函数总表」逐行对应：
#   uint32_t milp_ffi_abi_version(void)
#   int64_t  milp_ffi_last_error(char *buf, uint64_t cap)
#   int64_t  milp_ffi_format_number(double v, char *buf, uint64_t cap)
#   void*    milp_ffi_model_create(const char *name, int32_t sense)
#   void     milp_ffi_model_destroy(void *model)                     ← 无错误码
#   int32_t  milp_ffi_model_num_vars(void*, uint32_t *out)
#   int32_t  milp_ffi_model_num_constraints(void*, uint32_t *out)
#   int32_t  milp_ffi_add_var(void*, const char *name, int32_t vtype,
#                             double lb, double ub, uint32_t *index_out)
#   int32_t  milp_ffi_set_var_bounds(void*, uint32_t index, double lb, double ub)
#   int32_t  milp_ffi_add_constraint(void*, const char *name_or_null,
#                             const uint32_t *var_idx, const double *coeffs,
#                             uint32_t count, int32_t relop, double rhs,
#                             uint32_t *index_out)
#   int32_t  milp_ffi_set_objective(void*, const uint32_t *var_idx,
#                             const double *coeffs, uint32_t count,
#                             int32_t sense, double offset)
#   int32_t  milp_ffi_set_param_real(void*, const char *name, double value)
#   int32_t  milp_ffi_set_param_int(void*, const char *name, int64_t value)
#   int32_t  milp_ffi_set_param_str(void*, const char *name, const char *value)
#   int32_t  milp_ffi_solve(void*)
#   int32_t  milp_ffi_get_status(void*, int32_t *out)
#   int32_t  milp_ffi_get_objective(void*, double *out)
#   int32_t  milp_ffi_get_bound(void*, double *out)
#   int32_t  milp_ffi_get_gap_percent(void*, double *out)
#   int32_t  milp_ffi_get_solution(void*, double *out, uint32_t cap,
#                             uint32_t *written_or_null)
#   int32_t  milp_ffi_get_var_value(void*, uint32_t index, double *out)
#   int32_t  milp_ffi_get_stat(void*, int32_t stat_id, int64_t *out)
_SIGNATURES = {
    'milp_ffi_abi_version': ([], _c_uint32),
    'milp_ffi_last_error': ([_p_char, _c_uint64], _c_int64),
    'milp_ffi_format_number': ([_c_double, _p_char, _c_uint64], _c_int64),
    'milp_ffi_model_create': ([_c_char_p, _c_int32], _c_void_p),
    'milp_ffi_model_destroy': ([_c_void_p], None),
    'milp_ffi_model_num_vars': ([_c_void_p, _p_uint32], _c_int32),
    'milp_ffi_model_num_constraints': ([_c_void_p, _p_uint32], _c_int32),
    'milp_ffi_add_var': ([_c_void_p, _c_char_p, _c_int32, _c_double, _c_double,
                          _p_uint32], _c_int32),
    'milp_ffi_set_var_bounds': ([_c_void_p, _c_uint32, _c_double, _c_double],
                                _c_int32),
    'milp_ffi_add_constraint': ([_c_void_p, _c_char_p, _p_uint32, _p_double,
                                 _c_uint32, _c_int32, _c_double, _p_uint32],
                                _c_int32),
    'milp_ffi_set_objective': ([_c_void_p, _p_uint32, _p_double, _c_uint32,
                                _c_int32, _c_double], _c_int32),
    'milp_ffi_set_param_real': ([_c_void_p, _c_char_p, _c_double], _c_int32),
    'milp_ffi_set_param_int': ([_c_void_p, _c_char_p, _c_int64], _c_int32),
    'milp_ffi_set_param_str': ([_c_void_p, _c_char_p, _c_char_p], _c_int32),
    'milp_ffi_solve': ([_c_void_p], _c_int32),
    'milp_ffi_get_status': ([_c_void_p, _p_int32], _c_int32),
    'milp_ffi_get_objective': ([_c_void_p, _p_double], _c_int32),
    'milp_ffi_get_bound': ([_c_void_p, _p_double], _c_int32),
    'milp_ffi_get_gap_percent': ([_c_void_p, _p_double], _c_int32),
    'milp_ffi_get_solution': ([_c_void_p, _p_double, _c_uint32, _p_uint32],
                              _c_int32),
    'milp_ffi_get_var_value': ([_c_void_p, _c_uint32, _p_double], _c_int32),
    'milp_ffi_get_stat': ([_c_void_p, _c_int32, _p_int64], _c_int32),
}

# 导入期自检：符号名表与签名表必须一一对应（不一致即本文件内部错误）。
assert set(_SIGNATURES) == set(EXPECTED_SYMBOLS), \
    '_SIGNATURES 与 EXPECTED_SYMBOLS 不一致'


# ===========================================================================
# 装载（惰性单例）
# ===========================================================================

_lib_cache = None


def _default_candidate_paths():
    """默认候选：<本包目录上一级>/milp-ffi/target/release/ 下按序探测。"""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    release = os.path.join(root, 'milp-ffi', 'target', 'release')
    return [os.path.join(release, name) for name in _LIB_FILE_NAMES]


def _bind_symbols(lib):
    """逐一校验 22 个符号存在并声明 argtypes/restype。"""
    missing = [name for name in EXPECTED_SYMBOLS if not hasattr(lib, name)]
    if missing:
        raise FfiLoadError('动态库缺少必需符号（%d 个）：%s'
                           % (len(missing), ', '.join(missing)))
    for name in (_SIGNATURES):
        argtypes, restype = _SIGNATURES[name]
        fn = getattr(lib, name)
        fn.argtypes = argtypes
        fn.restype = restype


def load_library():
    """按 DLL 装载规则定位并装载 milp-ffi，返回已绑定符号的 CDLL。

    ① ``MIC_FFI_LIB`` 非空 → 视为库文件路径（须存在且可加载）；
    ② 否则探测默认候选（dll → so → dylib，依次 exists()）；
    ③ 全失败抛 ``FfiLoadError``（列全部候选路径与失败原因 + 构建提示）；
    ④ 校验 ``milp_ffi_abi_version() == 1``；
    ⑤ 校验并绑定全部 22 个符号。
    """
    env_path = os.environ.get('MIC_FFI_LIB', '')
    if env_path:
        candidates = [env_path]
    else:
        candidates = _default_candidate_paths()
    failures = []
    lib = None
    loaded_path = None
    for path in candidates:
        if not os.path.isfile(path):
            failures.append('%s：文件不存在' % path)
            continue
        try:
            lib = ctypes.CDLL(path)
            loaded_path = path
            break
        except OSError as exc:
            failures.append('%s：%s' % (path, exc))
            lib = None
    if lib is None:
        raise FfiLoadError('无法装载 milp-ffi 动态库。尝试过的候选路径与失败原因：\n'
                           '  - %s\n%s'
                           % ('\n  - '.join(failures), _BUILD_HINT))
    _bind_symbols(lib)
    version = lib.milp_ffi_abi_version()
    if version != ABI_VERSION:
        raise FfiLoadError('ABI 版本不匹配：%s 返回 %d，micpython 期望 %d。%s'
                           % (loaded_path, version, ABI_VERSION, _BUILD_HINT))
    return lib


def _lib():
    """模块级惰性单例：首次调用时装载，之后复用；失败每次都抛 FfiLoadError。"""
    global _lib_cache
    if _lib_cache is None:
        _lib_cache = load_library()
    return _lib_cache


def _reset_cache():
    """仅供测试：清空单例，强制下次重新装载（不卸载 OS 层的 DLL 句柄）。"""
    global _lib_cache
    _lib_cache = None


# ===========================================================================
# 错误消息读取与错误码分发
# ===========================================================================

def _copy_out_text(fn, *args):
    """按缓冲契约执行一次「写入 + 需要时扩容重取」，返回 UTF-8 文本。

    fn 满足 milp_ffi 缓冲契约：末参为 (buf, cap)，返回所需字节数（不含
    NUL）；buf 写 min(cap-1, len) 字节并补 NUL。返回负值（panic 兜底）按
    空文本处理。
    """
    buf = ctypes.create_string_buffer(_ERROR_BUFFER_SIZE)
    need = fn(*(list(args) + [buf, _ERROR_BUFFER_SIZE]))
    if isinstance(need, int) and 0 <= _ERROR_BUFFER_SIZE - 1 < need:
        size = need + 1
        buf = ctypes.create_string_buffer(size)
        fn(*(list(args) + [buf, size]))
    return buf.value.decode('utf-8', errors='replace')


def last_error_message():
    """读取线程局部粘滞错误消息（UTF-8，errors='replace'）；未装载返回 ''。"""
    lib = _lib_cache
    if lib is None:
        return ''
    return _copy_out_text(lib.milp_ffi_last_error)


def check(rc):
    """rc != 0 时取粘滞消息并按 _EXC_BY_CODE 抛对应异常；rc == 0 原样返回。"""
    if rc == 0:
        return rc
    message = last_error_message()
    if not message:
        message = 'FFI 返回错误码 %d（未携带错误消息）' % rc
    exc_class = _EXC_BY_CODE.get(rc)
    if exc_class is None:
        exc_class = MicError
    raise exc_class(message, rc)


def format_number(value):
    """经 FFI 格式化数值：与 milp::search_ctrl::format_number / CLI objective
    行逐字节一致（0<|v|<1e-6 或 |v|≥1e15 用科学记数，否则定点 6 位）。"""
    return _copy_out_text(_lib().milp_ffi_format_number, float(value))
