# -*- coding: utf-8 -*-
"""micpython._ffi 装载层测试（全套件中唯一允许 import _ffi 的测试文件）。

覆盖（对应 WP3 契约）：
- 默认路径装载成功且 ``milp_ffi_abi_version() == 1``；
- ``EXPECTED_SYMBOLS`` 共 22 个且在 CDLL 上全部存在；
- ``_EXC_BY_CODE`` 覆盖 1..=11 每码且值均为 MicError 子类；
- ``MIC_FFI_LIB`` 指向不存在路径 → ``FfiLoadError`` 且消息含构建提示
  ``cargo build --manifest-path milp-ffi/Cargo.toml --release``；
- ``MIC_FFI_LIB`` 指向真实 DLL → 装载成功（环境变量 + ``_reset_cache``，
  tearDown 恢复环境并重置单例）。
"""

import ctypes
import os
import tempfile
import unittest

from micpython import FfiLoadError, MicError, _ffi

BUILD_HINT = 'cargo build --manifest-path milp-ffi/Cargo.toml --release'


class LoaderTests(unittest.TestCase):

    def setUp(self):
        self._old_env = os.environ.get('MIC_FFI_LIB')
        _ffi._reset_cache()

    def tearDown(self):
        if self._old_env is None:
            os.environ.pop('MIC_FFI_LIB', None)
        else:
            os.environ['MIC_FFI_LIB'] = self._old_env
        _ffi._reset_cache()

    def test_default_path_loads_with_abi_version_1(self):
        os.environ.pop('MIC_FFI_LIB', None)
        lib = _ffi._lib()
        self.assertIsInstance(lib, ctypes.CDLL)
        self.assertEqual(lib.milp_ffi_abi_version(), 1)

    def test_expected_symbols_count_and_presence(self):
        self.assertEqual(len(_ffi.EXPECTED_SYMBOLS), 22)
        self.assertEqual(len(set(_ffi.EXPECTED_SYMBOLS)), 22,
                         '符号名不得重复')
        lib = _ffi._lib()
        for name in _ffi.EXPECTED_SYMBOLS:
            self.assertTrue(hasattr(lib, name), 'CDLL 缺少符号 %s' % name)

    def test_exc_by_code_covers_1_to_11_all_micerror(self):
        for code in range(1, 12):
            self.assertIn(code, _ffi._EXC_BY_CODE,
                          '错误码 %d 未映射异常类' % code)
            cls = _ffi._EXC_BY_CODE[code]
            self.assertTrue(isinstance(cls, type),
                            '错误码 %d 的映射不是类：%r' % (code, cls))
            self.assertTrue(issubclass(cls, MicError),
                            '错误码 %d 的映射 %r 不是 MicError 子类'
                            % (code, cls))

    def test_mic_ffi_lib_nonexistent_path_raises_with_build_hint(self):
        missing = os.path.join(tempfile.gettempdir(),
                               'micpython_test_no_such_dir', 'milp_ffi.dll')
        self.assertFalse(os.path.exists(missing),
                         '测试前提：路径必须不存在：%s' % missing)
        os.environ['MIC_FFI_LIB'] = missing
        with self.assertRaises(FfiLoadError) as ctx:
            _ffi.load_library()
        message = str(ctx.exception)
        self.assertIn(missing, message, '错误消息应含被尝试的路径')
        self.assertIn(BUILD_HINT, message, '错误消息应含构建提示')

    def test_mic_ffi_lib_real_dll_loads(self):
        candidates = _ffi._default_candidate_paths()
        real = [p for p in candidates if os.path.isfile(p)]
        self.assertTrue(real, '默认候选路径均不存在：%r' % (candidates,))
        os.environ['MIC_FFI_LIB'] = real[0]
        lib = _ffi._lib()
        self.assertIsInstance(lib, ctypes.CDLL)
        self.assertEqual(lib.milp_ffi_abi_version(), 1)
        # 强制路径分支同样完成 22 符号绑定。
        for name in _ffi.EXPECTED_SYMBOLS:
            self.assertTrue(hasattr(lib, name))


if __name__ == '__main__':
    unittest.main()
