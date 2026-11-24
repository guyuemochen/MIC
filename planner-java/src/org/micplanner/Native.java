package org.micplanner;

import java.io.File;
import java.net.URISyntaxException;

final class Native {
    static {
        String ffi = System.getenv("MIC_PLANNER_LIB");
        if (ffi != null && !ffi.isEmpty()) System.load(new File(ffi).getAbsolutePath());
        String override = System.getenv("MIC_PLANNER_JNI");
        if (override != null && !override.isEmpty()) {
            loadBridge(new File(override), ffi);
        } else {
            File sibling;
            try {
                File location = new File(Native.class.getProtectionDomain().getCodeSource()
                        .getLocation().toURI());
                sibling = new File(location.getParentFile(), System.mapLibraryName("micplanner_jni"));
            } catch (URISyntaxException e) {
                throw new ExceptionInInitializerError(e);
            }
            if (sibling.isFile()) loadBridge(sibling, ffi);
            else System.loadLibrary("micplanner_jni");
        }
    }

    private static void loadBridge(File bridge, String ffiOverride) {
        // Windows 的依赖库搜索路径不一定包含 JAR 所在目录，先显式装入同目录的 Rust 库。
        if (File.separatorChar == '\\' && (ffiOverride == null || ffiOverride.isEmpty())) {
            File ffi = new File(bridge.getAbsoluteFile().getParentFile(), System.mapLibraryName("planner_ffi"));
            if (ffi.isFile()) System.load(ffi.getAbsolutePath());
        }
        System.load(bridge.getAbsolutePath());
    }

    private Native() {}
    static native byte[] run(byte[] model);
}
