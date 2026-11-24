#include <jni.h>
#include <limits.h>
#include <stdlib.h>
#include "planner.h"

static void throw_error(JNIEnv *env, const char *type, const char *message) {
    jclass cls = (*env)->FindClass(env, type);
    if (cls != NULL) (*env)->ThrowNew(env, cls, message);
}

/* Rust 复制业务模型后自行搜索，JNI 不长期占用或固定 JVM 的堆内存。
   每次求解临时登记一个句柄，错误、异常和正常完成共用同一释放路径。 */
JNIEXPORT jbyteArray JNICALL Java_org_micplanner_Native_run(
        JNIEnv *env, jclass klass, jbyteArray model) {
    (void)klass;
    if (planner_abi_version() != 1) {
        throw_error(env, "java/lang/UnsatisfiedLinkError", "planner ABI version 1 required");
        return NULL;
    }
    if (model == NULL) {
        throw_error(env, "java/lang/NullPointerException", "model");
        return NULL;
    }
    jsize len = (*env)->GetArrayLength(env, model);
    uint8_t *input = malloc(len ? (size_t)len : 1);
    uint8_t *output = NULL;
    uint64_t handle = 0;
    jbyteArray result = NULL;
    if (input == NULL) goto no_memory;
    if (len) (*env)->GetByteArrayRegion(env, model, 0, len, (jbyte *)input);
    if ((*env)->ExceptionCheck(env)) goto done;
    handle = planner_create();
    if (!handle) goto no_memory;
    int32_t rc = planner_run(handle, input, (size_t)len);
    size_t size = planner_size(handle);
    if (size >= INT_MAX) {
        throw_error(env, "java/lang/IllegalStateException", "planner result too large or handle invalid");
        goto done;
    }
    output = malloc(size + 1);
    if (output == NULL) goto no_memory;
    output[0] = rc != 0;
    if (planner_read(handle, output + 1, size) != 0) {
        throw_error(env, "java/lang/IllegalStateException", "cannot read planner result");
        goto done;
    }
    result = (*env)->NewByteArray(env, (jsize)size + 1);
    if (result != NULL) (*env)->SetByteArrayRegion(env, result, 0, (jsize)size + 1, (jbyte *)output);
    goto done;
no_memory:
    throw_error(env, "java/lang/OutOfMemoryError", "planner native allocation");
done:
    if (handle) planner_destroy(handle);
    free(output);
    free(input);
    return result;
}
