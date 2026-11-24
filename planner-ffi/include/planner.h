#ifndef MIC_PLANNER_H
#define MIC_PLANNER_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

uint32_t planner_abi_version(void);
/* 0 表示创建失败。不同句柄可并行求解。 */
uint64_t planner_create(void);
void planner_destroy(uint64_t handle);
/* data 指向 len 个可读字节。返回 0 表示保存了结果，1 表示粘滞错误，
   2 表示句柄失效或内部错误。对同一句柄，调用方应将求解、量取长度和复制串行化。 */
int32_t planner_run(uint64_t handle, const uint8_t *data, size_t len);
/* 读取最近一次结果或错误的 UTF-8 字节。无效句柄的长度为 SIZE_MAX。 */
size_t planner_size(uint64_t handle);
/* out 是调用方提供的可写缓冲区，容量至少为 planner_size 返回的字节数。 */
int32_t planner_read(uint64_t handle, uint8_t *out, size_t capacity);

#ifdef __cplusplus
}
#endif
#endif
