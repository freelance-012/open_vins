/*
 * OpenVINS 在新 glibc (>=2.34) 上的编译修复。
 *
 * 背景：
 * glibc 的 _DYNAMIC_STACK_SIZE_SOURCE 会把 PTHREAD_STACK_MIN 定义为
 * __sysconf(__SC_THREAD_STACK_MIN_VALUE) 函数调用，导致 boost/thread.hpp 中
 * "#if PTHREAD_STACK_MIN > 0" 预处理失败（函数调用不是数值常量）。
 *
 * 用法：
 * 通过 CMake 的 -include 选项在所有编译单元最早期包含此文件，
 * 强制把 PTHREAD_STACK_MIN 重新定义为数值常量。
 *
 * 参考：/usr/include/x86_64-linux-gnu/bits/pthread_stack_min.h
 */
#pragma once

#ifdef __linux__
#pragma push_macro("PTHREAD_STACK_MIN")
#undef PTHREAD_STACK_MIN
#define PTHREAD_STACK_MIN 16384
#endif
