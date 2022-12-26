#ifndef MEMORY_H
#define MEMORY_H

#if defined(_MSC_VER)
    #include <malloc.h>
#elif defined(__INTEL_COMPILER)
    #include <malloc.h>
#else
    #include <stdlib.h>
#endif

#define MEM_ALIGN_SIZE 64

#ifdef __cplusplus
extern "C" {
#endif

inline void *aligned_malloc(size_t size, size_t alignment)
{
    void *res = NULL;
#if defined(_MSC_VER)
    res = _aligned_malloc(size, alignment);
#elif defined(__INTEL_COMPILER)
    res = _mm_malloc(size, alignment);
#else

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
    res = aligned_alloc(alignment, size);
#else
    if (posix_memalign(&res, alignment, size)) res = NULL;
#endif

#endif
    return res;
}

inline void *aligned_malloc_64byte(size_t size)
{
    return aligned_malloc(size, 64);
}

inline void aligned_free(void *ptr)
{
#ifdef _MSC_VER
    _aligned_free(ptr);
#elif defined(__INTEL_COMPILER)
    _mm_free(ptr);
#else
    free(ptr);
#endif
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif