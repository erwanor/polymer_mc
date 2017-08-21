#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>
#include <stdio.h>

void* xmalloc(size_t size);
void xfree(void* ptr);
FILE* xfopen(const char* __restrict filename,
            const char* __restrict mode);
void xfclose(FILE* fp);

#if _MSC_VER
#define xfscanf(fin, format, ...) fscanf_s(fin, format, __VA_ARGS__);
#else
#define xfscanf(fin, format, ...) fscanf(fin, format, __VA_ARGS__);
#endif

#define UNUSED_PARAMETER(val) (void)(val);

#endif
