#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void* xmalloc(size_t size)
{
	void* ptr = malloc(size);
	if (!ptr) {
		perror("xmalloc");
		exit(1);
	}
	return ptr;
}

void xfree(void* ptr)
{
	free(ptr);
}

FILE* xfopen(const char* __restrict filename,
            const char* __restrict mode)
{
	FILE* fp;
#if _MSC_VER
	fopen_s(&fp, filename, mode);
#else
	FILE* fp = fopen(filename, mode);
#endif
	if (!fp) {
		perror(filename);
		exit(1);
	}
	return fp;
}

void xfclose(FILE* fp)
{
	fclose(fp);
}