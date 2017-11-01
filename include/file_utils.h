#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <stdio.h>

#include "utils.h"
#include "string_c.h"

FILE* xfopen(const char* __restrict fname,
             const char* __restrict mode);
void xfclose(FILE* fp);
size_t get_file_size(FILE* fp);
vector_ptr_string* read_lines(FILE* fp);

#endif
