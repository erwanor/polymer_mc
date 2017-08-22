#include "file_utils.h"

#include <stdio.h>
#include <stdlib.h>

FILE* xfopen(const char* __restrict fname,
             const char* __restrict mode) {
  FILE* fp;
#if _MSC_VER
  errno_t error;
  if ((error = fopen_s(&fp, fname, mode)) != 0) {
#else
  if ((fp = fopen(fname, mode)) == NULL) {
#endif
    fprintf(stderr, "Error occurs\n");
	exit(1);
  } else {
    return fp;
  }
}

void xfclose(FILE* fp) {
  fclose(fp);
}

size_t get_file_size(FILE* fp) {
  fseek(fp, 0, SEEK_END);
  const size_t size = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  return size;
}

vector_ptr_string* read_lines(FILE* fp) {
  string* buffer = new_string();
  const size_t file_size = get_file_size(fp) + 1; // include '\0'
  resize_noinit_string(buffer, file_size);
  fread(string_to_char(buffer), sizeof(char), file_size, fp);
  vector_ptr_string* ret = split_string(buffer, "\n");
  delete_string(buffer);
  return ret;
}
