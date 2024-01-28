#include "file_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

FILE *xfopen(const char *__restrict fname,
             const char *__restrict mode)
{
  FILE *fp;
  printf("fname: %s\n", fname);
  if ((fp = fopen(fname, mode)) == NULL)
  {
    fprintf(stderr, "Error occurs\n");
    exit(1);
  }
  else
  {
    return fp;
  }
}

void xfclose(FILE *fp)
{
  fclose(fp);
}

size_t get_file_size(FILE *fp)
{
  struct stat stbuf;
  const int fd = fileno(fp);
  if (fstat(fd, &stbuf) == -1)
  {
    perror("Error occurs at get_file_size");
    exit(1);
  }
  return stbuf.st_size;
}

vector_ptr_string *read_lines(FILE *fp)
{
  string *buffer = new_string();
  const size_t file_size = get_file_size(fp) + 1; // include '\0'
  resize_noinit_string(buffer, file_size);
  const size_t num_read = fread(string_to_char(buffer), sizeof(char), file_size, fp);
  UNUSED_PARAMETER(num_read);
  vector_ptr_string *ret = split_string(buffer, "\n");
  delete_string(buffer);
  return ret;
}
