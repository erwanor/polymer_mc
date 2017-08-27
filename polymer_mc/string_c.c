#include "string_c.h"

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

// NOTE: size <- include '\0'
//       len  <- do not include '\0'
struct string_t {
  char* data;
  size_t capacity;
  size_t size;
};

static string* allocate(void);
static size_t get_ptr_from_org(const int addr, const size_t len);
static void expand_capacity_if_needed(string* self, const size_t new_size);

#undef INIT_MAX_DATA_SIZE
#define INIT_MAX_DATA_SIZE 4

static string* allocate(void) {
  string* self   = (string*) xmalloc(sizeof(string));
  self->data     = (char*) xmalloc(sizeof(char) * INIT_MAX_DATA_SIZE);
  self->capacity = INIT_MAX_DATA_SIZE;
  self->size     = 1;
  self->data[0] = '\0';
  return self;
}

string* new_string(void) {
  string* self = allocate();
  memset(self->data, '\0', sizeof(char) * INIT_MAX_DATA_SIZE);
  return self;
}

string* new_string_from_char(const char* src) {
  return new_string_from_char_n(src, strlen(src));
}

string* new_string_from_char_n(const char* src, const size_t len) {
  string* self = allocate();
  append_char_n(self, src, len);
  return self;
}

string* new_string_from_string(const string* src) {
  return new_string_from_char_n(src->data, strlen_string(src));
}

void delete_string(string* self) {
  xfree(self->data);
  self->data     = NULL;
  self->capacity = 0;
  self->size     = 0;
  xfree(self);
}

void reserve_string(string* self, const size_t new_capacity) {
  self->capacity = new_capacity;
  char* buffer_new = (char*) xmalloc(self->capacity * sizeof(char));
#if _MSC_VER
  strncpy_s(buffer_new, self->size, self->data, self->size - 1);
#else
  strncpy(buffer_new, self->data, self->size);
#endif
  xfree(self->data);
  self->data = buffer_new;
}

void resize_string(string* self, const size_t new_size) {
  resize_noinit_string(self, new_size);
  memset(self->data, '\0', self->size);
}

void resize_noinit_string(string* self, const size_t new_size) {
  reserve_string(self, new_size);
  self->size = new_size;
}

void clear_string(string* self) {
  self->size = 1; // only '\0'
}

void append_string(string* self, const string* str) {
  append_char_n(self, str->data, strlen_string(str));
}

void append_char(string* self, const char* src) {
  append_char_n(self, src, strlen(src));
}

void append_char_n(string* self, const char* src, const size_t len) {
  if (len == 0) return;
  const size_t new_str_size = self->size + len;
  expand_capacity_if_needed(self, new_str_size);
#if _MSC_VER
  strncpy_s(self->data + strlen_string(self), len + 1, src, len);
#else
  strncpy(self->data + strlen_string(self), src, len);
#endif
  self->size = new_str_size;
  self->data[self->size - 1] = '\0';
}

string* string_to_string(const string* self) {
  return new_string_from_string(self);
}

char* string_to_char(const string* self) {
  return self->data;
}

size_t capacity_string(const string* self) {
  return self->capacity;
}

size_t size_string(const string* self) {
  return self->size;
}

size_t strlen_string(const string* self) {
  return self->size - 1;
}

static size_t get_ptr_from_org(const int addr, const size_t len) {
  if (addr < 0) {
    return (size_t) (addr + (int) len);
  } else {
    return (size_t) addr;
  }
}

static void expand_capacity_if_needed(string* self, const size_t new_size) {
  if (new_size > self->capacity) {
    reserve_string(self, new_size * 2);
  }
}

string* slice_string(string* self, const int beg, const int end) {
  const size_t len = strlen_string(self);
  const size_t beg_from_org = get_ptr_from_org(beg, len);
  const size_t end_from_org = get_ptr_from_org(end, len);

  if (beg_from_org >= end_from_org) {
    DEBUG_PRINT("Cannot slice at %s:%s", __FILE__, __LINE__);
    return NULL;
  }
  if (end_from_org >= len) {
    fprintf(stderr, "%s: at %s %d\n", err_msgs[OUT_OF_RANGE], __FILE__, (int32_t)__LINE__);
    return NULL;
  }

  return new_string_from_char_n(self->data + beg_from_org, end_from_org - beg_from_org);
}

bool eq_string(const string* str0,
               const string* str1) {
  if (size_string(str0) != size_string(str1)) return false;
  const size_t len_str = strlen_string(str0);
  return (strncmp(string_to_char(str0), string_to_char(str1), len_str) == 0);
}

DEFINE_VECTOR_METHODS(ptr_string)

vector_ptr_string* split_string(const string* self, const char* delimiters) {
  string* buffer = new_string_from_string(self);
  vector_ptr_string* vptr_string = vector_ptr_string_new();
  char* splitted;
#if _MSC_VER
  char* cxt;
  splitted = strtok_s(string_to_char(buffer), delimiters, &cxt);
#else
  splitted = strtok(string_to_char(buffer), delimiters);
#endif
  ptr_string splitted_string = new_string_from_char(splitted);
  vector_ptr_string_push_back(vptr_string, splitted_string);
#if _MSC_VER
  while ((splitted = strtok_s(NULL, delimiters, &cxt)) != NULL) {
#else
  while ((splitted = strtok(NULL, delimiters)) != NULL) {
#endif
    splitted_string = new_string_from_char(splitted);
    vector_ptr_string_push_back(vptr_string, splitted_string);
  }
  delete_string(buffer);
  return vptr_string;
}

void delete_splitted_strings(vector_ptr_string* vptr_string) {
  const size_t size = vector_ptr_string_size(vptr_string);
  for (size_t i = 0; i < size; i++) {
    delete_string(vector_ptr_string_at_nocheck(vptr_string, i));
  }
  vector_ptr_string_delete(vptr_string);
}

#define DEFINE_STR_TO_DTYPE(Dtype, Converter)                           \
  DECL_STR_TO_DTYPE(Dtype) {                                            \
    char *endp = NULL, *endp_expected = NULL;                           \
    char* ptr_str_data = string_to_char(str);                           \
    const size_t size_str = size_string(str);                           \
    endp_expected = ptr_str_data + size_str - 1;                        \
    const Dtype ret = (Dtype)Converter;                                 \
    if (endp != endp_expected) {                                        \
      fprintf(stderr, "Error occurs during %s at %s:%d\n", STR(Converter), __FILE__, __LINE__); \
      exit(EXIT_FAILURE);                                               \
    }                                                                   \
    return ret;                                                         \
  }
DEFINE_STR_TO_DTYPE(double, strtod(ptr_str_data, &endp))
DEFINE_STR_TO_DTYPE(float, strtof(ptr_str_data, &endp))
DEFINE_STR_TO_DTYPE(int32_t, strtol(ptr_str_data, &endp, 0))
DEFINE_STR_TO_DTYPE(int, strtol(ptr_str_data, &endp, 0))
DEFINE_STR_TO_DTYPE(uint32_t, strtoul(ptr_str_data, &endp, 0))
DEFINE_STR_TO_DTYPE(int64_t, strtoll(ptr_str_data, &endp, 0))
DEFINE_STR_TO_DTYPE(uint64_t, strtoull(ptr_str_data, &endp, 0))
