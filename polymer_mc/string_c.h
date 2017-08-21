#ifndef STRING_C_H
#define STRING_C_H

#include <stddef.h>

struct string_t;
typedef struct string_t string;

string* new_string(void);
string* new_string_from_char(const char* src);
string* new_string_from_char_n(const char* src, const size_t len);
string* new_string_from_string(const string* src);
void delete_string(string* self);

void reserve_string(string* self, const size_t new_capacity);
void clear_string(string* self);

void append_string(string* self, const string* str);
void append_char(string* self, const char* src);
void append_char_n(string* self, const char* src, const size_t len);

size_t capacity_string(const string* self);
size_t size_string(const string* self);
size_t strlen_string(const string* self);
char* string_to_char(string* self);

string* slice_string(string* self, const int beg, const int end);

#endif
