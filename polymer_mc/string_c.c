#include "string_c.h"

#include <string.h>
#include <stdio.h>

#include "utils.h"

// NOTE: size <- include '\0'
//       len  <- do not include '\0'
struct string_t {
	char* data;
	size_t capacity;
	size_t size;
};

static size_t get_ptr_from_org(const int addr, const size_t len);
static void expand_capacity_if_needed(string* self, const size_t new_size);
static string* allocate(void);

#undef INIT_MAX_DATA_SIZE
#define INIT_MAX_DATA_SIZE 16

static string* allocate(void)
{
	string* self = (string*)xmalloc(sizeof(string));
	self->data = (char*)xmalloc(sizeof(char) * INIT_MAX_DATA_SIZE);
	self->capacity = INIT_MAX_DATA_SIZE;
	self->size = 1;
	return self;
}

string* new_string(void)
{
	string* self = allocate();
	memset(self->data, '\0', sizeof(char) * INIT_MAX_DATA_SIZE);
	return self;
}

string* new_string_from_char(const char* src)
{
	return new_string_from_char_n(src, strlen(src));
}

string* new_string_from_char_n(const char* src,
	const size_t len)
{
	string* self = allocate();
	append_char_n(self, src, len);
	return self;
}

string* new_string_from_string(const string* src)
{
	return new_string_from_char_n(src->data, strlen_string(src));
}

void delete_string(string* self)
{
	xfree(self->data);
	self->data = NULL;
	self->capacity = 0;
	self->size = 0;
	xfree(self);
}

void reserve_string(string* self,
	const size_t new_capacity)
{
	self->capacity = new_capacity;
	char* buffer_new = (char*)xmalloc(self->capacity * sizeof(char));
#if _MSC_VER
	strncpy_s(buffer_new, self->size, self->data, self->size - 1);
#else
	strncpy(buffer_new, self->data, self->size);
#endif
	xfree(self->data);
	self->data = buffer_new;
}

void clear_string(string* self)
{
	self->size = 1; // only '\0'
}

void append_string(string* self, const string* str)
{
	append_char_n(self, str->data, strlen_string(str));
}

void append_char(string* self, const char* src)
{
	append_char_n(self, src, strlen(src));
}

void append_char_n(string* self, const char* src, const size_t len)
{
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

char* string_to_char(string* self)
{
	return self->data;
}

size_t capacity_string(const string* self)
{
	return self->capacity;
}

size_t size_string(const string* self)
{
	return self->size;
}

size_t strlen_string(const string* self)
{
	return self->size - 1;
}

static size_t get_ptr_from_org(const int addr,
	const size_t len)
{
	if (addr < 0) {
		return (size_t)(addr + (int)len);
	} else {
		return (size_t)addr;
	}
}

static void expand_capacity_if_needed(string* self,
	const size_t new_size)
{
	if (new_size > self->capacity) {
		reserve_string(self, new_size * 2);
	}
}

string* slice_string(string* self,
	const int beg,
	const int end)
{
	const size_t len = strlen_string(self);
	const size_t beg_from_org = get_ptr_from_org(beg, len);
	const size_t end_from_org = get_ptr_from_org(end, len);

#pragma warning(push)
#pragma warning(disable:4777)
	if (beg_from_org >= end_from_org) {
		fprintf(stderr, "Cannot slice! Error occurs at %s:%d\n", __FILE__, __LINE__);
		return NULL;
	}
	if (end_from_org >= len) {
		fprintf(stderr, "Out of range! Error occurs at %s:%d\n", __FILE__, __LINE__);
		return NULL;
	}
#pragma warning(pop)

	return new_string_from_char_n(self->data + beg_from_org, end_from_org - beg_from_org);
}