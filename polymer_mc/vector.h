#ifndef VECTOR_H
#define VECTOR_H

#include <stdbool.h>
#include <string.h>
#include "utils.h"

#define VECTOR_TYPE(T) CONCAT(vector_, T)
#define VECTOR_REF(T) CONCAT(vectorref_, T)
#define VECTOR_REF_CONST(T) CONCAT(vectorrefc_, T)
#define VECTOR_TYPE_DUMMY(T) CONCAT(vector_t_, T)

#define VECTOR_METHOD(T, f) CONCAT(CONCAT(VECTOR_TYPE(T), _), f)

#undef INIT_MAX_DATA_SIZE
#define INIT_MAX_DATA_SIZE 128

#define DECL_VECTOR_TYPE(T)                           \
  struct VECTOR_TYPE_DUMMY(T);                        \
  typedef struct VECTOR_TYPE_DUMMY(T) VECTOR_TYPE(T); \
  typedef struct VECTOR_TYPE_DUMMY(T)* VECTOR_REF(T); \
  typedef const struct VECTOR_TYPE_DUMMY(T)* VECTOR_REF_CONST(T);

#define DEFINE_VECTOR_TYPE(T)                   \
  DECL_VECTOR_TYPE(T)                           \
  struct VECTOR_TYPE_DUMMY(T) {                 \
    T* data;                                    \
    size_t capacity;                            \
    size_t size;                                \
  };

// public members
#define DECL_VECTOR_NEW(T) VECTOR_REF(T) VECTOR_METHOD(T, new)(void)
#define DEFINE_VECTOR_NEW(T)                                            \
  DECL_VECTOR_NEW(T) {                                                  \
    VECTOR_REF(T) self = (VECTOR_REF(T)) xmalloc(sizeof(VECTOR_TYPE(T))); \
    self->capacity = INIT_MAX_DATA_SIZE;                                \
    self->data = (T*) xmalloc(self->capacity * sizeof(T));              \
    self->size = 0;                                                     \
    return self;                                                        \
  }

#define DECL_VECTOR_DELETE(T) void VECTOR_METHOD(T, delete)(VECTOR_REF(T) self)
#define DEFINE_VECTOR_DELETE(T)                 \
  DECL_VECTOR_DELETE(T) {                       \
    xfree(self->data);                          \
    self->data = NULL;                          \
    self->size = 0;                             \
    self->capacity = 0;                         \
    xfree(self);                                \
  }

#define DECL_VECTOR_PUSH_BACK(T) void VECTOR_METHOD(T, push_back)(VECTOR_REF(T) self, const T in)
#define DEFINE_VECTOR_PUSH_BACK(T)                              \
  DECL_VECTOR_PUSH_BACK(T) {                                    \
    if ((self->size + 1) >= self->capacity) {                   \
      VECTOR_METHOD(T, reserve)(self, self->capacity * 2);      \
    }                                                           \
    self->data[self->size] = in;                                \
    self->size++;                                               \
  }

#define DECL_VECTOR_POP_BACK(T) void VECTOR_METHOD(T, pop_back)(VECTOR_REF(T) self)
#define DEFINE_VECTOR_POP_BACK(T)              \
  DECL_VECTOR_POP_BACK(T) {                    \
    if (!VECTOR_METHOD(T, empty)(self)) {      \
      self->size--;                            \
    }                                          \
  }

#define DECL_VECTOR_CLEAR(T) void VECTOR_METHOD(T, clear)(VECTOR_REF(T) self)
#define DEFINE_VECTOR_CLEAR(T)               \
  DECL_VECTOR_CLEAR(T) {                     \
    self->size = 0;                          \
  }

#define DECL_VECTOR_AT(T) T VECTOR_METHOD(T, at)(VECTOR_REF_CONST(T) self, const size_t idx)
#define DEFINE_VECTOR_AT(T)                                             \
  DECL_VECTOR_AT(T) {                                                   \
    RANGE_CHECK(self, idx);                                             \
    return self->data[idx];                                             \
  }

#define DECL_VECTOR_AT_NOCHECK(T) T VECTOR_METHOD(T, at_nocheck)(VECTOR_REF_CONST(T) self, const size_t idx)
#define DEFINE_VECTOR_AT_NOCHECK(T)                                   \
  DECL_VECTOR_AT_NOCHECK(T) {                                         \
    return self->data[idx];                                           \
  }

#define DECL_VECTOR_ASSIGN(T) void VECTOR_METHOD(T, assign)(VECTOR_REF(T) self, const size_t idx, const T val)
#define DEFINE_VECTOR_ASSIGN(T)                 \
  DECL_VECTOR_ASSIGN(T) {                       \
    RANGE_CHECK(self, idx);                     \
    self->data[idx] = val;                      \
  }

#define DECL_VECTOR_ASSIGN_NOCHECK(T) void VECTOR_METHOD(T, assign_nocheck)(VECTOR_REF(T) self, const size_t idx, const T val)
#define DEFINE_VECTOR_ASSIGN_NOCHECK(T)                 \
  DECL_VECTOR_ASSIGN_NOCHECK(T) {                       \
    self->data[idx] = val;                              \
  }

#define DECL_VECTOR_FRONT(T) T VECTOR_METHOD(T, front)(VECTOR_REF_CONST(T) self)
#define DEFINE_VECTOR_FRONT(T)                          \
  DECL_VECTOR_FRONT(T) {                                \
    return VECTOR_METHOD(T, at_nocheck)(self, 0);               \
  }

#define DECL_VECTOR_BACK(T) T VECTOR_METHOD(T, back)(VECTOR_REF_CONST(T) self)
#define DEFINE_VECTOR_BACK(T)                           \
  DECL_VECTOR_BACK(T) {                                 \
    return VECTOR_METHOD(T, at_nocheck)(self, self->size - 1);  \
  }

#define DECL_VECTOR_RESERVE(T) void VECTOR_METHOD(T, reserve)(VECTOR_REF(T) self, size_t new_capacity)
#define DEFINE_VECTOR_RESERVE(T)                              \
  DECL_VECTOR_RESERVE(T) {                                    \
    self->capacity = new_capacity;                            \
    T* buffer_new = (T*) xmalloc(self->capacity * sizeof(T)); \
    memcpy(buffer_new, self->data, self->size * sizeof(T));   \
    xfree(self->data);                                        \
    self->data = buffer_new;                                  \
  }

#define DECL_VECTOR_RESIZE(T) void VECTOR_METHOD(T, resize)(VECTOR_REF(T) self, size_t new_size, T val)
#define DEFINE_VECTOR_RESIZE(T)                                   \
  DECL_VECTOR_RESIZE(T) {                                         \
    if (new_size > self->capacity) {                              \
      VECTOR_METHOD(T, reserve)(self, new_size);                  \
    }                                                             \
    self->size = new_size;                                        \
    for (size_t i = 0; i < self->size; i++) self->data[i] = val;  \
  }

#define DECL_VECTOR_SIZE(T) size_t VECTOR_METHOD(T, size)(VECTOR_REF_CONST(T) self)
#define DEFINE_VECTOR_SIZE(T)                 \
  DECL_VECTOR_SIZE(T) {                       \
    return self->size;                        \
  }

#define DECL_VECTOR_CAPACITY(T) size_t VECTOR_METHOD(T, capacity)(VECTOR_REF_CONST(T) self)
#define DEFINE_VECTOR_CAPACITY(T)                   \
  DECL_VECTOR_CAPACITY(T) {                         \
    return self->capacity;                          \
  }

#define DECL_VECTOR_EMPTY(T) bool VECTOR_METHOD(T, empty)(VECTOR_REF_CONST(T) self)
#define DEFINE_VECTOR_EMPTY(T)               \
  DECL_VECTOR_EMPTY(T) {                     \
    return (bool)(0 == self->size);          \
  }

#define DEFINE_VECTOR_METHODS(T)                \
  DEFINE_VECTOR_NEW(T)                          \
  DEFINE_VECTOR_DELETE(T)                       \
  DEFINE_VECTOR_PUSH_BACK(T)                    \
  DEFINE_VECTOR_POP_BACK(T)                     \
  DEFINE_VECTOR_CLEAR(T)                        \
  DEFINE_VECTOR_AT(T)                           \
  DEFINE_VECTOR_AT_NOCHECK(T)                   \
  DEFINE_VECTOR_ASSIGN(T)                       \
  DEFINE_VECTOR_ASSIGN_NOCHECK(T)               \
  DEFINE_VECTOR_FRONT(T)                        \
  DEFINE_VECTOR_BACK(T)                         \
  DEFINE_VECTOR_RESERVE(T)                      \
  DEFINE_VECTOR_RESIZE(T)                       \
  DEFINE_VECTOR_SIZE(T)                         \
  DEFINE_VECTOR_CAPACITY(T)                     \
  DEFINE_VECTOR_EMPTY(T)

#define DECL_VECTOR_METHODS(T)                  \
  DECL_VECTOR_NEW(T);                           \
  DECL_VECTOR_DELETE(T);                        \
  DECL_VECTOR_PUSH_BACK(T);                     \
  DECL_VECTOR_POP_BACK(T);                      \
  DECL_VECTOR_CLEAR(T);                         \
  DECL_VECTOR_AT(T);                            \
  DECL_VECTOR_AT_NOCHECK(T);                    \
  DECL_VECTOR_ASSIGN(T);                        \
  DECL_VECTOR_ASSIGN_NOCHECK(T);                \
  DECL_VECTOR_FRONT(T);                         \
  DECL_VECTOR_BACK(T);                          \
  DECL_VECTOR_RESERVE(T);                       \
  DECL_VECTOR_RESIZE(T);                        \
  DECL_VECTOR_SIZE(T);                          \
  DECL_VECTOR_CAPACITY(T);                      \
  DECL_VECTOR_EMPTY(T);

#endif
