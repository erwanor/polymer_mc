#ifndef ERR_MSGS_H
#define ERR_MSGS_H

enum error_code {
  OUT_OF_RANGE = 0,

  NUM_OF_ERR_CODE,
};

extern const char err_msgs[NUM_OF_ERR_CODE][128];

#endif
