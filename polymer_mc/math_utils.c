#include "math_utils.h"

#include <math.h>

bool isSquareNumber(const int32_t a) {
	return ((int32_t)sqrt(a) * (int32_t)sqrt(a) == a);
}
