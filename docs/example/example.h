/* copyright notice, etc */
#ifndef _EXAMPLE_H_
#define _EXAMPLE_H_

#include "profit.h"

typedef struct _profit_example_profile {
	profit_profile profile;
	double param1;
	double param2;
	int param3;
} profit_example_profile;

profit_profile *profit_create_example(void);

#endif
