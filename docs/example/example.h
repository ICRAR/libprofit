/* copyright notice, etc */
#ifndef _EXAMPLE_H_
#define _EXAMPLE_H_

#include <vector>

#include "profit.h"

namespace profit
{

class ExampleProfile : public Profile {

public:

	ExampleProfile();
	void validate();
	void evaluate(std::vector<double> &image);

	double param1;
	double param2;
	int param3;

};

} /* namespace profit */

#endif
