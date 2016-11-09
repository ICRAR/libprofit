/* copyright notice, etc */
#ifndef _EXAMPLE_H_
#define _EXAMPLE_H_

#include <string>
#include <vector>

#include "profit/profile.h"

namespace profit
{

class ExampleProfile : public Profile {

public:
	ExampleProfile(const Model &model, const std::string &name);
	void validate();
	void evaluate(std::vector<double> &image);

protected:
	bool parameter(const std::string &name, double value);
	bool parameter(const std::string &name, unsigned int value);

private:
	double param1;
	double param2;
	unsigned int param3;

};

} /* namespace profit */

#endif
