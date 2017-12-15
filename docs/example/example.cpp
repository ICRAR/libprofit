/* copyright statement, etc */

#include <cmath>
#include "example.h"

#include "profit/exceptions.h"
#include "profit/model.h"

namespace profit {

ExampleProfile::ExampleProfile(const Model &model, const std::string &name) :
    Profile(model, name),
	 param1(1.),
	 param2(2.),
	 param3(3)
{
	// no-op
}

bool ExampleProfile::parameter(const std::string &name, double value) {

	if( Profile::parameter(name, value) ) {
		return true;
	}

	if( name == "param1" )      { param1 = value; }
	else if( name == "param2" ) { param2 = value; }
	else {
		return false;
	}

	return true;
}

bool ExampleProfile::parameter(const std::string &name, unsigned int value) {
	if( Profile::parameter(name, value) ) {
		return true;
	}

	if( name == "param3" ) { param3 = value; }
	else {
		return false;
	}

	return true;
}

void ExampleProfile::validate() {

	if ( this->param1 < 0 ) {
		throw invalid_parameter("param1 is negative");
	}
	if ( this->param1 < 0 ) {
		throw invalid_parameter("param2 is negative");
	}
	if ( this->param3 < 0 ) {
		throw invalid_parameter("param3 is negative");
	}

	/*
	if ( this->model->width < 20 || this->model->height < 20 ) {
		throw invalid_parameter("can't apply example profile to images less than 20x20");
	}
	*/
}

void ExampleProfile::evaluate(Image &image, const Mask &mask, const PixelScale &scale, double magzero) {

	double x, y;
	unsigned int i, j;
	auto width = image.getWidth();
	double half_xbin = scale.first/2.;
	double half_ybin = scale.second/2.;

	x = 0;
	for (i=0; i < width; i++) {
		x += half_xbin;

		y = 0;
		for (j=0; j < image.getHeight(); j++) {
			y += half_ybin;

			if ( not mask or mask[i + j * width] ) {
				double val = std::abs( (this->param1 - this->param2) * this->param3 * (x - y) );
				image[i + j * width] = val;
			}

			y += half_ybin;
		}
		x += half_xbin;
	}
}

} /* namespace profit */
