/**
 * Header file for utility routines
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
 *
 * This file is part of libprofit.
 *
 * libprofit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libprofit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PROFIT_UTILS_H
#define PROFIT_UTILS_H

#include <string>
#include <vector>

#include "profit/common.h"

namespace profit
{

/**
 * Checks whether values `x` and `y` are equals,
 * with a difference of almost `e`.
 * @param x The first value to compare
 * @param y The second value to compare
 * @param e The maximum allowed difference
 * @return Whether the values are almost equals
 */
PROFIT_API bool almost_equals(double x, double y, double e = 1e-10);


/**
 * Computes the quantile of the gamma distribution for ``p`` and ``shape``
 */
PROFIT_API double qgamma(double p, double shape);

/**
 * Computes the probability of the gamma distribution for ``q`` and ``shape``
 */
PROFIT_API double pgamma(double q, double shape);

/**
 * Computes the gamma function for ``x``.
 *
 * @param x The parameter of the gamma function
 * @returns The value of `gamma(x)`. If `x` is a negative integer
 * `NaN` is returned. If underflow occurs, `0` is returned.
 * If x is positive and overflow occurs, `+Inf` is returned.
 */
PROFIT_API double gammafn(double x);

/**
 * Computes the beta function for ``a`` and ``b``.
 * This function follows R's beta function semantics;
 * that is, the result is undefined if either is negative,
 * and +Inf if either parameter is zero.
 *
 * @param a The first parameter of the beta function
 * @param b The second parameter of the beta function
 * @returns The value of `beta(a,b)`.
 * If either is negative `NaN` is returned. If either is 0,
 * `+Inf` is returned.
 */
PROFIT_API double beta(double a, double b);

/**
 * A function that can be integrated
 *
 * @param x The domain value used to evaluate the function
 * @param params Additional parameters used to calculate the function
 * @return The value of the integration function at `x`.
 */
typedef double (*integration_func_t)(double x, void *params);

/**
 * Integrates the function `f` on the semi-infinite interval (a, +Inf) using the
 * QAG algorithm (originally from QUADPACK).
 *
 * @param f The function to integrate
 * @param a The beginning of the integration interval
 * @param params A void pointer to any extra data needed by `f`
 * @return The integration result
 */
PROFIT_API double integrate_qagi(integration_func_t f, double a, void *params);

/**
 * Integrates the function `f` on the defined interval (a, b) using the
 * QAG algorithm (originally from QUADPACK).
 *
 * @param f The function to integrate
 * @param a The beginning of the integration interval
 * @param b The end of the integration interval
 * @param params A void pointer to any extra data needed by `f`
 * @return The integration result
 */
PROFIT_API double integrate_qags(integration_func_t f, double a, double b, void *params);

/**
 * Returns whether the named directory exists or not
 * @param dname The directory name
 * @return Whether the directory with the given name exists or not
 */
PROFIT_API bool dir_exists(const std::string &dname);

/**
 * Returns whether the named file exists or not
 * @param fname The file name
 * @return Whether the file with the given name exists or not
 */
PROFIT_API bool file_exists(const std::string &fname);

/**
 * Creates the hierarchy of directories given in @p parts at @p at
 * @param at The directory where the hierarchy of new directories will be created
 * @param parts The hierarchy of directories to create. The first element will
 * be created directly under @p at, the second element directly under the first,
 * and so on
 * @return The path to the last created directory
 */
PROFIT_API std::string create_dirs(const std::string &at, const std::vector<std::string> &parts);

/**
 * Recursively remove the given path
 * @param path The path to remove. If it doesn't exist, an error is issued
 */
PROFIT_API void recursive_remove(const std::string &path);

/**
 * Returns the name of the directory where libprofit stores its internal data.
 * If the directory doesn't exist it is created first.
 *
 * @return The name of the libprofit home directory.
 */
PROFIT_API std::string get_profit_home();

/**
 * Set the environment variable @a name to @a value. If a variable with that
 * name exists its value is replaced. If @a value is empty the variable is
 * removed from the environment
 *
 * @param name The environment variable's name
 * @param value The new value for the environment variable. An empty string
 * causes the variable to be removed from the environment.
 */
PROFIT_API void setenv(const std::string &name, const std::string &value);

/**
 * Breaks down @p s into substrings delimited by @p delims
 * @param s The string to split
 * @param delims The individual delimiters
 * @return The individual substrings after splitting
 */
PROFIT_API std::vector<std::string> split(const std::string &s, const std::string &delims);

/**
 * Trims @p s on both ends
 * @param s The string to trim
 * @return A reference to @p s
 */
PROFIT_API std::string &trim(std::string &s);

/**
 * Trims @p s on both ends
 * @param s The string to trim
 * @return The new trimmed string
 */
PROFIT_API std::string trim(const std::string &s);

/**
 * Like std::stoul, but returns an unsigned int
 * @param s The string to convert
 * @return The unsigned int
 */
inline unsigned int stoui(const std::string &s)
{
	return static_cast<unsigned int>(std::stoul(s));
}

/**
 * Returns the ceiling of the division between @p x and @p y
 * @param x The dividend
 * @param y The divisor
 * @return The ceiling of divind @p x by @p y
 */
inline unsigned int ceil_div(unsigned int x, unsigned int y)
{
	return (x + y - 1) / y;
}

} /* namespace profit */

#endif /* PROFIT_UTILS_H */
