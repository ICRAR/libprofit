# Generates OpenCL .h files out from the corresponding .cl
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2016
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of libprofit.
#
# libprofit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libprofit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libprofit.  If not, see <http://www.gnu.org/licenses/>.


file(READ ${KRN_FULL_FNAME} KRN_SRC)

set(_contents
"/**
 * C++-compatible OpenCL kernel source code from ${KRN_NAME}
 *
 * THIS FILE HAS BEEN AUTOMATICALLY GENERATED FROM ${KRN_FNAME}
 * DO NOT EDIT
 */

#include <string>

namespace profit {

const std::string ${KRN_VNAME} = R\"===(
${KRN_SRC}
)===\";

} // namespace profit")

file(WRITE ${KRN_HEADER_FNAME} "${_contents}")