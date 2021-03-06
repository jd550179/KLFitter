/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */

#include "KLFitter/ResSingleGaussMET.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::ResSingleGaussMET::ResSingleGaussMET(const char * filename) : KLFitter::ResSingleGaussLinearBase(filename) { }

// ---------------------------------------------------------
KLFitter::ResSingleGaussMET::ResSingleGaussMET(std::vector<double> const& parameters) : KLFitter::ResSingleGaussLinearBase(parameters) { }

// ---------------------------------------------------------
KLFitter::ResSingleGaussMET::~ResSingleGaussMET() = default;

// ---------------------------------------------------------
double KLFitter::ResSingleGaussMET::GetMean(double /*x*/) {
  return 0;
}

// ---------------------------------------------------------
double KLFitter::ResSingleGaussMET::GetSigma(double x) {
  return fParameters[0] + fParameters[1] * x;
}
