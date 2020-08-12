/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef BOGUS_SECOND_ORDER_IMPL_HPP
#define BOGUS_SECOND_ORDER_IMPL_HPP

#include "SOC/HBLaw.impl.hpp"
#include "SOC/SOCLaw.impl.hpp"
#include "SecondOrder.hpp"
#ifndef BOGUS_WITHOUT_EIGEN
#include "../Core/Eigen/EigenLinearSolvers.hpp"
#endif
#include "SOC/FischerBurmeister.impl.hpp"
#include "SOC/LocalHBSolver.impl.hpp"
#include "SOC/LocalSOCSolver.impl.hpp"

#endif
