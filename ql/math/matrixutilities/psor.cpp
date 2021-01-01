/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2021 Klaus Spanderen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file psor.cpp
    \brief projected successive over-relaxation
*/

#include <ql/math/matrixutilities/psor.hpp>

#if !defined(QL_NO_UBLAS_SUPPORT)

namespace QuantLib {

	PSOR::PSOR(const SparseMatrix& A,
			   Real omega, Size maxIter, Real relTol,
			   const PSOR::Projection& proj)
	: A_(A),
	  omega_(omega),
	  maxIter_(maxIter),
	  relTol_(relTol),
	  proj_(proj) {	}

	PSORResult PSOR::solve(const Array& b, const Array& x0) const {
        const Real bnorm2 = Norm2(b);
        if (bnorm2 == 0.0) {
        	PSORResult result = { 0, 0.0, b};
            return result;
        }

		const Size n = b.size();
        Array x = (!x0.empty()) ? x0 : Array(n, 0.0);
        Array x_old = x;


        const Array y = (proj_)
        	? proj_(Array(n, 0.0))
        	: Array(n, std::numeric_limits<Real>::lowest());

        Real error = Norm2(b - prod(A_, x))/bnorm2;

        Size i;
        for (i=0; i < maxIter_ && error >= relTol_; ++i) {
        	for (Size k=0; k < n; ++k) {
        		const Real a = A_(k,k);
        		const Real s = prod_i(A_, x, k);
        		x[k] = std::max(y[k],
        				(1-omega_)*x[k] + omega_/a*(b[k] - s + a*x[k]));
        	}
        	error = Norm2(x_old - x)/bnorm2;
        	x_old = x;
        }

        QL_REQUIRE(i < maxIter_, "max number of iterations exceeded");
        QL_REQUIRE(error < relTol_, "could not converge");

        PSORResult result = {i, error, x};
		return result;
	}
}

#endif
