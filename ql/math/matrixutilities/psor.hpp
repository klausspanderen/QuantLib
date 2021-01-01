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

/*! \file psor.hpp
    \brief projected successive over-relaxation
*/

#ifndef quantlib_psor_hpp
#define quantlib_psor_hpp

#if !defined(QL_NO_UBLAS_SUPPORT)

#include <ql/functional.hpp>
#include <ql/math/matrixutilities/sparsematrix.hpp>

namespace QuantLib {
	struct PSORResult {
		Size iterations;
		Real error;
		Array x;
	};

    class PSOR  {
      public:
        typedef ext::function<Disposable<Array>(const Array&)> Projection;

        PSOR(const SparseMatrix& A,
        	 Real omega, Size maxIter, Real relTol,
             const Projection& proj = Projection());

        PSORResult solve(const Array& b, const Array& x0 = Array()) const;

      protected:
        const SparseMatrix A_;
        const Real omega_;
        const Size maxIter_;
        const Real relTol_;
        const Projection proj_;
    };
}

#endif
#endif
