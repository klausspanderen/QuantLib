/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2009 Andreas Gaida
 Copyright (C) 2009 Ralph Schreyer
 Copyright (C) 2009, 2017, 2021 Klaus Spanderen

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

#include <ql/math/matrixutilities/psor.hpp>
#include <ql/math/matrixutilities/gmres.hpp>
#include <ql/math/matrixutilities/bicgstab.hpp>
#include <ql/methods/finitedifferences/schemes/impliciteulerscheme.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmamericanstepcondition.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmbermudanstepcondition.hpp>
#include <ql/functional.hpp>

namespace QuantLib {

    namespace {
        ext::shared_ptr<FdmStepConditionComposite> filterExerciseConditions(
            const ext::shared_ptr<FdmStepConditionComposite>& stepConditions) {

            typedef FdmStepConditionComposite::Conditions Conditions;

            const Conditions conditions =
                (stepConditions) != 0 ? stepConditions->conditions() : Conditions();

            Conditions exerciseConditions;
            for (Conditions::const_iterator iter = conditions.begin();
                 iter != conditions.end(); ++iter) {

                const ext::shared_ptr<FdmStepConditionComposite> scc =
                    ext::dynamic_pointer_cast<FdmStepConditionComposite>(*iter);
                if (scc != 0) {
                    const Conditions c =
                        filterExerciseConditions(scc)->conditions();
                    exerciseConditions.insert(
                        exerciseConditions.end(), c.begin(), c.end());
                }

                if ((ext::dynamic_pointer_cast<FdmAmericanStepCondition>(*iter) != 0) ||
                    (ext::dynamic_pointer_cast<FdmBermudanStepCondition>(*iter) != 0))
                    exerciseConditions.push_back(*iter);
            }

            return ext::make_shared<FdmStepConditionComposite>(
                std::list<std::vector<Time> >(),
                exerciseConditions);
        }

        Disposable<Array> applyTo(
            const ext::shared_ptr<FdmStepConditionComposite>& conditions,
            const Array& u, Time t) {

            Array v(u);
            conditions->applyTo(v, t);

            return v;
        }
    }

    ImplicitEulerScheme::ImplicitEulerScheme(
        const ext::shared_ptr<FdmLinearOpComposite>& map,
        const bc_set& bcSet,
        Real relTol,
        SolverType solverType,
        const ext::shared_ptr<FdmStepConditionComposite>& stepConditions)
        : dt_(Null<Real>()),
          iterations_(ext::make_shared<Size>(0U)),
          relTol_(relTol),
          map_(map),
          bcSet_(bcSet),
          solverType_(solverType),
          stepConditions_(stepConditions) {}

    Disposable<Array> ImplicitEulerScheme::apply(const Array& r, Real theta) const {
        return r - (theta*dt_)*map_->apply(r);
    }

#if !defined(QL_NO_UBLAS_SUPPORT)
    Disposable<SparseMatrix>
    ImplicitEulerScheme::applyToSparseMatrix(Real theta) const {
        SparseMatrix m((-theta*dt_)*map_->toMatrix());
        for (Size i=0; i < m.size1(); ++i)
            m(i,i) += 1.0;

        return m;
    }
#endif

    void ImplicitEulerScheme::step(array_type& a, Time t) {
        step(a, t, 1.0);
    }

    void ImplicitEulerScheme::step(array_type& a, Time t, Real theta) {
        using namespace ext::placeholders;
        QL_REQUIRE(t-dt_ > -1e-8, "a step towards negative time given");
        map_->setTime(std::max(0.0, t-dt_), t);
        bcSet_.setTime(std::max(0.0, t-dt_));

        bcSet_.applyBeforeSolving(*map_, a);


        if (solverType_ == PSOR) {
#if !defined(QL_NO_UBLAS_SUPPORT)
            const ext::shared_ptr<FdmStepConditionComposite> exerciseConditions(
                filterExerciseConditions(stepConditions_));

            const PSOR::Projection proj(
                ext::bind(&applyTo, exerciseConditions, _1, t));

            const PSORResult result =
                QuantLib::PSOR(applyToSparseMatrix(theta),
                    1.5, 100000, relTol_, proj).solve(a, a);

            (*iterations_) += result.iterations;
            a = result.x;
#else
            QL_FAIL("PSOR is not supported due to missing boost ublas");
#endif
        }
        else if (map_->size() == 1) {
            a = map_->solve_splitting(0, a, -theta*dt_);
        }
        else {
            const ext::function<Disposable<Array>(const Array&)> applyF(
                ext::bind(&ImplicitEulerScheme::apply, this, _1, theta));

            const ext::function<Disposable<Array>(const Array&)>
                preconditioner(ext::bind(
                    &FdmLinearOpComposite::preconditioner, map_, _1, -theta*dt_));

            if (solverType_ == BiCGstab) {
                const BiCGStabResult result =
                    QuantLib::BiCGstab(applyF, std::max(Size(10), a.size()),
                        relTol_, preconditioner).solve(a, a);

                (*iterations_) += result.iterations;
                a = result.x;
            }
            else if (solverType_ == GMRES) {
                const GMRESResult result =
                    QuantLib::GMRES(applyF, std::max(Size(10), a.size() / 10U), relTol_,
                                    preconditioner)
                        .solve(a, a);

                (*iterations_) += result.errors.size();
                a = result.x;
            }
            else
                QL_FAIL("unknown/illegal solver type");
        }
        bcSet_.applyAfterSolving(a);
    }

    void ImplicitEulerScheme::setStep(Time dt) {
        dt_=dt;
    }

    Size ImplicitEulerScheme::numberOfIterations() const {
        return *iterations_;
    }
}
