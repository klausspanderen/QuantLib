/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2025 Klaus Spanderen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file mpanalytichestonengine.hpp
    \brief multi-precision analytic Heston-model engine
*/


#ifndef quantlib_mp_analytic_heston_engine_hpp
#define quantlib_mp_analytic_heston_engine_hpp

#include <ql/math/functional.hpp>
#include <ql/models/equity/hestonmodel.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/genericmodelengine.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/distributions/normal.hpp>

#include <tuple>
#include <complex>
#include <iostream>

using namespace boost::multiprecision;

namespace QuantLib {

    template <class T>
    class MultiPrecisionHestonEngine
        : public GenericModelEngine<HestonModel, VanillaOption::arguments, VanillaOption::results> {

      public:
        enum ContourIntegral {
            Plain, AngledContourShift
        };

        enum ControlVariate {
            Without, BlackScholes, BlackScholes2
        };

        enum Quadrature {
            ExpSinh, TanhSinh, SinhSinh, GaussLaguerre
        };

        MultiPrecisionHestonEngine(
           const ext::shared_ptr<HestonModel> hestonModel,
           const T& precision = T(QL_EPSILON),
           ContourIntegral ci = Plain,
           ControlVariate cv = Without,
           Quadrature quad = ExpSinh)
        : GenericModelEngine<HestonModel,
              VanillaOption::arguments, VanillaOption::results>(hestonModel),
          precision_(precision),
          ci_(ci),
          cv_(cv),
          quad_(quad) {
            QL_REQUIRE(std::numeric_limits<T>::has_infinity,
                "value type does not support infinity");

            update();
        }

        void update() override {
            v0 = model_->v0();
            kappa = model_->kappa();
            theta = model_->theta();
            sigma = model_->sigma();
            rho = model_->rho();
        }

        std::complex<T> lnChF(std::complex<T> z, T t) const {
            const T sigma2 = sigma*sigma;

            const std::complex<T> g
                = kappa + T(rho*sigma)*std::complex<T>(z.imag(), -z.real());

            const std::complex<T> D = sqrt(
                g*g + (z*z + std::complex<T>(-z.imag(), z.real()))*sigma2);

            // reduce cancelation errors, see. L. Andersen and M. Lake
            std::complex<T> r(g-D);
            if (g.real()*D.real() + g.imag()*D.imag() > T(0)) {
                r = sigma2*z*std::complex<T>(-z.real(), -z.imag()-T(1))/(g+D);
            }

            std::complex<T> y;
            if (D.real() != T(0) || D.imag() != T(0)) {
                y = expm1(-D*t)/(T(2)*D);
            }
            else
                y = T(-0.5)*t;

            const std::complex<T> A
                = T(kappa*theta/sigma2)*(r*t - T(2)*log1p(-r*y));
            const std::complex<T> B
                = z*std::complex<T>(z.real(), z.imag() + T(1))*y/(T(1)-r*y);

            return (A + v0*B);
        }

        T operator()(
            const T& u, const T& freq,
            const T& tanPhi, const T& alpha, const T& t) const {

            const std::complex<T> hu(u, u*tanPhi - alpha);
            const std::complex<T> hPrime(hu.real(), hu.imag() - T(1));

            const std::complex<T> a(-u*tanPhi*freq, u*freq);

            std::complex<T> chfCv(0);
            if (cv_ == BlackScholes || cv_ == BlackScholes2)
                chfCv = Exp(
                    a - T(0.5)*vAvg*t*(hPrime*hPrime +
                            std::complex<T>(-hPrime.imag(), hPrime.real()))
                );

            return (
                std::complex<T>(T(1), tanPhi) *
                   (Exp(a + lnChF(hPrime, t)) - chfCv) / (hu * hPrime)
            ).real();
        }

        void calculate() const override {
            const ext::shared_ptr<HestonProcess>& process = model_->process();

            QL_REQUIRE(arguments_.exercise->type() == Exercise::European,
                       "not an European option");
            const Date maturityDate = arguments_.exercise->lastDate();
            const T t = process->time(maturityDate);

            ext::shared_ptr<PlainVanillaPayoff> payoff =
                ext::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
            QL_REQUIRE(payoff, "non plain vanilla payoff given");
            const Option::Type optionType = payoff->optionType();
            const T strike = payoff->strike();

            const DiscountFactor dr = process->riskFreeRate()->discount(maturityDate);
            const T fwd = process->s0()->value()
                 * process->dividendYield()->discount(maturityDate) / dr;

            const std::tuple<T, Size, Real> retVal =
                _calculate(optionType, fwd, strike, t, dr, v0, kappa, theta, sigma, rho);

            results_.value = std::get<2>(retVal);
            results_.additionalResults["value"] = std::get<0>(retVal);
            results_.additionalResults["function_calls"] = std::get<1>(retVal);
        }

        std::tuple<T, Size, Real> calculate(
            const Option::Type optionType,
            const std::string& spot, const std::string& strike,
            const std::string& t, const std::string& r, const std::string& q,
            const std::string& _v0, const std::string& _kappa,
            const std::string& _theta, const std::string& _sigma, const std::string& _rho) const {

            const T dr = exp(-parse_string(r)*parse_string(t));
            const T dq = exp(-parse_string(q)*parse_string(t));
            const T fwd = parse_string(spot) * dq / dr;

            return _calculate(
                optionType,
                fwd, parse_string(strike), parse_string(t), dr,
                parse_string(_v0), parse_string(_kappa), parse_string(_theta),
                parse_string(_sigma), parse_string(_rho)
            );
        }

        std::tuple<T, Size, Real> _calculate(
            const Option::Type optionType,
            const T& fwd, const T& strike, const T& t, const T& dr,
            const T& _v0, const T& _kappa, const T& _theta, const T& _sigma, const T& _rho) const {

            const T alpha(-0.5);
            const T freq = log(fwd/strike);

            v0 = _v0; kappa = _kappa; theta = _theta; sigma = _sigma; rho = _rho;
            if (cv_ == BlackScholes)
                vAvg = (T(1)-exp(-kappa*t))*(v0 - theta)/(kappa*t) + theta;
            else if (cv_ == BlackScholes2) {
                const T b = T(1) + alpha;
                vAvg = -T(2)*lnChF(std::complex<T>(0, -b), t).real()/(t*(b - b*b));
            }

            T tanPhi(0);
            if (ci_ == AngledContourShift) {
              const T r = rho - sigma*freq / (v0 + kappa*theta*t);
              tanPhi = tan((r*freq < 0.0)? m_pi/12*boost::math::sign(freq) : T(0));
            }

            Size fCalls(0);
            const auto integrant = [&](T u) -> T {
              ++fCalls;
              return this->operator()(u, freq, tanPhi, alpha, t);
            };

            T integral;
            switch (quad_) {
              case ExpSinh:
                integral = boost::math::quadrature::exp_sinh<T>(20)
                    .integrate(integrant, precision_);
                break;
              case TanhSinh:
                integral = boost::math::quadrature::tanh_sinh<T>(20)
                    .integrate(integrant, T(0), std::numeric_limits<T>::max(), precision_);
                break;
              case SinhSinh:
                integral = 0.5*boost::math::quadrature::sinh_sinh<T>(20)
                  .integrate([&](T u) -> T { return integrant(abs(u)); }, precision_);
                break;
              case GaussLaguerre:
                if (!laguerreQuadrature_)
                    laguerreQuadrature_ = ext::make_shared<GaussLaguerreQuadrature>(Size(precision_));

                integral = laguerreQuadrature_->integrate(integrant);
            }

            T cvValue;
            switch (cv_) {
              case Without:
                cvValue = ((alpha <=  T(0))
                          ? (optionType == Option::Call)? fwd : strike
                          : (optionType == Option::Call)? T(0) : strike - fwd)
                  -    ((alpha <= T(-1))? strike : T(0))
                  - T(0.5)*((alpha == T(0))? fwd : T(0))
                  + T(0.5)*((alpha == T(-1))? strike : T(0));
              break;
              case BlackScholes:
              case BlackScholes2:
              {
                const T v = sqrt(vAvg*t);
                const T d1 = (log(fwd/strike) + T(0.5)*vAvg*t)/v;
                const T d2 = d1 - v;
                const boost::math::normal_distribution<T> n(T(0), T(1));
                if (optionType == Option::Call)
                  cvValue = fwd*boost::math::cdf(n, d1) - strike*boost::math::cdf(n, d2);
                else
                  cvValue = strike*boost::math::cdf(n, -d2) - fwd*boost::math::cdf(n, -d1);
              }
              break;
            default:
              QL_FAIL("unknown control variate");
            }

            const T npv = dr
              * ( cvValue - fwd * integral * exp(alpha*freq) / boost::math::constants::pi<T>() );

            return std::make_tuple(npv, fCalls, convert_to_real(npv));
      }

      private:

        static Real convert_to_real(T x) {
            return x.template convert_to<Real>();
        }

        static T parse_string(const std::string& x) {
            return T(x);
        }

        const T m_pi = boost::math::constants::pi<T>();
        const T precision_;
        const ContourIntegral ci_;
        const ControlVariate cv_;
        const Quadrature quad_;
        mutable T v0, kappa, theta, sigma, rho;
        mutable T vAvg;

        static std::complex<T> Exp(const std::complex<T>& x) {
            if (x.real() != -std::numeric_limits<T>::infinity())
                return exp(x);
            else
                return std::complex<T>(0);
        };

        static std::complex<T> expm1(const std::complex<T>& z) {
            if (abs(z) < T(1)) {
                const T a = z.real();
                const T b = z.imag();
                const T exp_1 = boost::math::expm1(a);
                const T cos_1 = T(-2)*squared(sin(T(0.5)*b));

                return std::complex<T>(
                    exp_1*cos_1 + exp_1 + cos_1, sin(b)*exp(a)
                );
            }
            else {
                return exp(z)-T(1);
            }
        }

        static std::complex<T> log1p(const std::complex<T>& z) {
            const T a = z.real();
            const T b = z.imag();
            if (abs(a) < T(0.5) && abs(b) < T(0.5)) {
                return std::complex<T>(
                    T(0.5)*boost::math::log1p(a*a + T(2)*a + b*b),
                    std::arg(T(1) + z)
                );
            }
            else {
                return log(T(1)+z);
            }
        }

        class TqrEigenDecomposition {
          public:
            TqrEigenDecomposition(const std::vector<T>& diag,
                                  const std::vector<T>& sub)
            : iter_(0),
              d_(diag),
              ev_(d_.size(), 0.0) {
                Size n = diag.size();

                QL_REQUIRE(n == sub.size()+1, "Wrong dimensions");

                std::vector<T> e(n, 0.0);
                std::copy(sub.begin(),sub.end(),e.begin()+1);

                ev_[0] = 1.0;

                for (Size k=n-1; k >=1; --k) {
                    while (!offDiagIsZero(k, e)) {
                        Size l = k;
                        while (--l > 0 && !offDiagIsZero(l,e));
                        iter_++;

                        T q = d_[l];
                        // calculated eigenvalue of 2x2 sub matrix of
                        // [ d_[k-1] e_[k] ]
                        // [  e_[k]  d_[k] ]
                        // which is closer to d_[k+1].
                        const T t1 = sqrt(
                                              0.25*(d_[k]*d_[k] + d_[k-1]*d_[k-1])
                                              - 0.5*d_[k-1]*d_[k] + e[k]*e[k]);
                        const T t2 = 0.5*(d_[k]+d_[k-1]);

                        const T lambda =
                            (abs(t2+t1 - d_[k]) < abs(t2-t1 - d_[k]))?
                            T(t2+t1) : T(t2-t1);

                        q-=((k==n-1)? 1.25 : 1.0)*lambda;

                        // the QR transformation
                        T sine = 1.0;
                        T cosine = 1.0;
                        T u = 0.0;

                        bool recoverUnderflow = false;
                        for (Size i=l+1; i <= k && !recoverUnderflow; ++i) {
                            const T h = cosine*e[i];
                            const T p = sine*e[i];

                            e[i-1] = sqrt(p*p+q*q);
                            if (e[i-1] != T(0)) {
                                sine = p/e[i-1];
                                cosine = q/e[i-1];

                                const T g = d_[i-1]-u;
                                const T t = (d_[i]-g)*sine+2*cosine*h;

                                u = sine*t;
                                d_[i-1] = g + u;
                                q = cosine*t - h;

                                const T tmp = ev_[i-1];
                                ev_[i-1] = sine*ev_[i] + cosine*tmp;
                                ev_[i] = cosine*ev_[i] - sine*tmp;
                            } else {
                                // recover from underflow
                                d_[i-1] -= u;
                                e[l] = 0.0;
                                recoverUnderflow = true;
                            }
                        }

                        if (!recoverUnderflow) {
                            d_[k] -= u;
                            e[k] = q;
                            e[l] = 0.0;
                        }
                    }
                }

                // sort (eigenvalues, eigenvectors),
                // code taken from symmetricSchureDecomposition.cpp
                std::vector<std::pair<T, T> > temp(n);
                for (Size i=0; i<n; i++) {
                    temp[i] = std::make_pair(d_[i], ev_[i]);
                }
                std::sort(temp.begin(), temp.end(),
                          std::greater<std::pair<T, T> >());
                // first element is positive
                for (Size i=0; i<n; i++) {
                    d_[i] = temp[i].first;
                    if (temp[i].second < T(0))
                        ev_[i] = - temp[i].second;
                    else
                        ev_[i] = temp[i].second;
                }
            }

            const std::vector<T>& eigenvalues()  const {
                return d_;
            }
            const std::vector<T>& eigenvectors() const {
                return ev_;
            }

            Size iterations() const { return iter_; }

          private:

            // see NR for abort assumption as it is
            // not part of the original Wilkinson algorithm
            bool offDiagIsZero(Size k, const std::vector<T>& e) const {
                return abs(d_[k-1])+abs(d_[k])
                    == abs(d_[k-1])+abs(d_[k])+abs(e[k]);
            }

            Size iter_;
            std::vector<T> d_;
            std::vector<T> ev_;
        };

        class GaussLaguerreQuadrature {
          public:
            GaussLaguerreQuadrature(Size n)
            : x_(n), w_(n) {
                std::vector<T> e(n-1);

                for (Size i=1; i < n; ++i) {
                    e[i-1] = T(i);
                    x_[i] = 2*e[i-1]+1;
                }
                x_[0] = T(1);

                TqrEigenDecomposition tqr(x_, e);

                x_ = tqr.eigenvalues();
                const std::vector<T>& ev = tqr.eigenvectors();

                for (Size i=0; i<n; ++i)
                    w_[i] = squared(T(ev[i] * exp(0.5*x_[i])));
            }

            const std::vector<T>& weights() const { return w_; }
            const std::vector<T>& x()       const { return x_; }

            template<class F>
            T integrate(F f) const {
                T s(0);
                for (Integer i = x_.size()-1; i >= 0; --i)
                    s += w_[i] * f(x_[i]);

                return s;
            }

          protected:
            std::vector<T> x_, w_;
        };

        mutable ext::shared_ptr<GaussLaguerreQuadrature> laguerreQuadrature_;
    };

    template <> inline
    Real MultiPrecisionHestonEngine<Real>::convert_to_real(Real x) {
        return x;
    }

    template <> inline
    Real MultiPrecisionHestonEngine<Real>::parse_string(const std::string& x) {
        return std::stod(x);
    }

}
#endif
