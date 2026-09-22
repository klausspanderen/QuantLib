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

#include "toplevelfixture.hpp"
#include "utilities.hpp"

#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/daycounters/yearfractiontodate.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/pricingengines/vanilla/exponentialfittinghestonengine.hpp>
#include <ql/pricingengines/vanilla/multiprecisionhestonengine.hpp>

#include <ql/math/integrals/gaussianquadratures.hpp>
#include <ql/math/integrals/gausslaguerrecosinepolynomial.hpp>


#include <sstream>

using namespace QuantLib;
using namespace boost::unit_test_framework;


BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(MultiPrecisionHestonEngineTests)


BOOST_AUTO_TEST_CASE(testMultiPrecisionHestonPricing) {
    BOOST_TEST_MESSAGE("Testing mutli-precision Heston engine...");

    const auto dbl2str = [](Real x) -> std::string {
        std::stringstream ss;
        ss << std::setprecision(8) << x;
        return ss.str();
    };

    DayCounter dayCounter = Actual365Fixed();
    Date settlementDate = Settings::instance().evaluationDate();

    const Rate r = 0.01, q = 0.02;
    Handle<YieldTermStructure> riskFreeTS(flatRate(r, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(q, dayCounter));

    Handle<Quote> s0(ext::make_shared<SimpleQuote>(100.0));

    const Real kappa = 4.0;
    const Real theta = 0.25;
    const Real sigma = 1.0;
    const Real rho = -0.5;

    // https://financepress.com/2019/02/15/heston-model-reference-prices/
    const std::vector<std::tuple<std::string, Real, Real, Real, Option::Type> > testCases = {
        { "7.958878113256768285213263077598987193482161301733", 80, 0.04, 1.0, Option::Put},
        { "26.774758743998854221382195325726949201687074848341", 80, 0.04, 1.0, Option::Call},
        { "12.017966707346304987709573290236471654992071308187", 90, 0.04, 1.0, Option::Put},
        { "20.933349000596710388139445766564068085476194042256", 90, 0.04, 1.0, Option::Call},
        { "17.055270961270109413522653999411000974895436309183", 100, 0.04, 1.0, Option::Put},
        { "16.070154917028834278213466703938231827658768230714", 100, 0.04, 1.0, Option::Call},
        { "23.017825898442800538908781834822560777763225722188", 110, 0.04, 1.0, Option::Put},
        { "12.132211516709844867860534767549426052805766831181", 110, 0.04, 1.0, Option::Call},
        { "29.811026202682471843340682293165857439167301370697", 120, 0.04, 1.0, Option::Put},
        { "9.024913483457835636553375454092357136489051667150", 120, 0.04, 1.0, Option::Call},
        { "4.5183603586861772614990106188215872180542e-8", 90, 0.01, 0.01, Option::Put},
        { "9.989001595065276544935948045293485530832966049263", 90, 0.01, 0.01, Option::Call},
        { "0.000461954855653851579672612557018857858641926937", 95, 0.01, 0.01, Option::Put},
        { "4.989963479738160122154264702582719627807098780529", 95, 0.01, 0.01, Option::Call},
        { "0.477781171629504680023239655436072890669645669297", 100, 0.01, 0.01, Option::Put},
        { "0.467782671512844263098248405184095087949465507760", 100, 0.01, 0.01, Option::Call},
        { "5.009501052563650299130635110520904481889436667608", 105, 0.01, 0.01, Option::Put},
        { "2.527447823194706060519991248106500619490942e-6", 105, 0.01, 0.01, Option::Call},
        { "10.008998550115123724684210555728039829315964456261", 110, 0.01, 0.01, Option::Put},
        { "1.29932760052624920704881258510264466e-13", 110, 0.01, 0.01, Option::Call}
    };

    using MP_Real = typename boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100> >;
    typedef MultiPrecisionHestonEngine<MP_Real> MultiPrecisonEngineType;

    const ext::shared_ptr<HestonModel> hestonModel =
        ext::make_shared<HestonModel>(
            ext::make_shared<HestonProcess>(
                riskFreeTS, dividendTS, s0, theta, kappa, theta, sigma, rho
            )
        );

    const ext::shared_ptr<MultiPrecisonEngineType> mp_engine
        = ext::make_shared<MultiPrecisonEngineType>(
            hestonModel,
            1000,
            MultiPrecisonEngineType::Plain,
            MultiPrecisonEngineType::Without,
            MultiPrecisonEngineType::GaussLaguerre);

    const ext::shared_ptr<AnalyticHestonEngine> real_engine =
        ext::make_shared<AnalyticHestonEngine>(
            hestonModel,
            AnalyticHestonEngine::AngledContour,
            //AnalyticHestonEngine::Integration::gaussLobatto(1e-8, 1e-8)
            AnalyticHestonEngine::Integration::gaussLaguerre(32),
            1e-25, -0.5
        );

    for (const auto& p: testCases) {
        const Real strike = std::get<1>(p);
        const Real v0 = std::get<2>(p);
        const Real t = std::get<3>(p);
        const Option::Type optionType = std::get<4>(p);

        const auto retVal = mp_engine->calculate(
            optionType,
            dbl2str(s0->value()), dbl2str(strike),
            dbl2str(t), dbl2str(r), dbl2str(q),
            dbl2str(v0), dbl2str(kappa), dbl2str(theta), dbl2str(sigma), dbl2str(rho)
        );

        const MP_Real result = std::get<0>(retVal);
        const MP_Real mpExpected = MP_Real(std::get<0>(p));
        const Size fctCalls = std::get<1>(retVal);

        Date exerciseDate = yearFractionToDate(dayCounter, settlementDate, t);

        const ext::shared_ptr<Exercise> exercise(
            ext::make_shared<EuropeanExercise>(exerciseDate));

        const ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::make_shared<PlainVanillaPayoff>(optionType, strike);

        VanillaOption option(payoff, exercise);

        Array params = hestonModel->params();
        params[4] = v0;
        hestonModel->setParams(params);

        option.setPricingEngine(real_engine);
        const Real calculated = option.NPV();
        std::cout << std::setprecision(8)
                << (mpExpected - result) << " "
                << fctCalls << " "
                << mpExpected - calculated << " "
                << real_engine->numberOfEvaluations() << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(testGaussLaguerreCosineIntegration) {
    BOOST_TEST_MESSAGE("Testing Gauss-Laguerre-Cosine integration...");

    const DayCounter dc = Actual365Fixed();
    const Date today = Date(17, September, 2025);
    const Date maturity = today + Period(2, Years);

    const Time tau = dc.yearFraction(today, maturity);

    const Real K = 170;
    const Real fwd = 100;
    const Real vol = 0.01;
    const Real kappa = 0.001;
    const Real theta = 0.04;
    const Real sigma = 0.33147;
    const Real rho = 0.5258519;
    const Rate r = 0;

    const auto ts = Handle<YieldTermStructure>(flatRate(r, dc));
    const auto model = ext::make_shared<HestonModel>(
        ext::make_shared<HestonProcess>(
            ts, ts,
            Handle<Quote>(ext::make_shared<SimpleQuote>(fwd)),
            vol, kappa, theta, sigma, rho
        )
    );

    const auto engine = ext::make_shared<AnalyticHestonEngine>(
        model,
        AnalyticHestonEngine::AndersenPiterbarg,
        AnalyticHestonEngine::Integration::gaussLaguerre(192)
    );

    const auto opt_alpha = ext::make_shared<AnalyticHestonEngine::OptimalAlpha>(tau, engine.get());


}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
