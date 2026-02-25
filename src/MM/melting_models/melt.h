#pragma once

#include <cmath>
#include <utility>

namespace VeMel25_melting_models {
  const double cSol [3] = { 1085.7, 132.9e-9, -5.1e-18 };
  const double cLhe [3] = { 1475.0, 80.0e-9,  -3.2e-18 };
  const double cLiq [3] = { 1780.0, 45.0e-9,  -2.0e-18 };
  const double CtoK = 273.14;
  namespace McKenzie1988nmspc {
    const double a = 0.5;
    const double b = 0.25;
    const double c = 0.4259;
    const double d = 2.988;
  }
  namespace Katz2003nmspc {
    const double r1 = 0.5;
    const double r2 = 0.08e-9;
    const double b1 = 1.5;
    const double b2 = 1.5;
  }

  inline double T_Solidus ( const double p ) {
    return cSol[0] + cSol[1]*p + cSol[2]*p*p + CtoK;
  }
  inline double T_Lherzolite ( const double p ) {
    return cLhe[0] + cLhe[1]*p + cLhe[2]*p*p + CtoK;
  }
  inline double T_Liquidus ( const double p ) {
    return cLiq[0] + cLiq[1]*p + cLiq[2]*p*p + CtoK;
  }

  inline double McKenzie1988 ( const double T, const double T_sol, const double T_liq ) {
    using namespace VeMel25_melting_models::McKenzie1988nmspc;
    const double Tss = std::max ( -0.5, std::min ( 0.5, ( T - (T_sol + T_liq)/2.0 ) / ( T_liq - T_sol ) ) );
    return a + Tss + ( Tss*Tss - b ) * ( c + d*Tss );
  }

  inline double F_cpxout ( const double p, const double Cpx_frac ) {
    using namespace VeMel25_melting_models::Katz2003nmspc;
    return Cpx_frac/(r1 + r2*p);
  }
  inline double T_cpxout ( const double F_cpxout, const double T_solidus, const double T_lherzolite ) {
    using namespace VeMel25_melting_models::Katz2003nmspc;
    return pow ( F_cpxout, 1.0/b1 ) * ( T_lherzolite - T_solidus ) + T_solidus;
  }

  inline double Katz2003 ( const double T, const double T_sol, const double T_lhe, const double T_liq, const double F_cpxout, const double T_cpxout ) {
    using namespace VeMel25_melting_models::Katz2003nmspc;
    if ( T <= T_sol )
      return 0.0;
    else if ( T_liq <= T )
      return 1.0;

    const double F_cpx = pow ( ( T - T_sol ) / ( T_lhe - T_sol ), b1 );
    if ( F_cpx < F_cpxout )
      return F_cpx;
    else {
      const double F_opx = F_cpxout + ( 1.0 - F_cpxout ) * pow (  (T - T_cpxout) / (T_liq - T_cpxout), b2 );
      return F_opx;
    }
  }

  inline std::pair<double,double> Katz2003der ( const double T, const double T_sol, const double T_lhe, const double T_liq, const double F_cpxout, const double T_cpxout ) {
    using namespace VeMel25_melting_models::Katz2003nmspc;
    if ( T <= T_sol )
      return std::make_pair(0.0, 0.0);
    else if ( T_liq <= T )
      return std::make_pair(1.0, 0.0);

    const double F_cpx = pow ( ( T - T_sol ) / ( T_lhe - T_sol ), b1 );
    if ( F_cpx < F_cpxout ) {
      const double der_F_cpx = 1.5*b1 * pow ( ( T - T_sol ) / ( T_lhe - T_sol ), b1-1.0 ) / ( T_lhe - T_sol );
      return std::make_pair ( F_cpx, der_F_cpx );
    }
    else {
      const double F_opx = F_cpxout + ( 1.0 - F_cpxout ) * pow (  (T - T_cpxout) / (T_liq - T_cpxout), b2 );
      const double der_F_opx = 1.5*b2 * ( 1.0 - F_cpxout ) * pow (  (T - T_cpxout) / (T_liq - T_cpxout), b2-1.0 ) / (T_liq - T_cpxout);
      return std::make_pair ( F_opx, der_F_opx );
    }
  }
}
