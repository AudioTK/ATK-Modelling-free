/**
 * \file StaticDiode.h
 */

#ifndef ATK_MODELLING_STATICDIODE_H
#define ATK_MODELLING_STATICDIODE_H

#include <ATK/Modelling/SPICE/Utilities.h>
#include <ATK/Utility/fmath.h>

namespace ATK
{
  /// Diode component
  template<typename DataType_, unsigned int direct = 1, unsigned int indirect = 0>
  class StaticDiode
  {
  public:
    using DataType = DataType_;

    StaticDiode(DataType Is = 1e-14, DataType N = 1.24, DataType Vt = 26e-3): Is(Is), Vt(N * Vt)
    {
    }
    
    /**
     * Get current
     */
    DataType get_current() const
    {
      return Is * ((precomp - 1) - (invprecomp - 1));
    }
    
    /**
     * Get current gradient
     */
    DataType get_gradient() const
    {
      if constexpr(indirect != 0)
      {
        return Is / Vt * (precomp / direct + invprecomp / indirect);
      }
      else
      {
        return Is / Vt * precomp / direct;
      }
    }
    
    /**
     * Precompute internal value before asking current and gradients
     */
    void precompute(DataType V0, DataType V1) const
    {
      precomp = fmath::exp((V1 - V0) / (Vt * direct));
      if constexpr(direct == indirect)
      {
        invprecomp = 1 / precomp;
      }
      else if constexpr(indirect != 0)
      {
        invprecomp = fmath::exp((V0 - V1) / (Vt * indirect));
      }
    }

  private:
    DataType Is;
    DataType Vt;
    mutable DataType precomp{0};
    mutable DataType invprecomp{1};
  };

#define DIODE_SEQ ((vt, 26e-3))((is, 1e-14))((n, 1.24))
  HELPER(DiodeHelper, DIODE_SEQ, Diode)
}

#endif
