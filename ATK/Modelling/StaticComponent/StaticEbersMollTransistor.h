/**
 * \file StaticEbersMollTransistor.h
 */

#ifndef ATK_MODELLING_STATICEBERSMOLLTRANSISTOR_H
#define ATK_MODELLING_STATICEBERSMOLLTRANSISTOR_H

#include <ATK/Modelling/SPICE/Utilities.h>
#include <ATK/Utility/fmath.h>

namespace ATK
{
/// Transistor NPN component
template <typename DataType_>
class StaticEBNPN
{
public:
  using DataType = DataType_;

  StaticEBNPN(DataType Is = 1e-12, DataType Vt = 26e-3, DataType Ne = 1, DataType Br = 1, DataType Bf = 100)
    : Is(Is), Vt(Vt * Ne), Br(Br), Bf(Bf)
  {
  }

  /**
   * Precompute internal value before asking current and gradients
   */
  void precompute(DataType V0, DataType V1, DataType V2) const
  {
    expVbe = fmath::exp((V0 - V2) / Vt);
    expVbc = fmath::exp((V0 - V1) / Vt);
  }

private:
  const DataType Is;
  const DataType Vt;
  const DataType Br;
  const DataType Bf;
  mutable DataType expVbe;
  mutable DataType expVbc;

public:
  DataType ib() const
  {
    return Is * ((expVbe - 1) / Bf + (expVbc - 1) / Br);
  }
  DataType ic() const
  {
    return Is * ((expVbe - expVbc) - (expVbc - 1) / Br);
  }
  DataType ib_Vbc() const
  {
    return Is * expVbc / Vt / Br;
  }
  DataType ib_Vbe() const
  {
    return Is * expVbe / Vt / Bf;
  }
  DataType ic_Vbc() const
  {
    return Is * (-expVbc - expVbc / Br) / Vt;
  }
  DataType ic_Vbe() const
  {
    return Is * expVbe / Vt;
  }
};

/// Transistor PNP component
template <typename DataType_>
class StaticEBPNP
{
public:
  using DataType = DataType_;

  StaticEBPNP(DataType Is = 1e-12, DataType Vt = 26e-3, DataType Ne = 1, DataType Br = 1, DataType Bf = 100)
    : Is(Is), Vt(Vt * Ne), Br(Br), Bf(Bf), expVbe(0), expVbc(0)
  {
  }

  /**
   * Precompute internal value before asking current and gradients
   */
  void precompute(DataType V0, DataType V1, DataType V2) const
  {
    expVbe = fmath::exp(-(V0 - V2) / Vt);
    expVbc = fmath::exp(-(V0 - V1) / Vt);
  }

private:
  const DataType Is;
  const DataType Vt;
  const DataType Br;
  const DataType Bf;
  mutable DataType expVbe{0};
  mutable DataType expVbc{0};

public:
  DataType ib() const
  {
    return -Is * ((expVbe - 1) / Bf + (expVbc - 1) / Br);
  }
  DataType ic() const
  {
    return -Is * ((expVbe - expVbc) - (expVbc - 1) / Br);
  }
  DataType ib_Vbc() const
  {
    return Is * expVbc / Vt / Br;
  }
  DataType ib_Vbe() const
  {
    return Is * expVbe / Vt / Bf;
  }
  DataType ic_Vbc() const
  {
    return Is * (-expVbc - expVbc / Br) / Vt;
  }
  DataType ic_Vbe() const
  {
    return Is * expVbe / Vt;
  }
};

#define EBERMOLLSBIPOLAR_SEQ ((is, 1e-12))((vt, 26e-3))((ne, 1))((br, 1))((bf, 100))
HELPER(EbersMollBipolarHelper, EBERMOLLSBIPOLAR_SEQ, StaticEB)
} // namespace ATK

#endif
