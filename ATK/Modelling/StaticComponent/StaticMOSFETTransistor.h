/**
 * \file StaticMOSFETTransistor.h
 */

#ifndef ATK_MODELLING_STATICMOSFETTRANSISTOR_H
#define ATK_MODELLING_STATICMOSFETTRANSISTOR_H

#include <ATK/Modelling/SPICE/Utilities.h>
#include <ATK/Utility/fmath.h>

namespace ATK
{
/// Transistor NPN component
template <typename DataType_>
class StaticNMOS
{
public:
  using DataType = DataType_;

  StaticNMOS(DataType Kp = 8e-6, DataType Vto = 2, DataType lambda = 0, DataType W = 250e-6, DataType L = 10e-6)
    : Kp(Kp * W / L), Vto(Vto), lambda(lambda)
  {
  }

private:
  const DataType Kp;
  const DataType Vto;
  const DataType lambda;

public:
  DataType is(DataType V0, DataType V1, DataType V2) const
  {
    DataType Vgs = (V0 - V2);
    DataType Vds = (V1 - V2);
    if(Vgs < Vto)
    {
      return 0;
    }
    else
    {
      if(Vds < Vgs - Vto)
      {
        return Kp * (1 + lambda * Vds) * ((Vgs - Vto) * Vds - Vds * Vds / 2);
      }
      else
      {
        return Kp * (1 + lambda * Vds) / 2 * (Vgs - Vto) * (Vgs - Vto);
      }
    }
  }

  DataType is_Vgs(DataType V0, DataType V1, DataType V2) const
  {
    DataType Vgs = (V0 - V2);
    DataType Vds = (V1 - V2);
    if(Vgs < Vto)
    {
      return 0;
    }
    else
    {
      if(Vds < Vgs - Vto)
      {
        return Kp * (1 + lambda * Vds) * Vds;
      }
      else
      {
        return Kp * (1 + lambda * Vds) * (Vgs - Vto);
      }
    }
  }
  DataType is_Vds(DataType V0, DataType V1, DataType V2) const
  {
    DataType Vgs = (V0 - V2);
    DataType Vds = (V1 - V2);
    if(Vgs < Vto)
    {
      return 0;
    }
    else
    {
      if(Vds < Vgs - Vto)
      {
        return Kp * (lambda * ((Vgs - Vto) * Vds - Vds * Vds / 2) + (1 + lambda * Vds) * ((Vgs - Vto) - Vds));
      }
      else
      {
        return Kp * lambda / 2 * (Vgs - Vto) * (Vgs - Vto);
      }
    }
  }
};

/// Transistor PNP component
template <typename DataType_>
class StaticPMOS
{
public:
  using DataType = DataType_;

  StaticPMOS(DataType Kp = 8e-6, DataType Vto = -2, DataType lambda = 0, DataType W = 250e-6, DataType L = 10e-6)
    : Kp(Kp * W / L), Vto(Vto), lambda(lambda)
  {
  }

private:
  const DataType Kp;
  const DataType Vto;
  const DataType lambda;

public:
  DataType is(DataType V0, DataType V1, DataType V2) const
  {
    DataType Vsg = -(V0 - V2);
    DataType Vsd = -(V1 - V2);
    if(Vsg < Vto)
    {
      return 0;
    }
    else
    {
      if(Vsd < Vsg + Vto)
      {
        return -Kp * (1 + lambda * Vsd) * ((Vsg + Vto) * Vsd - Vsd * Vsd / 2);
      }
      else
      {
        return -Kp * (1 + lambda * Vsd) / 2 * (Vsg + Vto) * (Vsg + Vto);
      }
    }
  }
  DataType is_Vgs(DataType V0, DataType V1, DataType V2) const
  {
    DataType Vsg = -(V0 - V2);
    DataType Vsd = -(V1 - V2);
    if(Vsg < Vto)
    {
      return 0;
    }
    else
    {
      if(Vsd < Vsg + Vto)
      {
        return Kp * (1 + lambda * Vsd) * Vsd;
      }
      else
      {
        return Kp * (1 + lambda * Vsd) * (Vsg + Vto);
      }
    }
  }
  DataType is_Vds(DataType V0, DataType V1, DataType V2) const
  {
    DataType Vsg = -(V0 - V2);
    DataType Vsd = -(V1 - V2);
    if(Vsg < Vto)
    {
      return 0;
    }
    else
    {
      if(Vsd < Vsg + Vto)
      {
        return Kp * (lambda * ((Vsg + Vto) * Vsd - Vsd * Vsd / 2) + (1 + lambda * Vsd) * ((Vsg + Vto) - Vsd));
      }
      else
      {
        return Kp * lambda / 2 * (Vsg + Vto) * (Vsg + Vto);
      }
    }
  }
};

#define NMOSFET_SEQ ((Kp, 8e-6))((Vto, 2))((lambda, 0))((W, 250e-6))((L, 10e-6))
#define PMOSFET_SEQ ((Kp, 8e-6))((Vto, -2))((lambda, 0))((W, 250e-6))((L, 10e-6))
HELPER(NMOSFETHelper, NMOSFET_SEQ, Static)
HELPER(PMOSFETHelper, PMOSFET_SEQ, Static)
} // namespace ATK

#endif
