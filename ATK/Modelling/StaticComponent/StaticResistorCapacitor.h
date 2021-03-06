/**
 * \file StaticStaticResistorCapacitor.h
 */

#ifndef ATK_MODELLING_STATICRESISTORCAPACITOR_H
#define ATK_MODELLING_STATICRESISTORCAPACITOR_H

namespace ATK
{
/// Capacitor component
template <typename DataType_>
class StaticResistorCapacitor
{
public:
  using DataType = DataType_;

  StaticResistorCapacitor(DataType R, DataType C): R(R), C(C)
  {
  }

  /**
   * Update the component for its steady state condition
   * @param dt is the delat that will be used in following updates
   */
  void update_steady_state(DataType dt, DataType V0, DataType V1)
  {
    c2t = (2 * C) / dt;
    iceq = c2t * (V1 - V0);
  }

  /**
   * Update the component for its current state condition
   */
  void update_state(DataType V0, DataType V1) const
  {
    iceq = 2 * (iceq + ((V1 - V0) * c2t - iceq) / (1 + c2t * R)) - iceq;
  }

  /**
   * Get current for the given voltages
   */
  DataType_ get_current(DataType V0, DataType V1) const
  {
    return ((V1 - V0) * c2t - iceq) / (1 + c2t * R);
  }

  /**
   * Get current gradient for the given voltages
   */
  DataType_ get_gradient() const
  {
    return c2t / (1 + c2t * R);
  }

private:
  DataType R;
  DataType C;
  DataType c2t{0};
  mutable DataType iceq{0};
};
} // namespace ATK

#endif
