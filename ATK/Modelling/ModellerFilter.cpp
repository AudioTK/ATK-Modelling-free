/**
 * \file ModellerFilter.cpp
 */

#include "ModellerFilter.h"

#include <ATK/Modelling/Component.h>

#include <ATK/Core/Utilities.h>

#if ENABLE_LOG
#define BOOST_LOG_DYN_LINK
#include <boost/log/trivial.hpp>
#endif

namespace ATK
{
template <typename DataType_>
ModellerFilter<DataType_>::ModellerFilter(gsl::index nb_dynamic_pins, gsl::index nb_input_pins)
  : TypedBaseFilter<DataType_>(nb_input_pins, nb_dynamic_pins)
{
}

template <typename DataType_>
gsl::index ModellerFilter<DataType_>::find_dynamic_pin(const std::string& name)
{
  for(gsl::index i = 0; i < get_nb_dynamic_pins(); ++i)
  {
    if(get_dynamic_pin_name(i) == name)
    {
      return i;
    }
  }
  throw ATK::RuntimeError("Unknown dynamic pin " + name);
}

template <typename DataType_>
gsl::index ModellerFilter<DataType_>::find_input_pin(const std::string& name)
{
  for(gsl::index i = 0; i < get_nb_input_pins(); ++i)
  {
    if(get_input_pin_name(i) == name)
    {
      return i;
    }
  }
  throw ATK::RuntimeError("Unknown input pin " + name);
}

template <typename DataType_>
gsl::index ModellerFilter<DataType_>::find_static_pin(const std::string& name)
{
  for(gsl::index i = 0; i < get_nb_static_pins(); ++i)
  {
    if(get_static_pin_name(i) == name)
    {
      return i;
    }
  }
  throw ATK::RuntimeError("Unknown static pin " + name);
}

template class ModellerFilter<double>;
}
