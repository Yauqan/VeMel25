#pragma once

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class VeMel25 : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        enum MeltingModel {
          none,
          McKenzie1988,
          Katz2003
        };

        bool is_compressible () const override;
        void evaluate ( const MaterialModel::MaterialModelInputs<dim> & in, MaterialModel::MaterialModelOutputs<dim> & out ) const override;
        static void declare_parameters ( ParameterHandler & prm );
        void parse_parameters ( ParameterHandler & prm ) override;

        double depletion_increment ( const double T, const double p, const double F_0, const double cp, const double lat_heat, const int maximum_melting_nonlinear_iterations = 1000, const double required_melting_precision = 1e-6 ) const;

      private:
        double volatile_partition;
        double viscosity_minimum;
        double viscosity_maximum;
        double thermal_conductivity;
        MeltingModel melting_model;
        double Cpx;
        double maxPressure = 10e9; // Pa s
        Rheology::DiffusionCreep<dim> diffusion_creep;
        EquationOfState::MulticomponentIncompressible<dim> equation_of_state;
    };

  }
}
