#pragma once

#include <aspect/simulator_access.h>
#include <aspect/heating_model/interface.h>

namespace aspect {
  namespace HeatingModel {
    template <int dim>
    class VeMel25 : public Interface<dim>, public ::aspect::SimulatorAccess<dim> {
      public:
        void evaluate ( const MaterialModel::MaterialModelInputs<dim> & material_model_inputs,
                        const MaterialModel::MaterialModelOutputs<dim> & material_model_outputs,
                        HeatingModel::HeatingModelOutputs & heating_model_outputs ) const override;
        static void declare_parameters ( ParameterHandler & prm );
        void parse_parameters ( ParameterHandler & prm ) override;
        
        double latent_heat;
        double required_melting_precision;
        int maximum_melting_nonlinear_iterations;
      private:
        mutable std::shared_ptr<const InitialComposition::Manager<dim>> initcomp;
    }
  }
}
