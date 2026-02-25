#include "HM.h"
#include "../MM/MM.h"

namespace aspect {
  namespace HeatingModel {
    template <int dim>
    void VeMel25<dim>::evaluate ( const MaterialModel::MaterialModelInputs<dim> & min,
                                  const MaterialModel::MaterialModelOutputs<dim> & mout,
                                  HeatingModel::HeatingModelOutputs & hout ) const {
      AssertThrow ( this->introspection().compositional_name_exists("mantle_depletion") == true,
        ExcMessage ( "mantle_depletion composition field not found." ) );
      for ( unsigned int q = 0; q < hout.heating_source_terms.size(); ++q ) {
        const unsigned int mantle_depletion_indx = this->introspection().compositional_index_for_name ( "mantle_depletion" );
        const double dF = mout.reaction_terms[q][mantle_depletion_indx];
        if ( this->get_timestep_number() > 1 ) {
          hout.heating_source_terms[q] - out.densities[q]*latent_heat*dF/this->get_timestep();
        }
        else if ( this->get_timestep_number() == 1 ) {
          const double F0 = initcomp.get()->initial_composition ( in.position[q], mantle_depletion_indx );
          hout.heating_source_terms[q] = -out.densities[q]*latent_heat*(in.composition[q][mantle_depletion_indx]-F0 + dF)/this->get_timestep();
        }
        else {
          hout.heating_source_terms[q] = 0.0;
          if ( initcomp.get() == nullptr )
            initcomp = this->get_initial_composition_manager_pointer();
        }
        hout.lhs_latent_heat_terms[q] = 0.0;
      }
    }

    template <int dim>
    void VeMel25<dim>::declare_parameters ( ParameterHandler & prm ) {
      prm.enter_subsection ( "Heating model" );
        prm.enter_subsection ( "VeMel25" );
          prm.declare_entry ( "Latent heat of melting", "325e3", Patterns::Double(0.0),
                              "The value of latent heat of melting L = c_p \\cdot \\delta T, units: J/kg" );
          prm.declare_entry ( "Required melting precision", "1e-6", Patterns::Double(0.),
                              "Required melting precision when computing melting fraction with the latent heat effects in nonlinear iterations. No units." );
          prm.declare_entry ( "Maximum melting nonlinear iterations", "10", Patterns::Integer(0),
                              "Maximum number of nonlinear iterations when computing melt fraction with latent heat effects. No units." );
        prm.leave_subsection ();
      prm.leave_subsection ();
    }

    template <int dim>
    void VeMel25::parse_parameters ( ParameterHandler & prm ) {
      prm.enter_subsection ( "Heating model" );
        prm.enter_subsection ( "VeMel25" );
          latent_heat                             = prm.get_double ( "Latent heat of melting" );
          required_melting_precision              = prm.get_double ( "Required melting precision" );
          maximum_melting_nonlinear_iterations    = prm.get_double ( "Maximum melting nonlinear iterations" );
        prm.leave_subsection ();
      prm.leave_subsection ();
    }
  }
}

namespace aspect {
  namespace HeatingModel {
    ASPECT_REGISTER_HEATING_MODEL ( VeMel25,
                                    "VeMel25",
                                    "Heating model by A. Dizov for the VeMel25 material model" )
  }
}
