
#include "MM.h"
#include "../HM/HM.h"
#include "melting_models/melt.h"

#include <deal.II/numerics/fe_field_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim> double VeMel25<dim>::depletion_increment ( const double T, const double p, const double F_0, const double cp, const double lat_heat, const int maximum_melting_nonlinear_iterations, const double required_melting_precision ) const {
      using namespace VeMel25_melting_models;
      if ( p > maxPressure )
        return 0.0;
      
      if ( melting_model == MeltingModel::none )
        return 0.0;
      
      const double Tsol = T_Solidus(p);
      if ( T < Tsol )
        return 0.0;
      const double Tliq = T_Liquidus(p);

      if ( melting_model == MeltingModel::Katz2003 ) {
        const double Tlhe = T_Lherzolite(p);
        const double dTf = lat_heat / cp;
        const double Fcpxout = F_cpxout ( p, Cpx );
        const double Tcpxout = T_cpxout ( Fcpxout, Tsol, Tlhe );
        double dF = 0.0;
        for ( int k = 0; k <= maximum_melting_nonlinear_iterations; ++k ) {
          const std::pair<double, double> melting = Katz2003der ( T - dTf*dF, Tsol, Tlhe, Tliq, Fcpxout, Tcpxout );
          if ( abs ( melting.first - F_0 - dF ) < required_melting_precision )
            break;
          dF = dF - std::max(-0.1,std::min(0.1, ( melting.first - F_0 - dF ) / ( -dTf * melting.second - 1.0 )));
        }
        return std::max (dF, 0.0);
      }
      else if ( melting_model == MeltingModel::McKenzie1988 ) {
        AssertThrow ( lat_heat == 0.0, ExcMessage ( "McKenzie 1988 with latent heat was not implemented." ) );
        return std::max( 0.0, VeMel25_melting_models::McKenzie1988 ( T, Tsol, Tliq ) - F_0 );
      }
    }

    template <int dim> void VeMel25<dim>:: evaluate ( const MaterialModel::MaterialModelInputs<dim> & in, MaterialModel::MaterialModelOutputs<dim> & out ) const
    {
      AssertThrow ( this->introspection().compositional_name_exists("mantle_depletion") == true,
        ExcMessage ( "mantle_depletion composition field not found." ) );
      AssertThrow ( this->get_parameters().use_operator_splitting == false, 
        ExcMessage ( "operator splitting not supported." ) );

      // prepare outputs
      const unsigned int mantle_depletion_indx = this->introspection().compositional_index_for_name ( "mantle_depletion" );
      EquationOfStateOutputs<dim> eos_outputs ( this->introspection().n_chemical_composition_fields()+1 );
      std::vector<double> phase_function_values ( this->introspection().n_chemical_composition_fields()+1, 1.0 );
      std::vector<unsigned int> n_phase_transitions_per_composition ( this->introspection().n_chemical_composition_fields()+1, 0 );

      for ( unsigned int q = 0; q < in.n_evaluation_points(); ++q ) {
        equation_of_state.evaluate(in, q, eos_outputs);
        const double p = in.pressure[q];
        const double T = in.temperature[q];
        const std::vector<double> volume_fractions_old = MaterialUtilities::compute_only_composition_fractions ( in.composition[q], this->introspection().chemical_composition_field_indices() );
        const double spec_heat_old = MaterialUtilities::average_value ( volume_fractions_old, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic );
        
        const aspect::HeatingModel::VeMel25<dim> * HM = nullptr;
        for ( const std::unique_ptr< HeatingModel::Interface< dim > > & hm: this->get_heating_model_manager().get_active_heating_models() )
          if ( ( HM = dynamic_cast<const aspect::HeatingModel::VeMel25<dim>*>( hm.get() ) ) != nullptr )
            break;
        const double lat_heat = ( HM != nullptr ) ? ( HM->latent_heat ) : 0.0;
        const double mmni = ( HM != nullptr ) ? ( HM->maximum_melting_nonlinear_iterations ) : 1;
        const double rmp = ( HM != nullptr ) ? ( HM->required_melting_precision ) : 1e-6;

        const double Fold = in.composition[q][mantle_depletion_indx];
        const double dF = depletion_increment ( T, p, Fold, spec_heat_old, lat_heat, mmni, rmp );
        const double Fnew = Fold + dF;

        // Compute depletion source
        out.reaction_terms[q][mantle_depletion_indx] = dF;
        for ( unsigned int j = 0; j < this->n_compositional_fields()+1; ++j )
          if ( j != mantle_depletion_indx )
            out.reaction_terms[q][j] = 0.0;

        // Output material properties
        const std::vector<double> volume_fractions = MaterialUtilities::compute_only_composition_fractions ( in.composition[q], this->introspection().chemical_composition_field_indices() );
        
        std::vector<double> viscosities ( volume_fractions.size(), 0.0 );
        for ( unsigned int j = 0; j < volume_fractions.size(); ++j )
          viscosities[j] = diffusion_creep.compute_viscosity ( p, T, j, phase_function_values, n_phase_transitions_per_composition );
        out.viscosities[q]  = MaterialUtilities::average_value ( volume_fractions, viscosities, MaterialUtilities::harmonic );
        out.viscosities[q] *= std::pow ( 1.0 - Fnew, ( volatile_partition - 1.0 ) / volatile_partition );
        out.viscosities[q]  = std::min ( viscosity_maximum, std::max ( viscosity_minimum, out.viscosities[q] ) );
        out.densities[q] = MaterialUtilities::average_value ( volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic );
        out.thermal_expansion_coefficients[q] = MaterialUtilities::average_value ( volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic );
        out.specific_heat[q] = MaterialUtilities::average_value ( volume_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic );
        out.thermal_conductivities[q] = thermal_conductivity;
        out.entropy_derivative_pressure[q] = MaterialUtilities::average_value ( volume_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic );
        out.entropy_derivative_temperature[q] = MaterialUtilities::average_value ( volume_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic );
      }
    }



    template <int dim>
    bool
    VeMel25<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible ();
    }



    template <int dim>
    void
    VeMel25<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("VeMel25");
        {
          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters ( prm );
          Rheology::DiffusionCreep<dim>::declare_parameters ( prm );
          prm.declare_entry ( "Thermal conductivity", "3.96", Patterns::Double ( 0.0 ),
                              "The value of the thermal conductivity $k = \\kappa c_p \\alpha$, where $\\kappa$ is thermal diffusivity. Units: \\si{\\watt\\per\\meter\\per\\kelvin}." );
          prm.declare_entry ( "Minimum viscosity", "1.0e18", Patterns::Double ( 0.0 ),
                              "The value of minimum viscosity $\\eta_{min}$. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ( "Maximum viscosity", "1.0e25", Patterns::Double ( 0.0 ),
                              "The value of maximum viscosity $\\eta_{max}$. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ( "Volatile partition coefficient", "1.0", Patterns::Double ( 0.0 ),
                              "The volatile partition coefficient. Dimensionless, no units." );
          const std::string available_melting_models = "none|McKenzie 1988|Katz 2003";
          prm.declare_entry ( "Melting model", "McKenzie 1988", Patterns::Selection ( available_melting_models ),
                              "Melting model that is to be used. Models available: " + available_melting_models + "." );
          prm.declare_entry ( "Fertile Cpx wt%", "15", Patterns::Double ( 0.0, 100.0 ),
                              "The weigth percent of Clinopyroxene in a fertile mantle for Katz 2003 melting model. Value is expected to be between 0 and 100." );
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    VeMel25<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("VeMel25");
        {
          equation_of_state.initialize_simulator ( this->get_simulator() );
          diffusion_creep.initialize_simulator ( this->get_simulator() );
          equation_of_state.parse_parameters ( prm );
          diffusion_creep.parse_parameters ( prm );
          
          thermal_conductivity  = prm.get_double ( "Thermal conductivity" );
          viscosity_minimum     = prm.get_double ( "Minimum viscosity" );
          viscosity_maximum     = prm.get_double ( "Maximum viscosity" );
          volatile_partition    = prm.get_double ( "Volatile partition coefficient" );

          Cpx                   = prm.get_double ( "Fertile Cpx wt%" ) / 100.0;
          std::string MMS       = prm.get ( "Melting model" );
          if ( MMS == "none" )
            melting_model = MeltingModel::none;
          else if ( MMS == "McKenzie 1988" )
            melting_model = MeltingModel::McKenzie1988;
          else if ( MMS == "Katz 2003" )
            melting_model = MeltingModel::Katz2003;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      this->model_dependence.viscosity = NonlinearDependence::pressure | NonlinearDependence::temperature | NonlinearDependence::compositional_fields | NonlinearDependence::strain_rate;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(VeMel25,
                                   "VeMel25",
                                   "")
  }
}
