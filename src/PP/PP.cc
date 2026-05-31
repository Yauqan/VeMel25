
#include <aspect/simulator.h>
#include "PP.h"

#include <aspect/postprocess/boundary_pressures.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/mpi.h>
#include <limits>
#include <map>
#include <filesystem>
#include <fstream>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    VeMel25<dim>::execute (TableHandler & statistics)
    {
      const aspect::GeometryModel::Box<dim> * box = dynamic_cast<const aspect::GeometryModel::Box<dim>*> ( &this->get_geometry_model() );
      AssertThrow ( box != nullptr, ExcMessage("The postprocessor only works for Box geometry, sorry!") );

      const unsigned int mantle_depletion_indx = this->introspection().compositional_index_for_name ( "mantle_depletion" );
      const Quadrature<dim> & quadform = this->introspection().quadratures.compositional_fields[mantle_depletion_indx];
      const unsigned int n_q_points = quadform.size();
      MaterialModel::MaterialModelInputs<dim> in(n_q_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(n_q_points, this->n_compositional_fields());
      FEValues<dim> fe_values ( this->get_mapping(), this->get_fe(), quadform, update_values | update_gradients | update_quadrature_points | update_JxW_values );

      std::vector<double> mantle_depletion ( n_q_points );

      double partial_eruptions = 0.0;
      for ( const auto & cell : this->get_dof_handler().active_cell_iterators() )
        if ( cell->is_locally_owned() ) {
          fe_values.reinit(cell);
          in.reinit(fe_values, cell, this->introspection(), this->get_solution());
          this->get_material_model().evaluate(in, out);

          for ( unsigned int q = 0; q < n_q_points; ++q ) {
            const double depeltion_change = out.reaction_terms[q][mantle_depletion_indx];
            const double erupt = depeltion_change * fe_values.JxW(q) * rho_lit / rho_crust;
            partial_eruptions += erupt;
            if constexpr ( dim == 2 )
              eruptions ( fe_values.quadrature_point(q)[0] ) += erupt;
            else if ( dim == 3 )
              eruptions ( fe_values.quadrature_point(q)[0], fe_values.quadrature_point(q)[1] ) += erupt;
          }
        }
      
      for ( size_t i = 0; i < eruptions.size(); ++i )
        eruptions[i] = Utilities::MPI::sum ( eruptions[i], this->get_mpi_communicator() );
      partial_eruptions = Utilities::MPI::sum ( partial_eruptions, this->get_mpi_communicator() );
      statistics.add_value ( "Solidified eruption volume", partial_eruptions );
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
          last_output_timestep = this->get_timestep_number();
        }
      if ((this->get_time() < last_output_time + output_interval)
          && (this->get_timestep_number() < last_output_timestep + maximum_timesteps_between_outputs)
          && (this->get_timestep_number() != 0))
        return {"Computing crustal thickness increment", ""};

      double total_eruptions = 0.0;
      
      if ( Utilities::MPI::this_mpi_process ( this->get_mpi_communicator() ) == 0 ) {
        std::filesystem::path addrout ( this->get_output_directory() + "Eruptions/crustal_topography_increment." + Utilities::int_to_string(output_index, 5) );
        if ( addrout.has_parent_path() )
          std::filesystem::create_directories ( addrout.parent_path() );
        std::ofstream out ( addrout );
        if constexpr ( dim == 2 ) {
          out << "# x crustal_thickness_increment\n";
          for ( unsigned int ix = 0; ix < eruptions.nx(); ++ix ) {
            const double x = eruptions.dx() * ( 0.5 + ix );
            out << x << ' ' << eruptions ( ix ) / eruptions.dS() << std::endl; 
            total_eruptions += eruptions ( ix );
          }
        }
        else if ( dim == 3 ) {
          out << "# x y crustal_thickness_increment\n";
          for ( unsigned int ix = 0; ix < eruptions.nx(); ++ix )
            for ( unsigned int iy = 0; iy < eruptions.ny(); ++iy ) {
              const double x = eruptions.dx() * ( 0.5 + ix );
              const double y = eruptions.dy() * ( 0.5 + iy );
              out << x << ' ' << y << ' ' << eruptions ( ix, iy ) / eruptions.dS() << std::endl;
              total_eruptions += eruptions ( ix, iy );
            }
        }
      }

      total_eruptions = Utilities::MPI::sum ( total_eruptions, this->get_mpi_communicator() );
      eruptions = 0.0;
      set_last_output_time ( this->get_time() );
      last_output_timestep = this->get_timestep_number();
      ++output_index;

      return std::pair<std::string,std::string>("Outputting crustal thickness increment", "total erupted volume: " + std::to_string(total_eruptions) + " m^3");
    }

    template <int dim> void VeMel25<dim>::declare_parameters ( ParameterHandler & prm ) {
      prm.enter_subsection ( "Postprocess" );
        prm.enter_subsection ( "VeMel25" );
          prm.declare_entry ( "Gridpoints", "1",
                              Patterns::List(Patterns::Integer(0)),
                              "The number of gridpoints in individual directions, dim-1 values expected");
          prm.declare_entry ( "Mantle density", "3200.0",
                              Patterns::Double(0.0),
                              "The density of the melted material." );
          prm.declare_entry ( "Crustal density", "2900.0",
                              Patterns::Double(0.0),
                              "The density of the solidified crust." );
        prm.leave_subsection ();
      prm.leave_subsection ();
    }

    template <int dim> void VeMel25<dim>::parse_parameters ( ParameterHandler & prm ) {
      prm.enter_subsection ( "Postprocess" );
        prm.enter_subsection ( "VeMel25" );
          std::vector<std::string> gridpoints_field_names = (dim == 2) ? std::vector<std::string>{"x"} : std::vector<std::string>{"x", "y"};
          Utilities::MapParsing::Options options ( gridpoints_field_names, "Gridpoints" );
          options.list_of_allowed_keys = gridpoints_field_names;
          options.allow_multiple_values_per_key = false;
          std::vector<double> tmp = Utilities::MapParsing::parse_map_to_double_array ( prm.get("Gridpoints"), options );
          if constexpr ( dim == 2 )
            eruptions.set_size ( int(tmp[0]) );
          else if ( dim == 3 )
            eruptions.set_size ( int(tmp[0]), int(tmp[1]) );
          rho_lit = prm.get_double ( "Mantle density" );
          rho_crust = prm.get_double ( "Crustal density" );
        prm.leave_subsection ();
        prm.enter_subsection ( "Visualization" );
          output_interval = prm.get_double ( "Time between graphical output" );
          if ( this->convert_output_to_years() )
            output_interval *= year_in_seconds;
          maximum_timesteps_between_outputs = prm.get_integer ( "Time steps between graphical output" );
          if (output_interval > 0.0)
            AssertThrow ( this->get_parameters().run_postprocessors_on_nonlinear_iterations == false,
                          ExcMessage ( "Postprocessing nonlinear iterations is only supported if every time "
                                      "step is visualized, or in other words, if the 'Time between graphical "
                                      "output' in the Visualization postprocessor is set to zero." ) );
        prm.leave_subsection ();
      prm.leave_subsection ();
      prm.enter_subsection ( "Geometry model" );
        prm.enter_subsection ( "Box" );
          if constexpr ( dim == 2 )
            eruptions.set_extents ( prm.get_double ( "X extent" ) );
          else if ( dim == 3 )
            eruptions.set_extents ( prm.get_double ( "X extent" ), prm.get_double ( "Y extent" ) );
        prm.leave_subsection ();
      prm.leave_subsection ();
    }
  
    template <int dim> void VeMel25<dim>::set_last_output_time ( const double current_time )
    {
      if ( output_interval > 0 )
      {
        const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
        last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
      }
    }

    template <int dim> template <class Archive> void VeMel25<dim>::serialize ( Archive & ar, const unsigned int version )
    {
      eruptions.serialize ( ar, version );
      ar & last_output_time
      & last_output_timestep
      & output_index
      ;
    }

    template <int dim> void VeMel25<dim>::save ( std::map<std::string, std::string> & status_strings ) const
    {
      std::ostringstream os;
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }
      status_strings["VeMel25Postprocessor"] = os.str();
    }

    template <int dim> void VeMel25<dim>::load ( const std::map<std::string, std::string> & status_strings )
    {
      if (status_strings.find("VeMel25Postprocessor") != status_strings.end())
      {
        std::istringstream is (status_strings.find("VeMel25Postprocessor")->second);
        aspect::iarchive ia (is);
        ia >> (*this);
      }
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VeMel25,
                                  "VeMel25",
                                  "Postprocessor that outputs the map of heights of produced additional crustal material for every visualization timestep.")
  }
}
