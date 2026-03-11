
// Melt extraction postprocessor for VeMel25 model

#pragma once

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim> class SurfaceGrid {
      public:
        SurfaceGrid ( unsigned int nx ) requires (dim == 2) : sizes({nx}) {}
        SurfaceGrid ( unsigned int nx, unsigned int ny ) requires (dim == 3) : sizes({nx, ny}) {}
        double & operator () ( unsigned int ix ) requires (dim == 2) { return values[ix]; }
        double & operator () ( unsigned int ix, unsigned int iy ) requires (dim == 3) { return values[ix + sizes[0]*iy]; }
        void operator = ( double val ) { std::fill ( values.begin(), values.end(), val ); }
        void operator += ( double val ) { std::transform ( values.begin(), values.end(), [](double & v) { return v + val; } ); }
        void operator -= ( double val ) { std::transform ( values.begin(), values.end(), [](double & v) { return v - val; } ); }
        void operator *= ( double val ) { std::transform ( values.begin(), values.end(), [](double & v) { return v * val; } ); }
        void operator /= ( double val ) { std::transform ( values.begin(), values.end(), [](double & v) { return v / val; } ); }
        unsigned int nx () requires (dim == 2 || dim == 3) { return sizes[0]; }
        unsigned int ny () requires (dim == 3) { return sizes[1]; }
      private:
        std::array<unsigned int, dim-1> sizes;
        std::vector<double> values;
    };
    
    template <int dim>
    class VeMel25 : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string> execute ( TableHandler &statistics ) override;
        std::list<std::string> required_other_postprocessors() const override;
        void parse_parameters ( ParameterHandler & prm ) override;
        static void declare_parameters ( ParameterHandler & prm );

        template <class Archive> void serialize ( Archive & ar, const unsigned int version );
        void save ( std::map<std::string, std::string> & status_strings ) const override;
        void load ( const std::map<std::string, std::string> & status_strings ) override;

      private:
        SurfaceGrid<dim> eruptions;

        void set_last_output_time ( const double current_time );

        unsigned int maximum_timesteps_between_outputs;
        double output_interval;
        unsigned int last_output_timestep = numbers::invalid_unsigned_int;
        double last_output_time = std::numeric_limits<double>::quiet_NaN();
        int output_index;
    };
  }
}


