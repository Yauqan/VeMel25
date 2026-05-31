
// Melt extraction postprocessor for VeMel25 model

#pragma once

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <algorithm>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim> class SurfaceGrid {
      public:
        inline void set_size ( unsigned int nx ) requires (dim == 2) { sizes[0] = nx; values.resize ( nx ); }
        inline void set_size ( unsigned int nx, unsigned int ny ) requires (dim == 3) { sizes[0] = nx; sizes[1] = ny; values.resize ( nx*ny ); }
        inline void set_extents ( double wx ) requires (dim == 2) { extents[0] = wx; celldata[0] = wx/sizes[0]; celldata[1] = celldata[0]; }
        inline void set_extents ( double wx, double wy ) requires (dim == 3) { extents[0] = wx; extents[1] = wy; celldata[0] = wx/sizes[0]; celldata[1] = wy/sizes[1]; celldata[2] = celldata[0]*celldata[1]; }

        inline double & operator () ( unsigned int ix ) requires (dim == 2) { return values[ix]; }
        inline double & operator () ( unsigned int ix, unsigned int iy ) requires (dim == 3) { return values[ix + sizes[0]*iy]; }
        inline double & operator () ( double x ) requires (dim == 2) {
          unsigned int ix = std::clamp(static_cast<int>(x / celldata[0]), 0, static_cast<int>(sizes[0]-1));
          return values[ix];
        }
        inline double & operator () ( double x, double y ) requires (dim == 3) {
          unsigned int ix = std::clamp(static_cast<int>(x / celldata[0]), 0, static_cast<int>(sizes[0]-1));
          unsigned int iy = std::clamp(static_cast<int>(y / celldata[1]), 0, static_cast<int>(sizes[1]-1));
          return values[ix + sizes[0]*iy];
        }
        inline unsigned int nx () const requires (dim == 2 || dim == 3) { return sizes[0]; }
        inline unsigned int ny () const requires (dim == 3) { return sizes[1]; }
        inline double dx () const requires ( dim == 2 || dim == 3 ) { return celldata[0]; }
        inline double dy () const requires ( dim == 3 ) { return celldata[1]; }
        inline double dS () const { return celldata [ (dim == 2) ? 1 : 2 ]; }
        inline double & operator [] ( size_t indx ) { return values[indx]; }
        inline size_t size () const { return values.size(); }

        inline void operator = ( double val ) { std::fill ( values.begin(), values.end(), val ); }
        inline void operator += ( double val ) { std::transform ( values.begin(), values.end(), values.begin(), [val](double v) { return v + val; } ); }
        inline void operator -= ( double val ) { std::transform ( values.begin(), values.end(), values.begin(), [val](double v) { return v - val; } ); }
        inline void operator *= ( double val ) { std::transform ( values.begin(), values.end(), values.begin(), [val](double v) { return v * val; } ); }
        inline void operator /= ( double val ) { std::transform ( values.begin(), values.end(), values.begin(), [val](double v) { return v / val; } ); }

        template <class Archive> inline void serialize ( Archive & ar, const unsigned int ) { 
          if ( dim == 2 ) {
            ar & sizes[0]
            & extents[0]
            & celldata[0] & celldata[1];
            for ( size_t k = 0; k < values.size(); ++k )
              ar & values[k];
          }
          else if ( dim == 3 ) {
            ar & sizes[0] & sizes[1]
            & extents[0] & extents[1]
            & celldata[0] & celldata[1] & celldata[2];
            for ( size_t k = 0; k < values.size(); ++k )
              ar & values[k];
          }
        }
      private:
        std::array<unsigned int, dim-1> sizes;
        std::array<double, dim-1> extents;
        std::array<double, dim> celldata;   // individual cell sizes, last one is cell's surface area
        std::vector<double> values;
    };

    template <int dim>
    class VeMel25 : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string> execute ( TableHandler &statistics ) override;
        void parse_parameters ( ParameterHandler & prm ) override;
        static void declare_parameters ( ParameterHandler & prm );

        template <class Archive> void serialize ( Archive & ar, const unsigned int version );
        void save ( std::map<std::string, std::string> & status_strings ) const override;
        void load ( const std::map<std::string, std::string> & status_strings ) override;

      private:
        SurfaceGrid<dim> eruptions;
        double rho_lit;
        double rho_crust;

        void set_last_output_time ( const double current_time );

        unsigned int maximum_timesteps_between_outputs;
        double output_interval;
        unsigned int last_output_timestep = numbers::invalid_unsigned_int;
        double last_output_time = std::numeric_limits<double>::quiet_NaN();
        int output_index;
    };
  }
}


