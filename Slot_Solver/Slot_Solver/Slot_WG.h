#ifndef SLOT_WG_H
#define SLOT_WG_H

// Data type for an interval [xlower, xupper]
// R. Sheehan 8 - 3 - 2013
class interval {
public:
	// Constructor
	interval();
	interval(double xl, double xu);

	// Methods

	void set_xl_xu(double xl, double xu);

	inline bool has_bounds() { return interval_defined; }

	inline double get_x_lower() { return xlower; }
	inline double get_x_upper() { return xupper; }
	inline double get_length() { return length;  }

private:
	bool interval_defined;
	double xlower;
	double xupper;
	double length; 
};

// Implementation of the code needed to compute the solution in a 1D slot waveguide structure
// Formulae are given in the paper "Guiding and confining light in void nanostructure", Opt. Lett., 29 (11), 2004
// I will follow a similar structure to what I've used for my Slab_WG code, i.e. 
// base class containing all parameters, derived class for computing effective indices and separate derived class for computing mode profile
// RI values and wavelength will be input manually to maintain simple interface. 
// R. Sheehan 5 - 7 - 2019

class slot {
public:
	slot();

	~slot(); 

protected:
	inline bool params_status() { return params_defined;  }

	inline double prop_const() { return beta;  } // return computed propagation constant

protected:
	bool params_defined;
	bool beta_defined;
	bool coeffs_defined; 

	// work in nm, to ensure field calculation is less prone to numerical error
	double lambda; // wavelength of light in the slot, in units of nm
	double k; // wavenumber of light in the slot, in units of nm^{-1}
	
	double n_sl; // RI of material in the slot
	double n_cl; // RI of material in the cladding, may be the same as that in the slot
	double n_h; // RI of high index slab region

	double w_sl; // width of slot region in units of nm
	double w_h; // width of high-index slab region in units of nm
	double a; // half slot width in units of um
	double b; // location of high-index slab edge in units of nm

	double k_nh_sqr; // k_{0}^{2} n_{h}^{2}
	double k_nsl_sqr; // k_{0}^{2} n_{sl}^{2}
	double k_ncl_sqr; // k_{0}^{2} n_{cl}^{2}
	
	double nh_nsl_sqr; // ratio (n_{h} / n_{sl})^{2}
	double nh_ncl_sqr; // ratio (n_{h} / n_{cl})^{2}
	
	double k_inv; // 1 / k_{0}

	double nh_sqr_inv; // 1 / n_{h}^{2}
	double nsl_sqr_inv; // 1 / n_{sl}^{2}
	double ncl_sqr_inv; // 1 / n_{cl}^{2}

	double beta_high; // bounds on the search space for the slot effective index 
	double beta_low; 

	double beta; // TM propagation constant, in units of nm^{-1}
	double neff; // TM effective index, unitless

	// parameters associated with calculation of mode profile
	// can only be compute once slot WG beta is known
	double kah; // transverse wavenumber in high index slab region
	double gcl; // field decay coefficient in cladding
	double gsl; // field decay coefficient in slot region
	double ampl; // field amplitude, possibly equal to slot effective index

	double gsl_kah; // ratio g_{sl} / k_{ah}
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11; // various constants needed for mode profile calculation
};

// Use this class to compute the propagation constant of a slot waveguide structure

class slot_neff : protected slot {
public:
	slot_neff(); 
	slot_neff(double wavelength, double slot_width, double slab_width, double slot_index, double slab_index, double cladding_index); 

	void set_params(double wavelength, double slot_width, double slab_width, double slot_index, double slab_index, double cladding_index); 

	void neff_search(bool loud);

	void print_eigenequation();

private:
	double eigenequation(double x); 

	double zbrent(double x1, double x2, double tol);	

	void bracket_roots(bool loud = false); // function for bracketing the roots

private:
	std::vector<interval> sub_intervals; // an array of sub-intervals known to contain a roots
};

class slot_mode : public slot_neff {
public:
	slot_mode(); 

	slot_mode(double wavelength, double slot_width, double slab_width, double slot_index, double slab_index, double cladding_index);

	void output_statistics(std::string &storage_directory);

	void output_mode_profile(int N, double Lx, std::string &storage_directory);

private:
	void set_mode_params(); 

	double Ex(double x); 
};

#endif
