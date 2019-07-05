#ifndef SLOT_WG_H
#define SLOT_WG_H

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
	int nbeta(); // return no. of modes in waveguide

	double kah(int i); // transverse wavenumber in high-index slab region
	double gsl(int i); // transverse wavenumber in slot region
	double gcl(int i); // transverse wavenumber in cladding region
	double prop_const(int i); // return i^{th} computed propagation constant

protected:
	double wavelength; // wavelength of light in the slot, in units of um
	double k; // wavenumber of light in the slot, in units of um^{-1}
	double n_sl; // RI of material in the slot
	double n_cl; // RI of material in the cladding, may be the same as that in the slot
	double n_h; // RI of high index slab region
	double w_sl; // width of slot region in units of um
	double w_h; // width of high-index slab region in units of um

	double k_nh_sqr; // k_{0}^{2} n_{h}^{2}
	double k_nsl_sqr; // k_{0}^{2} n_{sl}^{2}
	double k_ncl_sqr; // k_{0}^{2} n_{cl}^{2}
	double nh_nsl_sqr; // ratio (n_{h} / n_{sl})^{2}
	double nh_ncl_sqr; // ratio (n_{h} / n_{cl})^{2}
	double k_inv; // 1 / k_{0}
	double nh_sqr_inv; // 1 / n_{h}^{2}
	double nsl_sqr_inv; // 1 / n_{sl}^{2}
	double ncl_sqr_inv; // 1 / n_{cl}^{2}

	std::vector<double> beta; // TM propagation constants, in units of um^{-1}
};

// Use this class to compute the propagation constant of a slot waveguide structure

class slot_neff : protected slot {
public:
	slot_neff(); 


	void set_params(); 

	void neff_search(); 

private:
	double eigenequation(); 
	double phi(); 

	double zbrent(); 
};

class slot_mode : public slot_neff {};

#endif
