#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the slot class
// R. Sheehan 5 - 7 - 2019

slot::slot()
{
	// Default constructor
	wavelength = k = n_sl = n_cl = n_h = w_sl = w_h = 0.0; 
	k_nh_sqr = k_nsl_sqr = k_ncl_sqr = nh_nsl_sqr = nh_ncl_sqr = 0.0; 
	k_inv = nh_sqr_inv = nsl_sqr_inv = ncl_sqr_inv = 0.0; 
}

slot::~slot()
{
	// Deconstructor
	beta.clear(); 
}

int slot::nbeta()
{
	// return the number of computed propagation constants

	return static_cast<int>(beta.size()); 
}

double slot::kah(int i)
{
	// transverse wavenumber in high-index slab region

	return 0.0; 
}

double slot::gsl(int i)
{
	// transverse wavenumber in slot region

	return 0.0;
}

double slot::gcl(int i)
{
	// transverse wavenumber in cladding region

	return 0.0;
}

double slot::prop_const(int i)
{
	// return i^{th} computed propagation constant

	return 0.0; 
}

// Definition of the member functions for the slot_neff class
// R. Sheehan 5 - 7 - 2019

slot_neff::slot_neff()
{
	
}

void slot_neff::set_params()
{

}

void slot_neff::neff_search()
{
	
}

double slot_neff::eigenequation()
{
	return 0.0; 
}

double slot_neff::phi()
{
	return 0.0; 
}

double slot_neff::zbrent()
{
	return 0.0;
}