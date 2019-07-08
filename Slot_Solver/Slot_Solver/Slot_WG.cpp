#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the slot class
// R. Sheehan 5 - 7 - 2019

slot::slot()
{
	// Default constructor
	params_defined = false; 
	lambda = k = n_sl = n_cl = n_h = w_sl = w_h = 0.0; 
	k_nh_sqr = k_nsl_sqr = k_ncl_sqr = nh_nsl_sqr = nh_ncl_sqr = 0.0; 
	k_inv = nh_sqr_inv = nsl_sqr_inv = ncl_sqr_inv = 0.0; 
}

slot::~slot()
{
	// Deconstructor
	params_defined = false; 
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

	try {
		bool c1 = i > -1 ? true : false; 
		bool c2 = static_cast<int>(beta.size()) > 0 ? true : false; 
		bool c3 = i < static_cast<int>(beta.size()) ? true : false; 
		bool c10 = c1 && c2 && c3 && params_defined; 

		if (c10) {
			double v = k_nh_sqr - template_funcs::DSQR(beta[i]);

			return v > 0.0 ? sqrt(v) : 0.0; 
		}
		else {
			std::string reason;
			reason = "Error: double slot::kah(int i)\n";
			if (!c1 || !c3) reason += "index value: " + template_funcs::toString(i) + " not valid\n";
			if (!c2) reason += "propagation constants not computed\n";
			if (!params_defined) reason += "No parameters defined for the object\n"; 
			throw std::invalid_argument(reason);

			return 0.0; 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slot::gsl(int i)
{
	// transverse wavenumber in slot region

	try {
		bool c1 = i > -1 ? true : false;
		bool c2 = static_cast<int>(beta.size()) > 0 ? true : false;
		bool c3 = i < static_cast<int>(beta.size()) ? true : false;
		bool c10 = c1 && c2 && c3 && params_defined;

		if (c10) {
			double v = template_funcs::DSQR(beta[i]) - k_nsl_sqr;

			return v > 0.0 ? sqrt(v) : 0.0;
		}
		else {
			std::string reason;
			reason = "Error: double slot::gsl(int i)\n";
			if (!c1 || !c3) reason += "index value: " + template_funcs::toString(i) + " not valid\n";
			if (!c2) reason += "propagation constants not computed\n";
			if (!params_defined) reason += "No parameters defined for the object\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slot::gcl(int i)
{
	// transverse wavenumber in cladding region

	try {
		bool c1 = i > -1 ? true : false;
		bool c2 = static_cast<int>(beta.size()) > 0 ? true : false;
		bool c3 = i < static_cast<int>(beta.size()) ? true : false;
		bool c10 = c1 && c2 && c3 && params_defined;

		if (c10) {
			double v = template_funcs::DSQR(beta[i]) - k_ncl_sqr;

			return v > 0.0 ? sqrt(v) : 0.0;
		}
		else {
			std::string reason;
			reason = "Error: double slot::gcl(int i)\n";
			if (!c1 || !c3) reason += "index value: " + template_funcs::toString(i) + " not valid\n";
			if (!c2) reason += "propagation constants not computed\n";
			if (!params_defined) reason += "No parameters defined for the object\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slot::prop_const(int i)
{
	// return i^{th} computed propagation constant

	try {
		bool c1 = i > -1 ? true : false;
		bool c2 = static_cast<int>(beta.size()) > 0 ? true : false;
		bool c3 = i < static_cast<int>(beta.size()) ? true : false;
		bool c10 = c1 && c2 && c3 && params_defined;

		if (c10) {
			return beta[i];
		}
		else {
			std::string reason;
			reason = "Error: double slot::prop_const(int i)\n";
			if (!c1 || !c3) reason += "index value: " + template_funcs::toString(i) + " not valid\n";
			if (!c2) reason += "propagation constants not computed\n";
			if (!params_defined) reason += "No parameters defined for the object\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definition of the member functions for the slot_neff class
// R. Sheehan 5 - 7 - 2019

slot_neff::slot_neff()
{
	// Default Constructor
}

slot_neff::slot_neff(double wavelength, double slot_width, double slab_width, double slot_index, double slab_index, double cladding_index)
{
	// Primary Constructor

	set_params(wavelength, slot_width, slab_width, slot_index, slab_index, cladding_index);
}

void slot_neff::set_params(double wavelength, double slot_width, double slab_width, double slot_index, double slab_index, double cladding_index)
{
	try {
	
		bool c1 = wavelength > 0.0 ? true : false; 
		bool c2 = slot_width > 0.0 ? true : false; 
		bool c3 = slab_width > 0.0 ? true : false; 
		bool c4 = slot_index > 0.0 && slot_index < slab_index ? true : false; 
		bool c5 = slab_index > slot_index ? true : false; 
		bool c6 = cladding_index > 0.0 && cladding_index < slab_index ? true : false; 
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6;

		if (c10) {
			
			lambda = wavelength; // wavelength of light in the slot, in units of um
			k = Two_PI / lambda; // wavenumber of light in the slot, in units of um^{-1}
			k_inv = 1.0 / k; // 1 / k_{0}

			w_sl = slot_width; // width of slot region in units of um
			w_h = slab_width; // width of high-index slab region in units of um

			n_sl = slot_index; // RI of material in the slot
			n_h = slab_index; // RI of high index slab region
			n_cl = cladding_index; // RI of material in the cladding, may be the same as that in the slot

			nh_sqr_inv = 1.0 / template_funcs::DSQR(n_h); // 1 / n_{h}^{2} 
			nsl_sqr_inv = 1.0 / template_funcs::DSQR(n_sl); // 1 / n_{sl}^{2}
			ncl_sqr_inv = 1.0 / template_funcs::DSQR(n_cl); // 1 / n_{cl}^{2}

			nh_nsl_sqr = template_funcs::DSQR(n_h / n_sl); // ratio (n_{h} / n_{sl})^{2}
			nh_ncl_sqr = template_funcs::DSQR(n_h / n_cl); // ratio (n_{h} / n_{cl})^{2}

			k_nh_sqr = template_funcs::DSQR(k*n_h); // k_{0}^{2} n_{h}^{2}
			k_nsl_sqr = template_funcs::DSQR(k*n_sl); // k_{0}^{2} n_{sl}^{2}
			k_ncl_sqr = template_funcs::DSQR(k*n_cl); // k_{0}^{2} n_{cl}^{2}

			params_defined = true; 
		}
		else {
			std::string reason; 
			reason = "Error: void slot_neff::set_params(double wavelength, double slot_width, double slab_width, double slot_index, double slab_index, double cladding_index)\n";
			if (!c1) reason += "wavelength: " + template_funcs::toString(wavelength, 2) + " is not valid input\n"; 
			if (!c2) reason += "slot_width: " + template_funcs::toString(slot_width, 2) + " is not valid input\n";
			if (!c3) reason += "slab_width: " + template_funcs::toString(slab_width, 2) + " is not valid input\n";
			if (!c4) reason += "slot_index: " + template_funcs::toString(slot_index, 2) + " is not valid input\n";
			if (!c5) reason += "slab_index: " + template_funcs::toString(slab_index, 2) + " is not valid input\n";
			if (!c6) reason += "cladding_index: " + template_funcs::toString(cladding_index, 2) + " is not valid input\n";
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slot_neff::neff_search()
{
	try {

		if (params_defined) {

		}
		else {
			std::string reason;
			reason = "Error: void slot_neff::void slot_neff::neff_search()\n";
			reason += "parameters not defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slot_neff::eigenequation(double x)
{
	// eigenequation whose roots define the propagation constant of the slot waveguide

	try {

		if (params_defined) {
			
			double xsqr, kah, gsl, gcl, phi, arg1, arg2, t1, t2, t3, tmp; 

			xsqr = template_funcs::DSQR(x); 
			
			tmp = k_nh_sqr - xsqr; 
			kah = tmp > 0 ? sqrt(tmp) : 0.0;

			tmp = xsqr - k_nsl_sqr; 
			gsl = tmp > 0 ? sqrt(tmp) : 0.0;

			tmp = xsqr - k_ncl_sqr; 
			gcl = tmp > 0.0 ? sqrt(tmp) : 0.0;

			if (kah > 0.0) {
				t2 = (gcl / kah); t3 = (gsl / kah); 
			}
			else {
				t2 = t3 = 0.0; 
			}
			
			//t2 = kah > 0.0 ? (gcl / kah) : 0.0; 
			//t3 = kah > 0.0 ? (gsl / kah) : 0.0; 
			
			phi = atan(nh_ncl_sqr*t2);
			
			arg1 = kah * w_sl - phi; 
			
			arg2 = 0.5*gsl*w_sl; 
			 
			t1 = nh_nsl_sqr * t3; 

			return tan(arg1) - t1 * tanh(arg2); 
		}
		else {
			std::string reason;
			reason = "Error: void slot_neff::double slot_neff::eigenequation()\n";
			reason += "parameters not defined\n";
			throw std::invalid_argument(reason);

			return 0.0; 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slot_neff::phi()
{
	try {

		if (params_defined) {

		}
		else {
			std::string reason;
			reason = "Error: void slot_neff::double slot_neff::phi()\n";
			reason += "parameters not defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}

	return 0.0; 
}

double slot_neff::zbrent()
{
	try {	

		if (params_defined) {

		}
		else {
			std::string reason;
			reason = "Error: void slot_neff::double slot_neff::zbrent()\n";
			reason += "parameters not defined\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}

	return 0.0;
}