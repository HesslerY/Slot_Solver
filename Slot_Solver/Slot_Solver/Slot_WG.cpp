#ifndef ATTACH_H
#include "Attach.h"
#endif

// Defintions of the members of interval
// R. Sheehan 8 - 3 - 2013

// interval object
// Constructors
interval::interval()
{
	// Default constructor
	interval_defined = false;
	xlower = xupper = 0.0;
}

interval::interval(double xl, double xu)
{
	// Constructor
	// construct an interval over the range [xl, xu]
	set_xl_xu(xl, xu);
}

//Methods
void interval::set_xl_xu(double xl, double xu)
{
	// set the values of xlower and xupper
	// ensure that xl < xu

	try {

		if (fabs(xl - xu) > 1.0e-12) {

			// xl != xu => interval can be created

			xlower = std::min(xl, xu);

			xupper = std::max(xl, xu);

			length = xupper - xlower; 

			interval_defined = true;
		}
		else {
			// xl == xu throw exception

			std::string reason = "Error: void interval::set_xl_xu(double xl, double xu)\n";
			reason += "Attempting to construct an interval in which xl == xu\n";
			reason += "xl = " + template_funcs::toString(xl, 4) + "\n";
			reason += "xu = " + template_funcs::toString(xl, 4) + "\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definition of the slot class
// R. Sheehan 5 - 7 - 2019

slot::slot()
{
	// Default constructor
	params_defined = false; 
	lambda = k = n_sl = n_cl = n_h = w_sl = w_h = 0.0; 
	k_nh_sqr = k_nsl_sqr = k_ncl_sqr = nh_nsl_sqr = nh_ncl_sqr = 0.0; 
	beta_high = beta_low = k_inv = nh_sqr_inv = nsl_sqr_inv = ncl_sqr_inv = 0.0; 
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

			beta_high = k * n_h; // search space bounds
			beta_low = k * std::min(n_sl, n_cl); 

			k_nh_sqr = template_funcs::DSQR(beta_high); // k_{0}^{2} n_{h}^{2}
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

void slot_neff::neff_search(bool loud)
{
	// Compute slot waveguide fundamental mode effective index
	// Only compute the fundamental mode effective index since the second solution may not be physically valid
	// since it is located at the position of a discontinuity in the eigenequation
	// R. Sheehan 12 - 7 - 2019

	try {

		if (params_defined) {

			bracket_roots(loud); // locate intervals on which roots can be found

			// if no solutions exist program will exit before this point
			double sol;
			int r = static_cast<int>(sub_intervals.size()) - 1; 

			sol = zbrent(sub_intervals[r].get_x_lower(), sub_intervals[r].get_x_upper(), 1.0e-6); 

			beta.push_back(sol); 

			if (loud) std::cout << "beta[" << r + 1 << "] = " << std::setprecision(6) << sol << " , n_{eff} = " << sol / k << "\n";
			
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
				t2 = (gcl / kah); 
				t3 = (gsl / kah); 
			}
			else {
				t2 = t3 = 0.0; 
			}
			
			//t2 = kah > 0.0 ? (gcl / kah) : 0.0; 
			//t3 = kah > 0.0 ? (gsl / kah) : 0.0; 
			
			phi = atan(nh_ncl_sqr*t2);
			
			arg1 = kah * w_h - phi; 
			
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

void slot_neff::print_eigenequation()
{
	// print the eigenequation data to a file for analysis
	// R. Sheehan 8 - 7 - 2019

	try {

		if (params_defined) {

			std::string name = "Slot_WG_Eigenequation.txt"; 

			std::ofstream write(name, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {
				int Nsteps = 151;
				double x, delta_x;

				delta_x = (beta_high - beta_low) / (Nsteps - 1);
				x = beta_low+delta_x;

				int j = 1;
				while (j < Nsteps-1) {
					write << x << " , " << eigenequation(x) << "\n"; 
					x += delta_x;
					j++;
				}

				write.close();
			}
			else {
				write.close();
				std::string reason;
				reason = "Error: void slot_neff::print_eigenequation()\n";
				reason += "Filename: " + name + " is not valid\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void slot_neff::print_eigenequation()\n";
			reason += "parameters not defined\n";
			throw std::invalid_argument(reason);
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

double slot_neff::zbrent(double x1, double x2, double tol)
{
	//Find the root of the function between x1 and x2 using brent's method
	//The root is refined to +/- tol
	//Seemingly one of the best methods to use

	//This will be used to compute the roots of eigeneq_3
	//R. Sheehan 28 - 5 - 2010

	try {	

		bool c1 = fabs(x1 - x2) > 1.0e-9 ? true : false; // cannot have x1 == x2
		bool c2 = tol > 1.0e-16 ? true : false;
		bool c10 = c1 && c2 && params_defined; 

		if (c10) {
			int iter;

			static const int ITMAX = 100;//Used in rtbis, rtflsp, rtsec, zriddr

			double a = std::min(x1, x2), b = std::max(x1, x2), c = std::max(x1, x2), d, e, min1, min2;
			double fc, p, q, r, s, tol1, xm;
			double fa = eigenequation(a), fb = eigenequation(b);

			if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
				std::cerr << "Root must be bracketed in zbrent\n";
			}
			fc = fb;
			for (iter = 1; iter <= ITMAX; iter++) {
				if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (fabs(fc)<fabs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}
				tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
				xm = 0.5*(c - b);
				if (fabs(xm) <= tol1 || fb == 0.0) return b;
				/*if(fabs(xm)<=tol1 || fb==0.0){
				std::cout<<"Brent's Method converged in "<<iter<<" iterations\n";
				return b;
				}*/
				if (fabs(e) >= tol1 && fabs(fa)>fabs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0*xm*s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
						q = (q - 1.0)*(r - 1.0)*(s - 1.0);
					}
					if (p>0.0) q = -q;
					p = fabs(p);
					min1 = 3.0*xm*q - fabs(tol1*q);
					min2 = fabs(e*q);
					if (2.0*p<std::min(min1, min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (fabs(d)>tol1) {
					b += d;
				}
				else {
					b += template_funcs::SIGN(tol1, xm);
				}
				fb = eigenequation(b);
			}
			std::cerr << "Maximum number of iterations exceeded in zbrent\n";
			return 0.0;
		}
		else {
			std::string reason;
			reason = "Error: void slot_neff::double slot_neff::zbrent(double x1, double x2, double tol)\n";
			if (!c1) reason += "Cannot have x1 = x2\nx1 = " + template_funcs::toString(x1) + ", x2 = " + template_funcs::toString(x1) + "\n";
			if (!c2) reason += "Desired tolerance is less than smallest allowable EPS\n";
			if(!params_defined) reason += "parameters not defined\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slot_neff::bracket_roots(bool loud)
{
	// Find the sub-intervals of [a, b] known to contain a root of the equation f(x) = 0
	// Requires that a search_space [a, b] be defined
	// R. Sheehan 8 - 3 - 2013

	try {

		if (params_defined) {

			int nroots = 0, nsub = 151;
			double dx;

			dx = (beta_high - beta_low) / (nsub - 1);
			
			interval search_space; // this will define the interval [a, b]

			// set the search space away from the exact endpoints as there seems to be discontinuitues in the eigenequations 
			// near the actual endpoints
			nsub -= 2; 
			search_space.set_xl_xu(beta_low + dx, beta_high - dx);

			if ( search_space.has_bounds() ) {
				// search_space is properly bounded, search for sub-intervals containing roots can proceed
				// nroots is not known a-priori, each of nsub sub-intervals must be tested to see if contains a root
				// in general user chooses some nsub > nroot, correct value can be found by experimentation

				double fl, fu;
				double xl, xu;
				
				xl = search_space.get_x_lower();
				xu = xl + dx;
				for (int i = 0; i < nsub; i++) {

					// Evaluate the function on the sub-interval
					if (i == 0) {
						fl = template_funcs::Signum( eigenequation(xl) );
					}
					else {
						fl = fu; // minimise the number of function evaluations
					}

					fu = template_funcs::Signum( eigenequation(xu) );

					// Perform the bisection test
					if (fl*fu < 0.0) {
						// the sub-interval contains a root so it is stored
						nroots++;
						sub_intervals.push_back( interval(xl, xu) );
					}

					// update the endpoints of the sub-interval
					xl = xu;
					xu += dx;
				}

				if (nroots == 0) {
					// eigenequation does not have any solutions on the search interval
					std::string reason = "Error: void find_root::bracket_roots()\n";
					reason += "eigenequation has no solutions in [beta_{low}, beta_{high}]\n";
					reason += "root search cannot proceed\n";

					throw std::invalid_argument(reason);
				}

				if (loud) {
					std::cout << "The function contains " << nroots << " roots on the interval [ " << search_space.get_x_lower() << " , " << search_space.get_x_upper() << " ]\n";
					if (nroots > 0) {
						std::cout << "The roots are located in \n";
						for (int i = 0; i < nroots; i++) {
							std::cout << "[ " << sub_intervals[i].get_x_lower() << " , " << sub_intervals[i].get_x_upper() << " ], dx = " << sub_intervals[i].get_length() << "\n";
						}
						std::cout << "\n";
					}
				}
			}
			else {
				// search_space is not bounded properly, throw exception
				std::string reason = "Error: void find_root::bracket_roots()\n";
				reason += "search_space is not properly bounded\n";
				reason += "root search cannot proceed\n";

				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void slot_neff::bracket_roots()\n";
			reason += "parameters not defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}