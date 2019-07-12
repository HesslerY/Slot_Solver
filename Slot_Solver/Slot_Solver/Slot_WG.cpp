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
	coeffs_defined = beta_defined = params_defined = false;
	lambda = k = n_sl = n_cl = n_h = w_sl = w_h = a = b = 0.0; 
	k_nh_sqr = k_nsl_sqr = k_ncl_sqr = nh_nsl_sqr = nh_ncl_sqr = 0.0; 
	beta_high = beta_low = k_inv = nh_sqr_inv = nsl_sqr_inv = ncl_sqr_inv = 0.0; 
	beta = neff = kah = gsl = gcl = ampl = 0.0; 
	c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = c9 = c10 = 0.0; 
}

slot::~slot()
{
	// Deconstructor
	coeffs_defined = beta_defined = params_defined = false;
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
			
			lambda = wavelength; // wavelength of light in the slot, in units of nm
			k = Two_PI / lambda; // wavenumber of light in the slot, in units of nm^{-1}
			k_inv = 1.0 / k; // 1 / k_{0}

			w_sl = slot_width; // width of slot region in units of nm
			a = 0.5*w_sl; // half slot width in units of nm
			w_h = slab_width; // width of high-index slab region in units of nm
			b = a + w_h; // location of high-index slab edge in units of nm

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
			int r = static_cast<int>(sub_intervals.size()) - 1; 

			beta = zbrent(sub_intervals[r].get_x_lower(), sub_intervals[r].get_x_upper(), 1.0e-6); 

			neff = beta / k; 

			beta_defined = true; 

			if (loud) std::cout << "beta = " << std::setprecision(6) << beta << ", n_{eff} = " << neff << "\n";
			
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
	// compute equation (3) of the source paper

	try {

		if (params_defined) {
			
			// use local values for kah, gsl, gcl to avoid confusion when computing field profile
			double xsqr, kkah, ggsl, ggcl, phi, arg1, arg2, t1, t2, t3, tmp; 

			xsqr = template_funcs::DSQR(x); 
			
			tmp = k_nh_sqr - xsqr; 
			kkah = tmp > 0 ? sqrt(tmp) : 0.0;

			tmp = xsqr - k_nsl_sqr; 
			ggsl = tmp > 0 ? sqrt(tmp) : 0.0;

			tmp = xsqr - k_ncl_sqr; 
			ggcl = tmp > 0.0 ? sqrt(tmp) : 0.0;

			if (kkah > 0.0) {
				t2 = (ggcl / kkah); 
				t3 = (ggsl / kkah); 
			}
			else {
				t2 = t3 = 0.0; 
			}
		
			phi = atan(nh_ncl_sqr*t2);
			
			arg1 = kkah * w_h - phi; 
			
			arg2 = ggsl * a;
			 
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

// Definitions for the members of slot_mode class
// R. Sheehan 12 - 7 - 2019

slot_mode::slot_mode()
{
	// default constructor
}

slot_mode::slot_mode(double wavelength, double slot_width, double slab_width, double slot_index, double slab_index, double cladding_index)
{
	// primary constructor
	set_params(wavelength, slot_width, slab_width, slot_index, slab_index, cladding_index);
}

void slot_mode::set_mode_params()
{
	// Compute various parameters required to compute the slot mode profile
	// R. Sheehan 12 - 7 - 2019

	try {
		if (beta_defined) {
			double xsqr, tmp; 

			xsqr = template_funcs::DSQR(beta);

			tmp = k_nh_sqr - xsqr;
			kah = tmp > 0 ? sqrt(tmp) : 0.0; // transverse wavenumber in high index slab region

			tmp = xsqr - k_nsl_sqr;
			gsl = tmp > 0 ? sqrt(tmp) : 0.0; // field decay coefficient in cladding

			tmp = xsqr - k_ncl_sqr;
			gcl = tmp > 0.0 ? sqrt(tmp) : 0.0; // field decay coefficient in slot region
		
			gsl_kah = gsl / kah; // ratio g_{sl} / k_{ah}

			c1 = gsl * a; // argument for the cosh and sinh terms in the profile

			c2 = cosh(c1); c7 = sinh(c1); // sinh and cosh terms in the profile

			c3 = gsl_kah * nsl_sqr_inv; // g_{sl} / (n_{sl}^{2} k_{ah})

			c4 = kah * w_h; // argument for cosine and sine terms in the profile

			c5 = cos(c4); c8 = sin(c4); // cosine and sine terms in the profile

			c6 = nh_nsl_sqr * gsl_kah; // // ratio (n_{h} / n_{sl} )^{2} * ( g_{sl} / k_{ah} )

			c9 = nh_sqr_inv * c2; // ( 1 / n_{h}^{2} ) cosh ( g_{sl} a  )

			c10 = ncl_sqr_inv * ( ( c2 * c5 ) + ( c6 * c7 * c8 ) ); // ( 1 / n_{cl}^{2} ) * various

			coeffs_defined = true; 
		}
		else {
			std::string reason;
			reason = "Error: void slot_mode::set_mode_params()\n";
			reason += "Cannot proceed with calculation as propagation constant is not defined\n"; 
			throw std::invalid_argument(reason);			
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slot_mode::Ex(double x)
{
	// return the value of the slot waveguide field
	// working in nm scale for field calculation
	// R. Sheehan 12 - 7 - 2019

	try {
	
		if (coeffs_defined) {
			
			double arg, t = fabs(x); 

			if (t < a || fabs(t-a) < 1.0e-9 ) {
				arg = gsl * x; 
				return nsl_sqr_inv * cosh(arg); 
			}
			else if (t > a && t < b) {
				arg = kah * (t - a); 
				return c9 * cos(arg) + c3*sin(arg);
			}
			else if (t > b || fabs(t-b) < 1.0e-9) {
				arg = -1.0*gcl*(t - b); 
				return c10 * exp(arg);
			}
			else {
				return 0.0; 
			}		
		}
		else {
			std::string reason;
			reason = "Error: double slot_mode::Ex(double x)\n";
			reason += "Cannot proceed with calculation as solution coefficients are not defined\n";
			throw std::invalid_argument(reason);
			return 0.0; 
		}

	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}