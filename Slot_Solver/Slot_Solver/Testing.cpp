#ifndef  ATTACH_H
#include "Attach.h"
#endif // ! ATTACH_H

void testing::example_1()
{
	// implementation of the code needed to test operation of slot solver
	// example taken from paper "Guiding and confining light in void nanostructure", Opt. Lett., 29 (11), 2004

	double ns, nh, wh, ws, wl; 

	wl = 1550; // wavelength in units of nm
	ns = 1.44; // slot index
	nh = 3.48; // slab index
	//wh = (180.0 / 1000.0); // slab thickness in units of nm
	//ws = (25 / 1000.0); // slot thickness in units of nm
	wh = 180.0; // slab thickness in units of nm
	ws = 25; // slot thickness in units of nm

	slot_neff waveguide; 

	waveguide.set_params(wl, ws, wh, ns, nh, ns); 
	waveguide.neff_search(true); 
}

void testing::example_2()
{
	// implementation of the code needed to test operation of slot solver
	// example taken from paper "Guiding and confining light in void nanostructure", Opt. Lett., 29 (11), 2004

	double ns, nh, ncl, wh, ws, wl;

	wl = 1550; // wavelength in units of nm
	ns = 1.44; // slot index
	ncl = 1.0; // slot index
	nh = 3.48; // slab index
	//wh = (180.0 / 1000.0); // slab thickness in units of nm
	//ws = (25 / 1000.0); // slot thickness in units of nm
	wh = 180.0; // slab thickness in units of nm
	ws = 75; // slot thickness in units of nm

	slot_mode waveguide;

	waveguide.set_params(wl, ws, wh, ns, nh, ncl);

	waveguide.print_eigenequation(); 

	waveguide.neff_search(true); 

	int nn = 501; 
	double lx = 1000.0; // length of solution domain in units of nm
	std::string storage_dir = empty_str; 

	waveguide.output_mode_profile(nn, lx, storage_dir);

	waveguide.output_statistics(storage_dir);
}
