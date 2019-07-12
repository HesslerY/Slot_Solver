#ifndef  ATTACH_H
#include "Attach.h"
#endif // ! ATTACH_H

void testing::example_1()
{
	// implementation of the code needed to test operation of slot solver
	// example taken from paper "Guiding and confining light in void nanostructure", Opt. Lett., 29 (11), 2004

	double ns, nh, wh, ws, wl; 

	wl = 1.55; // wavelength in units of um
	ns = 1.44; // slot index
	nh = 3.48; // slab index
	wh = (180.0 / 1000.0); // slab thickness in units of um
	ws = (25 / 1000.0); // slot thickness in units of um

	slot_neff waveguide; 

	waveguide.set_params(wl, ws, wh, ns, nh, ns); 
	waveguide.neff_search(true); 
}
