nuflux(#LARplus)<			/* NAME can be replaced, e.g. FLUXplus */
	@flux_file ="./flux/LARplus.dat"	/* loads the flux file */
	@time =	3.       		/* sets the runtime in units years */
  	@power = 2.2                   /* e20 POT per year */
        @norm = 5.39e-5		       /* sets a normalization factor */
>

nuflux(#LARminus)<			/* NAME can be replaced, e.g. FLUXplus */
	@flux_file ="./flux/LARminus.dat"	/* loads the flux file */
	@time =	3.		/* sets the runtime in units years */
  	@power = 2.2            /* e20 POT per year */
        @norm =	5.39e-5		/* sets a normalization factor */
>

/* I can put arbitrary units and renormalize after the running. */
