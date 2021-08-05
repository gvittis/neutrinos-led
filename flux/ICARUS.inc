nuflux(#ICARUSplus)<			/* NAME can be replaced, e.g. FLUXplus */
	@flux_file ="./flux/ICARUSplus.dat"	/* loads the flux file */
	@time =	3.       		/* sets the runtime in units years */
  	@power = 2.2                   /* e20 POT per year */
        @norm = 1.61e-3		       /* sets a normalization factor */
>

nuflux(#ICARUSminus)<			/* NAME can be replaced, e.g. FLUXplus */
	@flux_file ="./flux/ICARUSminus.dat"	/* loads the flux file */
	@time =	3.		/* sets the runtime in units years */
  	@power = 2.2            /* e20 POT per year */
        @norm =	1.61e-3		/* sets a normalization factor */
>

/* I can put arbitrary units and renormalize after the running. */
