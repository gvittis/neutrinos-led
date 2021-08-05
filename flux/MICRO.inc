nuflux(#MICROplus)<			/* NAME can be replaced, e.g. FLUXplus */
	@flux_file ="./flux/MICROplus.dat"	/* loads the flux file */
	@time =	6.       		/* sets the runtime in units years */
  	@power = 2.2                   /* e20 POT per year */
        @norm = 4.93e-4		       /* sets a normalization factor */
>

nuflux(#MICROminus)<			/* NAME can be replaced, e.g. FLUXplus */
	@flux_file ="./flux/MICROminus.dat"	/* loads the flux file */
	@time =	6.		/* sets the runtime in units years */
  	@power = 2.2            /* e20 POT per year */
        @norm =	4.93e-4		/* sets a normalization factor */
>

/* I can put arbitrary units and renormalize after the running. */

