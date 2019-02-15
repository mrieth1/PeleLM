#include <actual_Creactor_unit.h> 

/**********************************/

/* Global Variables */
  N_Vector y_u            = NULL;
  SUNLinearSolver LS_unit = NULL;
  SUNMatrix A_unit        = NULL;
  int NEQ_unit            = 0;
  int iJac_Creact_u         = 0;
  int iE_Creact_u           = 1;
  int iverbose_u            = 1;
  void *cvode_mem_unit      = NULL;
  double *rhoe_init_u       = NULL;
  double *rhoh_init_u      = NULL;
  double *rhoesrc_ext_u    = NULL;
  double *rhohsrc_ext_u    = NULL;
  double *rYsrc_u          = NULL;

/**********************************/
/* Definitions */
/* Initialization routine, called once at the begining of the problem */
int extern_cInit_unit(const int* cvode_meth,const int* cvode_itmeth, 
		const int* cvode_iJac, const int* cvode_iE){

	int flag;
	realtype reltol, time;
	N_Vector atol;
	realtype *ratol;
	double rwrk;
	int iwrk;
	int mm, ii, nfit;
	int neq_tot;

	ckindx_(&iwrk,&rwrk,&mm,&NEQ_unit,&ii,&nfit);

        if (iverbose_u > 0) {
	    printf("(UNIT) Nb of spec is %d \n", NEQ_unit);
	}

	iJac_Creact_u    = *cvode_iJac;
	iE_Creact_u      = *cvode_iE;
        neq_tot          = (NEQ_unit + 1);

	/* Definition of main vector */
	y_u = N_VNew_Serial(neq_tot);
	if (check_flag_unit((void *)y_u, "N_VNew_Serial", 0)) return(1);

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	if ((*cvode_meth == 2) && (*cvode_itmeth == 2))
	{
	    cvode_mem_unit = CVodeCreate(CV_BDF);
	    if (check_flag_unit((void *)cvode_mem_unit, "CVodeCreate", 0)) return(1);
	} else {
	    amrex::Abort("\n--> Weird inputs to CVodeCreate. Viable options are CV_BDF (=2), CV_NEWTON (=2)\n");
	}

        time = 0.0e+0;
	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function, the inital time, and 
	 * initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem_unit, cF_RHS_unit, time, y_u);
	if (check_flag_unit(&flag, "CVodeInit", 1)) return(1);
	
	/* Definition of tolerances: one for each species */
	reltol = 1.0e-10;
        atol  = N_VNew_Serial(neq_tot);
	ratol = N_VGetArrayPointer(atol);
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = 1.0e-10;
        }
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem_unit, reltol, atol);
	if (check_flag_unit(&flag, "CVodeSVtolerances", 1)) return(1);

	flag = CVodeSetInitStep(cvode_mem_unit, 1.0e-09);
	if (check_flag_unit(&flag, "CVodeSetInitStep", 1)) return(1);

        printf("\n--> (UNIT) Using a Direct Dense Solver \n");

        /* Create dense SUNMatrix for use in linear solves */
	A_unit = SUNDenseMatrix(neq_tot, neq_tot);
	if(check_flag_unit((void *)A_unit, "SUNDenseMatrix", 0)) return(1);

	/* Create dense SUNLinearSolver object for use by CVode */
	LS_unit = SUNDenseLinearSolver(y_u, A_unit);
	if(check_flag_unit((void *)LS_unit, "SUNDenseLinearSolver", 0)) return(1);

	/* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	flag = CVDlsSetLinearSolver(cvode_mem_unit, LS_unit, A_unit);
	if(check_flag_unit(&flag, "CVDlsSetLinearSolver", 1)) return(1);


	if (iJac_Creact_u == 0) {
            printf("\n--> (UNIT) Without Analytical J\n");
	} else {
            printf("\n--> (UNIT) With Analytical J\n");
	    /* Set the user-supplied Jacobian routine Jac */
            flag = CVodeSetJacFn(cvode_mem_unit, cJac_unit);
	    if(check_flag_unit(&flag, "CVodeSetJacFn", 1)) return(1);
	}

        /* Set the max order */ 
        flag = CVodeSetMaxOrd(cvode_mem_unit, 2);
	if(check_flag_unit(&flag, "CVodeSetMaxOrd", 1)) return(1);

	/* Define vectors to be used later in creact */
	if (iE_Creact_u == 1) { 
	    rhoe_init_u   = (double *) malloc(sizeof(double));
	    rhoesrc_ext_u = (double *) malloc(sizeof(double));
	} else {
	    rhoh_init_u   = (double *) malloc(sizeof(double));
	    rhohsrc_ext_u = (double *) malloc(sizeof(double));
	}
	rYsrc_u = (double *)  malloc(NEQ_unit*sizeof(double));

	N_VDestroy(atol);          /* Free the atol vector */

	/* Ok we're done ...*/
        if (iverbose_u > 0) {
	    printf(" --> (UNIT) DONE WITH INITIALIZATION (CPU) %d \n", iE_Creact_u);
	}

	return(0);
}


/* Main CVODE call routine */
int actual_cReact_unit(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
                realtype *dt_react, realtype *time){

	realtype time_init, time_out, dummy_time ;
	int flag;

        time_init = *time;
	time_out  = *time + (*dt_react);

        if (iverbose_u > 3) {
	    printf("BEG : time curr is %14.6e and dt_react is %14.6e and final time should be %14.6e \n", time_init, *dt_react, time_out);
	}

	/* Get Device MemCpy of in arrays */
	/* Get Device pointer of solution vector */
	realtype *yvec_d      = N_VGetArrayPointer(y_u);
	// rhoY,T
	std::memcpy(yvec_d, rY_in, sizeof(realtype) * (NEQ_unit+1));
	// rhoY_src_ext
	std::memcpy(rYsrc_u, rY_src_in, NEQ_unit*sizeof(double));
	// rhoE/rhoH
	if (iE_Creact_u == 1) { 
	    std::memcpy(rhoe_init_u, rX_in, sizeof(realtype));
	    std::memcpy(rhoesrc_ext_u, rX_src_in, sizeof(realtype));
	} else {
	    std::memcpy(rhoh_init_u, rX_in, sizeof(realtype));
	    std::memcpy(rhohsrc_ext_u, rX_src_in, sizeof(realtype));
	}

	/* Call CVODE: ReInit for convergence */
        if (iverbose_u > 1) {
            printf("\n -----------------(UNIT)---------------\n");
	}
	CVodeReInit(cvode_mem_unit, time_init, y_u);

	flag = CVode(cvode_mem_unit, time_out, y_u, &dummy_time, CV_NORMAL);
	if (check_flag_unit(&flag, "CVode", 1)) return(1);

	*dt_react = dummy_time - time_init;
        if (iverbose_u > 3) {
	    printf("END : time curr is %14.6e and actual dt_react is %14.6e \n", dummy_time, *dt_react);
	}

	/* Pack data to return in main routine external */
	std::memcpy(rY_in, yvec_d, (NEQ_unit+1)*sizeof(realtype));
        *rX_in = *rX_in + (*dt_react) * (*rX_src_in);

	if (iverbose_u > 2) {
	    printf("\nAdditional verbose info --\n");
	    PrintFinalStats_unit(cvode_mem_unit);
            printf(" -------------------------------------\n");
	}

        long int nfe;
	flag = CVodeGetNumRhsEvals(cvode_mem_unit, &nfe);
	return nfe;
}


/* RHS routine used in CVODE */
static int cF_RHS_unit(realtype t, N_Vector y_in, N_Vector ydot_in, 
		void *user_data){

	realtype *y_d      = N_VGetArrayPointer(y_in);
	realtype *ydot_d   = N_VGetArrayPointer(ydot_in);

        if (iE_Creact_u == 1) {
	    fKernelSpec_unit(&t, y_d, ydot_d, 
			    rhoe_init_u, rhoesrc_ext_u, rYsrc_u);
	} else {
	    fKernelSpec_unit(&t, y_d, ydot_d, 
			    rhoh_init_u, rhohsrc_ext_u, rYsrc_u);
	}

	return(0);
}

/*
 * CUDA kernels
 */
void fKernelSpec_unit(realtype *dt, realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs)
{
      realtype massfrac[NEQ_unit],activity[NEQ_unit];
      realtype Xi[NEQ_unit], cXi[NEQ_unit];
      realtype cdot[NEQ_unit], molecular_weight[NEQ_unit];
      realtype temp, energy;
      int lierr;
      int * iwrk;
      double * rwrk;

      /* MW CGS */
      ckwt_(iwrk, rwrk, molecular_weight);
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NEQ_unit; i++){
          rho = rho + yvec_d[i];
      }
      /* temp */
      temp = yvec_d[NEQ_unit];
      /* Yks, C CGS*/
      for (int i = 0; i < NEQ_unit; i++){
          massfrac[i] = yvec_d[i] / rho;
	  activity[i] = yvec_d[i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      energy = (*rhoX_init + (*rhoXsrc_ext)*(*dt)) /rho;

      /* Fuego calls on device */
      if (iE_Creact_u == 1){
          get_t_given_ey_(&energy, massfrac, iwrk, rwrk, &temp, &lierr);
          ckums_(&temp, iwrk, rwrk, Xi);
          ckcvms_(&temp, iwrk, rwrk, cXi);
      } else {
          get_t_given_hy_(&energy, massfrac, iwrk, rwrk, &temp, &lierr);
          ckhms_(&temp, iwrk, rwrk, Xi);
          ckcpms_(&temp, iwrk, rwrk, cXi);
      }
      ckwc_(&temp, activity, iwrk, rwrk, cdot);
      int cX = 0.0;
      for (int i = 0; i < NEQ_unit; i++){
          cX = cX + massfrac[i] * cXi[i];
      }

      /* Fill ydot vect */
      ydot_d[NEQ_unit] = *rhoXsrc_ext;
      for (int i = 0; i < NEQ_unit; i++){
          ydot_d[i] = cdot[i] * molecular_weight[i] + rYs[i];
          ydot_d[NEQ_unit] = ydot_d[NEQ_unit]  - ydot_d[i] * Xi[i];
      }
      ydot_d[NEQ_unit] = ydot_d[NEQ_unit] /(rho * cX);
}


/* Analytical Jacobian evaluation */
static int cJac_unit(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  realtype *ydata  = N_VGetArrayPointer(u);

      realtype *J_col_k;
      realtype  temp; 
      realtype activity[NEQ_unit], molecular_weight[NEQ_unit];
      realtype Jmat_tmp[(NEQ_unit+1)*(NEQ_unit+1)];
      double rwrk;
      int iwrk;

      /* MW CGS */
      ckwt_(&iwrk, &rwrk, molecular_weight);
      /* temp */
      temp = ydata[NEQ_unit];
      for (int i = 0; i < NEQ_unit; i++){
          activity[i] = ydata[i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      int consP;
      if (iE_Creact_u == 1) {
	  consP = 0;
          dwdot_(Jmat_tmp, activity, &temp, &consP);
      } else {
          consP = 1;
          dwdot_(Jmat_tmp, activity, &temp, &consP);
      }
      /* fill the sunMat */
      for (int k = 0; k < NEQ_unit; k++){
	  J_col_k = SM_COLUMN_D(J,k);
	  for (int i = 0; i < NEQ_unit; i++){
	        J_col_k[i] = Jmat_tmp[k*(NEQ_unit+1)+i] * molecular_weight[i] / molecular_weight[k]; 
          }
	  J_col_k[NEQ_unit] = Jmat_tmp[k*(NEQ_unit+1)+NEQ_unit] / molecular_weight[k]; 
      }
      J_col_k = SM_COLUMN_D(J,NEQ_unit);
      for (int i = 0; i < NEQ_unit; i++){
          J_col_k[i] = Jmat_tmp[NEQ_unit*(NEQ_unit+1)+i] * molecular_weight[i]; 
      }

  return(0);

}

void extern_cFree_unit(){

	// Free y and abstol vectors 
	CVodeFree(&cvode_mem_unit);
	SUNLinSolFree(LS_unit);
        SUNMatDestroy(A_unit);
}


/* 
 * Get and print some final statistics
 */
static void PrintFinalStats_unit(void *cvodeMem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  long int nli, npe, nps, ncfl, netfails;
  int flag;
  realtype hlast, hinused, hcur;

  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag_unit(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag_unit(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag_unit(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag_unit(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag_unit(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netfails);
  check_flag_unit(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetLastStep(cvodeMem, &hlast);
  check_flag_unit(&flag, "CVodeGetLastStep", 1);
  flag = CVodeGetActualInitStep(cvodeMem, &hinused);
  check_flag_unit(&flag, "CVodeGetActualInitStep", 1);
  flag = CVodeGetCurrentTime(cvodeMem, &hcur);
  check_flag_unit(&flag, "CVodeGetCurrentTime", 1);

  flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
  check_flag_unit(&flag, "CVDlsGetNumRhsEvals", 1);
  flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
  check_flag_unit(&flag, "CVDlsGetNumJacEvals", 1);

  printf("-- Final Statistics --\n");
  printf("NonLinear (Newton) related --\n");
  printf("    DT(dt, dtcur), RHS, Iterations, ErrTestFails, LinSolvSetups = %-6ld(%14.6e %14.6e) %-6ld %-6ld %-6ld %-6ld \n",
  nst, hlast, hcur, nfe, nni, netfails, nsetups);

  printf("Linear (Dense Direct Solve) related --\n");
  printf("    FD RHS, NumJacEvals                           = %-6ld %-6ld \n", nfeLS, nje);
}


/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag_unit(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}






