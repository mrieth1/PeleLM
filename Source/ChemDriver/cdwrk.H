      !
      ! Hardwired sizes.
      !
      integer, parameter :: maxspec  = 200
      integer, parameter :: maxspnml =  16


      ! SDC additions
      integer negative_Y_test
      double precision c_0(maxspec+1), c_1(maxspec+1), rhoH_INIT, T_cell, rhoY_INIT(maxspec)
      common /sdc_function/ c_0, c_1, rhoH_INIT, T_cell, rhoY_INIT, negative_Y_test
      save /sdc_function/
!$omp threadprivate(/sdc_function/)
