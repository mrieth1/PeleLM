
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>


module prob_2D_module

  use fuego_chemistry

  implicit none

  private
  
  public :: amrex_probinit, setupbc, init_data

contains

! ::: -----------------------------------------------------------
! ::: This routine is called at problem initialization time
! ::: and when restarting from a checkpoint file.
! ::: The purpose is (1) to specify the initial time value
! ::: (not all problems start at time=0.0) and (2) to read
! ::: problem specific data from a namelist or other input
! ::: files and possibly store them or derived information
! ::: in FORTRAN common blocks for later use.
! ::: 
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: init      => TRUE if called at start of problem run
! :::              FALSE if called from restart
! ::: strttime <=  start problem with this time variable
! ::: 
! ::: -----------------------------------------------------------

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
  
      use PeleLM_F,  only: pphys_getP1atm_MKS
      use mod_Fvar_def, only : pamb, dpdt_factor, closed_chamber
      use mod_Fvar_def, only : fuelID, domnhi, domnlo, dim
      use mod_Fvar_def, only : ac_hist_file, cfix, changemax_control, &
                               coft_old, controlvelmax, corr, dv_control, &
                               h_control, navg_pnts, scale_control, sest, &
                               tau_control, tbase_control, V_in, v_in_old, zbase_control, &
                               pseudo_gravity
      use probdata_module, only : standoff, rho_bc, Y_bc
      use probdata_module, only : flame_dir
      
      
      implicit none
      
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i,istemp
      REAL_T area

      namelist /fortin/ V_in, flame_dir, &
                        standoff
      namelist /heattransin/ pamb, dpdt_factor

      namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, &
          zbase_control, pseudo_gravity, istemp,corr,controlVelMax,navg_pnts

!
!      Build `probin' filename -- the name of file containing fortin namelist.
!
      integer maxlen, isioproc
      parameter (maxlen=256)
      character probin*(maxlen)

      call bl_pd_is_ioproc(isioproc)

      if (init.ne.1) then
!         call bl_abort('probinit called with init ne 1')
      end if

      if (namlen .gt. maxlen) then
         call bl_abort('probin file name too long')
      end if

      if (namlen .eq. 0) then
         namlen = 6
         probin(1:namlen) = 'probin'
      else
         do i = 1, namlen
            probin(i:i) = char(name(i))
         end do
      endif

      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      
!     Set defaults
      pamb = pphys_getP1atm_MKS()
      dpdt_factor = 0.3d0
      closed_chamber = 0

      zbase_control = 0.d0

!     Note: for setup with no coflow, set Ro=Rf+wallth
      standoff = zero

!     Initialize control variables
      tau_control = one
      sest = zero
      corr = one
      changeMax_control = .05
      coft_old = -one
      cfix = zero
      ac_hist_file = 'AC_History'
      h_control = -one
      dV_control = zero
      tbase_control = zero
      h_control = -one
      pseudo_gravity = 0
      istemp = 0
      navg_pnts = 10

      read(untin,fortin)
      
!     Initialize control variables that depend on fortin variables
      V_in_old = V_in
      
      read(untin,heattransin)
 
      read(untin,control)
      close(unit=untin)

!     Set up boundary functions
      call setupbc()
      
      area = 1.d0
      do i=1,dim
        if (flame_dir /= i) then
         area = area*(domnhi(i)-domnlo(i))
        endif
      enddo
      scale_control = Y_bc(fuelID-1) * rho_bc(1) * area

      if (h_control .gt. zero) then
         cfix = scale_control * h_control
      endif


      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
         write(6,control)
      end if

  end subroutine amrex_probinit
  
!------------------------------------
  
  subroutine setupbc()bind(C, name="setupbc")

    use network,   only: nspec
    use PeleLM_F, only: pphys_getP1atm_MKS
    use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, maxspec, maxspnml, V_in
    use probdata_module, only : standoff, Y_bc, T_bc, u_bc, v_bc, rho_bc, h_bc
    use probdata_module, only : bcinit, flame_dir
  
    implicit none

    REAL_T Patm, pmf_vals(maxspec+3)
    REAL_T Xt(maxspec), Yt(maxspec), loc
    
    integer n
    integer b(2)
    data  b / 1, 1 /
      
    Patm = pamb / pphys_getP1atm_MKS()
             
  !     Take fuel mixture from pmf file
        loc = (domnlo(2)-standoff)*100.d0
        call pmf(loc,loc,pmf_vals,n)
        if (n.ne.Nspec+3) then
          call bl_pd_abort('setupbc: n(pmf) .ne. Nspec+3')
        endif
              
        do n = 1,Nspec
          Xt(n) = pmf_vals(3+n)
        end do 
              
        CALL CKXTY (Xt, Yt)
  
        do n=1,Nspec
          Y_bc(n-1) = Yt(n)
        end do
        
        T_bc = pmf_vals(1)
        if (flame_dir == 1) then
          v_bc = zero
          if (V_in .lt. 0) then
            u_bc = pmf_vals(2)*1.d-2
          else
            u_bc = V_in
          endif
        else if (flame_dir == 2) then
          u_bc = zero
          if (V_in .lt. 0) then
            v_bc = pmf_vals(2)*1.d-2
          else
            v_bc = V_in
          endif
        else
          write(6,*) 'Flame in Z dir not yet implemented'
          call bl_pd_abort(' ')
        end if 

!     Set density and hmix consistent with data

      call pphys_RHOfromPTY(b, b, &
                           rho_bc(1), DIMARG(b), DIMARG(b), &
                           T_bc(1),   DIMARG(b), DIMARG(b), &
                           Y_bc(0), DIMARG(b), DIMARG(b), Patm)
      call pphys_HMIXfromTY(b, b, &
                           h_bc(1),   DIMARG(b), DIMARG(b), &
                           T_bc(1),   DIMARG(b), DIMARG(b), &
                           Y_bc(0), DIMARG(b), DIMARG(b))

    bcinit = .true.

  end subroutine setupbc

      


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  The velocity field you
! ::: provide does not have to be divergence free and the pressure
! ::: field need not be set.  A subsequent projection iteration
! ::: will define aa divergence free velocity field along with a
! ::: consistant pressure.
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nscal     => number of scalar quantities.  You should know
! :::		   this already!
! ::: vel      <=  Velocity array
! ::: scal     <=  Scalar array
! ::: press    <=  Pressure array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------

  subroutine init_data(level,time,lo,hi,nscal, &
                       vel,scal,DIMS(state),press,DIMS(press), &
                       delta,xlo,xhi) &
                       bind(C, name="init_data")
                              

      use network,   only: nspec
      use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
      use mod_Fvar_def, only : bathID, domnhi, domnlo, maxspec, maxspnml 
      use probdata_module, only : standoff, flame_dir

      implicit none
      integer    level, nscal
      integer    lo(dim), hi(dim)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     xlo(dim), xhi(dim)
      REAL_T     time, delta(dim)
      REAL_T     vel(DIMV(state),dim)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))
      integer nPMF

      integer i, j, n
      REAL_T x, y, Yl(maxspec), Xl(maxspec), Patm
      REAL_T pmf_vals(maxspec+3)
      REAL_T pert,Lx, pos1, pos2

!      write(6,*)" made it to initdata"
      if (bathID.lt.1 .or. bathID.gt.Nspec) then
         call bl_pd_abort()
      endif

      do j = lo(2), hi(2)
        y = (float(j)+.5d0)*delta(2)+domnlo(2)
        do i = lo(1), hi(1)
          x = (float(i)+.5d0)*delta(1)+domnlo(1)
               
          pert = 0.d0
          
          if (flame_dir == 1) then
             pos1 = (x - standoff - 0.5d0*delta(1) + pert)*100.d0
             pos2 = (x - standoff + 0.5d0*delta(1) + pert)*100.d0
          elseif (flame_dir == 2) then
             pos1 = (y - standoff - 0.5d0*delta(2) + pert)*100.d0
             pos2 = (y - standoff + 0.5d0*delta(2) + pert)*100.d0
          else
             write(6,*) 'Flame in Z dir not yet implemented'
             call bl_pd_abort(' ')
          end if
       
          call pmf(pos1,pos2,pmf_vals,nPMF) 

          if (nPMF.ne.Nspec+3) then
            call bl_abort('INITDATA: n .ne. Nspec+3')
          endif
               
          scal(i,j,Temp) = pmf_vals(1)
          do n = 1,Nspec
            Xl(n) = pmf_vals(3+n)
          end do 
               
          CALL CKXTY (Xl, Yl)
               
          do n = 1,Nspec
            scal(i,j,FirstSpec+n-1) = Yl(n)
          end do

          scal(i,j,Trac) = 0.d0

          if (flame_dir == 1) then
             vel(i,j,1) = pmf_vals(2)*1.d-2
             vel(i,j,2) = 0.d0
          elseif (flame_dir == 2) then
             vel(i,j,1) = 0.d0
             vel(i,j,2) = pmf_vals(2)*1.d-2
          else
             write(6,*) 'Flame in Z dir not yet implemented'
             call bl_pd_abort(' ')
           end if

        end do
      end do

      Patm = pamb / pphys_getP1atm_MKS()
!      write(6,*)"Patm",Patm

      call pphys_RHOfromPTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state), &
          Patm)

      call pphys_HMIXfromTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state)) 

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
              scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
            enddo
            scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
         enddo
      enddo
      
  end subroutine init_data

end module prob_2D_module
