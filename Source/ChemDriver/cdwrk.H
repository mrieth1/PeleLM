      !
      ! The only ChemKin stuff we need.
      !
      integer, parameter :: ckbi = 1
      integer, parameter :: ckbr = 1
      !
      ! EGLib package for transport properties.
      !
      ! For normal LMC usage:
      !
      integer, parameter :: eg_nodes=1, eg_IFLAG=3, eg_ITLS=1
      !
      ! If you're using the DD* routines:
      !
      ! integer, parameter :: eg_nodes=1, eg_IFLAG=7, eg_ITLS=3
      !
      ! DVODE & DVODE driver workspace requirements.
      !
      integer, parameter :: dvbr = 1
      integer, parameter :: dvbi = 1

      integer dvr, dvi, dvdr, dvdbr, dvder
      common / wrkidx / dvr, dvi, dvdr, dvdbr, dvder
      save   / wrkidx /

      integer,          pointer :: IWRK(:)
      double precision, pointer :: RWRK(:)
      common / cdwork / RWRK, IWRK
      save   / cdwork /

!$omp threadprivate(/CDWORK/)
      !
      ! EGlib workspace -- allocated once in FORT_INITCHEM()
      !
      integer,          pointer :: EGIWRK(:)
      double precision, pointer :: EGRWRK(:)
      common / EGWORK / EGRWRK, EGIWRK
      save   / EGWORK /
      !
      ! Tranlib workspace -- allocated once in FORT_INITCHEM()
      !
      integer,          pointer :: MCIWRK(:)
      double precision, pointer :: MCRWRK(:)
      common / MCWORK / MCRWRK, MCIWRK
      save   / MCWORK /

!$omp threadprivate(/EGWORK/,/MCWORK/)
      !
      ! I/O
      !
      integer          LLINKMC
      common / ckdio / LLINKMC
      save   / ckdio /

      integer verbose_vode, max_vode_subcycles, vode_itol
      common / ckdio1 / verbose_vode, max_vode_subcycles, vode_itol
      save   / ckdio1 /
      !            
      ! Transport library
      !
      double precision  TMIN_TRANS
      common / tran2 / TMIN_TRANS
      save   / tran2 /

      integer Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset,use_mc,use_eg
      common / tran3 / Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset,use_mc,use_eg
      save   / tran3 /


