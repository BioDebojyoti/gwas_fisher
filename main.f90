PROGRAM GWAS

   IMPLICIT NONE

   INTEGER, ALLOCATABLE :: allelecount(:), coverage(:)
   DOUBLE PRECISION, ALLOCATABLE :: frequency(:)
   CHARACTER(len=1), ALLOCATABLE :: sex(:)
   CHARACTER(len=11), ALLOCATABLE :: typeOfSNP(:)
   CHARACTER(len=10), ALLOCATABLE :: reference(:), allele(:)
   CHARACTER(len=16), ALLOCATABLE :: region(:)
   CHARACTER(len=258), ALLOCATABLE :: chromosome(:)

   LOGICAL :: I_OPEN, I_EXIST, II_OPEN, II_EXIST
   INTEGER :: I_NUMBER, II_NUMBER
   INTEGER :: i, j, currRow, matchindex
   INTEGER :: N
   CHARACTER(len=10) :: currRef
   CHARACTER(len=16) :: currRegion
   CHARACTER(len=258) :: currChr

!! Command line arguments
   INTEGER :: ii, iostat, error
   CHARACTER(len=100) :: filename, option1, option2
   CHARACTER(len=256) :: outputfile, conttablefile

!! Timing the execution
   REAL :: start, finish

!! Declaration for contingency table
    
   INTEGER :: iii, icr, jcr
   INTEGER :: ict, it, jt, ncol, nrow, sum
   DOUBLE PRECISION :: emin, expect, percnt, pre, prt, table(40), ctable(40)

   use fisher_exact_mod
  use io_utils
  use math_utils

!  TIME STARTS NOW
   CALL CPU_TIME(start)

!  File connections
   CALL GET_COMMAND_ARGUMENT(1,option1)
   CALL GET_COMMAND_ARGUMENT(2,option2)

   filename=option1
   READ(option2, '(I9)')N

   WRITE(*,*)
   WRITE(*,*)"File:",filename
   WRITE(*,*)"Number of lines to process:",N
   WRITE(*,*)

!  Dynamic allocation
   ALLOCATE(allelecount(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF
!  Dynamic allocation
   ALLOCATE(coverage(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF
!  Dynamic allocation
   ALLOCATE(frequency(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF
!  Dynamic allocation
   ALLOCATE(sex(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF
!  Dynamic allocation
   ALLOCATE(typeOfSNP(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF
!  Dynamic allocation
   ALLOCATE(reference(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF
!  Dynamic allocation
   ALLOCATE(allele(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF
!  Dynamic allocation
   ALLOCATE(region(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF
!  Dynamic allocation
   ALLOCATE(chromosome(N),stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't allocate memory, N = ",N
   STOP
   END IF

   OPEN(12, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=iostat)

   conttablefile="CONTINGENCY-TABLE-"//trim(filename)
   outputfile="OUTPUT-"//trim(filename)

   INQUIRE(FILE="subset.csv", OPENED=I_OPEN, EXIST=I_EXIST, NUMBER=I_NUMBER)

   IF (I_EXIST) THEN
    OPEN (UNIT=15, FILE=conttablefile, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
   ELSE
    OPEN (UNIT=15, FILE=conttablefile, STATUS='NEW', ACTION='WRITE')
   END IF

!  WRITING RESULTS OF FISHER EXACT TEST 

   INQUIRE(FILE="fisher.csv", OPENED=II_OPEN, EXIST=II_EXIST, NUMBER=II_NUMBER)

   IF (II_EXIST) THEN
    OPEN (UNIT=21, FILE=outputfile, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
   ELSE
    OPEN (UNIT=21, FILE=outputfile, STATUS='NEW', ACTION='WRITE')
   END IF

   WRITE(21,99)"Chromosome", "Region", "Reference", "p-value"
   99 format (3X,A15,36X,A15,3X,A15,3X,A15)

! Reading the input file
   DO i = 1, N
   READ(12,*)chromosome(i),region(i),typeOfSNP(i),reference(i),allele(i), &
   & allelecount(i),coverage(i),frequency(i),sex(i)
   END DO

   CLOSE(12)

   WRITE(*,1)"No. of Elements in chromosome: ",SIZE(chromosome)
   WRITE(*,1)"No. of Elements in region: ",SIZE(region)
   WRITE(*,1)"No. of Elements in typeOfSNP: ",SIZE(typeOfSNP)
   WRITE(*,1)"No. of Elements in reference: ",SIZE(reference)
   WRITE(*,1)"No. of Elements in allele: ",SIZE(allele)
   WRITE(*,1)"No. of Elements in allelecount: ",SIZE(allelecount)
   WRITE(*,1)"No. of Elements in coverage: ",SIZE(coverage)
   WRITE(*,1)"No. of Elements in frequency: ",SIZE(frequency)
   WRITE(*,1)"No. of Elements in sex: ",SIZE(sex)
   1 FORMAT(A30,I10)

   currRow = 1
   currChr = chromosome(1)
   currRegion = region(1)
   currRef = reference(1)

!  FISHER TEST PARAMETERS
   emin = 0.0
   expect = 0.0
   percnt = 0.0

   DO WHILE(currRow.LT.N) 

      matchindex = currRow

      DO WHILE((matchindex.NE.N).AND.(chromosome(matchindex) == currChr).AND.(region(matchindex) == currRegion))
         matchindex = matchindex + 1
      END DO

!  ict = index for contigency table
   ict = 1
!  To make sure all elements of the contingency table are not zero
   sum = 0

   DO j=currRow,(matchindex - 1)

      WRITE(15,15)chromosome(j),region(j),typeOfSNP(j),reference(j),allele(j), &
      & allelecount(j),coverage(j),frequency(j),sex(j)
      15 FORMAT (A60,3X,A15,3X,A10,3X,A5,3X,A5,3X,I5,3X,I5,3X,F10.5,3X,A2) 

!  Building the contigency table

      ctable(ict) = allelecount(j)
      sum = sum + ctable(ict) 
      ict = ict + 1

   END DO

!  Dimensions of the contingency table
   nrow = 2
   ncol = (ict -1)/ 2

!  REORDERING CONTINGENCY TABLE

      iii = 1
      DO icr=1,ncol
          DO jcr=1,nrow
             table(iii) = ctable(icr+(jcr-1)*ncol)
             iii = iii + 1
          END DO
      END DO
      
   WRITE(15,192)"CONTINGENCY TABLE WITH ",nrow," ROWS AND ", ncol," COLS"
   192 FORMAT (A23,I2,A9,I2,A5)
   DO  it=1, nrow
       write (15,197) (table(it+(jt-1)*nrow),jt=1,ncol)
   197 FORMAT (1x, 10f7.0)
   END DO

   IF(((NROW *NCOL).GE.4).AND.(sum>0)) THEN
   call fexact (nrow, ncol, table, nrow, expect, percnt, emin, prt, pre)
   WRITE(21,99999)currChr, currRegion, currRef, pre
 99999 format (A60,3X,A15,3X,A5,3X,F8.6)
   END IF

   currRow = matchindex
   currChr = chromosome(currRow)
   currRegion = region(currRow)
   currRef = reference(currRow)
   
   IF(currRow.GE.N) THEN 
      EXIT 
   END IF

   END DO

!  DO-WHILE LOOP ENDS HERE

   CLOSE(15)
   CLOSE(21)

!  Dynamic deallocation
   DEALLOCATE(allelecount,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF
!  Dynamic deallocation
   DEALLOCATE(coverage,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF
!  Dynamic deallocation
   DEALLOCATE(frequency,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF
!  Dynamic deallocation
   DEALLOCATE(sex,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF
!  Dynamic deallocation
   DEALLOCATE(typeOfSNP,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF
!  Dynamic deallocation
   DEALLOCATE(reference,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF
!  Dynamic deallocation
   DEALLOCATE(allele,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF
!  Dynamic deallocation
   DEALLOCATE(region,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF
!  Dynamic deallocation
   DEALLOCATE(chromosome,stat=error)
   IF(error.ne.0) THEN
   PRINT*,"error: couldn't deallocate memory, N = ",N
   STOP
   END IF

   CALL CPU_TIME(finish)

   WRITE(*,*)
   WRITE(*,*)"Time elapsed ",finish-start," seconds!!!"
   WRITE(*,*)

END PROGRAM GWAS

!!MODULE FISHER_EXACT
!!  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : int64

!!  CONTAINS
! ALGORITHM 643, COLLECTED ALGORITHMS FROM ACM.
! THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
! VOL. 19, NO. 4, DECEMBER, 1993, PP. 484-488.
! Driving program to read and 
! process the data file
!!      integer    i, j, ncol, nrow
!!      double precision emin, expect, percnt, pre, prt, table(200)
! 
!!      external   fexact
! Read table dimensions
!!   10 read (*,*) nrow, ncol, expect, percnt, emin
! Terminate on 0, 0
!!      if (nrow.eq.0 .and. ncol.eq.0) go to 9000
! Read and output TABLE
!!      write (*,99998)
!!      do i=1, nrow
!!         read (*,*) (table(i+(j-1)*nrow),j=1,ncol)
!!      end do
!
!!      do 20  i=1, nrow
!!         write (*,99997) (table(i+(j-1)*nrow),j=1,ncol)
!!99997    format (1x, 10f7.0)
!!   20 continue
!
!!      call fexact (nrow, ncol, table, nrow, expect, percnt, emin, prt, pre)
!
!!      write (*,99999) prt, pre
!!      go to 10
!
!!99998 format (/, 2x, 'The contingency table for this problem is:')
!!99999 format (2x, 'PRT =', f8.6, ' PRE = ', f8.6)
!! 9000 stop
!!      end
!-----------------------------------------------------------------------
! Name: FEXACT
!
! Purpose: Computes Fisher's exact test probabilities and a hybrid
! approximation to Fisher exact test probabilities for a
! contingency table using the network algorithm.
!
! Usage: CALL FEXACT (NROW, NCOL, TABLE, LDTABL, EXPECT, PERCNT,
! EMIN, PRT, PRE)
!
! Arguments:
! NROW - The number of rows in the table. (Input)
! NCOL - The number of columns in the table. (Input)
! TABLE - NROW by NCOL matrix containing the contingency table.
! (Input)
! LDTABL - Leading dimension of TABLE exactly as specified in the
! dimension statement in the calling program. (Input)
! EXPECT - Expected value used in the hybrid algorithm for
! deciding when to use asymptotic theory probabilities.
! (Input)
! If EXPECT .LE. 0.0 then asymptotic theory probabilities
! are not used and Fisher exact test probabilities are
! computed. Otherwise, if PERCNT or more of the cells in
! the remaining table have estimated expected values of
! EXPECT or more, with no remaining cell having expected
! value less than EMIN, then asymptotic chi-squared
! probabilities are used. See the algorithm section of the
! manual document for details. Use EXPECT = 5.0 to obtain
! the 'Cochran' condition.
! PERCNT - Percentage of remaining cells that must have estimated
! expected values greater than EXPECT before asymptotic
! probabilities can be used. (Input)
! See argument EXPECT for details. Use PERCNT = 80.0 to
! obtain the 'Cochran' condition.
! EMIN - Minimum cell estimated expected value allowed for
! asymptotic chi-squared probabilities to be used. (Input)
! See argument EXPECT for details. Use EMIN = 1.0 to
! obtain the 'Cochran' condition.
! PRT - Probability of the observed table for fixed marginal
! totals. (Output)
! PRE - Table p-value. (Output)
! PRE is the probability of a more extreme table, where
! 'extreme' is in a probabilistic sense.
! If EXPECT .LT. 0 then the Fisher exact probability
! is returned. Otherwise, an approximation to the
! Fisher exact probability is computed based upon
! asymptotic chi-squared probabilities for ``large''
! table expected values. The user defines ``large''
! through the arguments EXPECT, PERCNT, and EMIN.
!
! Remarks:
! 1. For many problems one megabyte or more of workspace can be 
! required. If the environment supports it, the user should begin 
! by increasing the workspace used to 200,000 units. 
!
! 2. In FEXACT, LDSTP = 30*LDKEY. The proportion of table space used 
! by STP may be changed by changing the line MULT = 30 below to 
! another value.
!
! 3. FEXACT may be converted to single precision by setting IREAL = 3,
! and converting all DOUBLE PRECISION specifications (except the 
! specifications for RWRK, IWRK, and DWRK) to REAL. This will 
! require changing the names and specifications of the intrinsic
! functions ALOG, AMAX1, AMIN1, EXP, and REAL. In addition, the
! machine specific constants will need to be changed, and the name
! DWRK will need to be changed to RWRK in the call to F2XACT.
!
! 4. Machine specific constants are specified and documented in F2XACT.
! A missing value code is specified in both FEXACT and F2XACT.
!
! 5. Although not a restriction, is is not generally practical to call
! this routine with large tables which are not sparse and in
! which the 'hybrid' algorithm has little effect. For example,
! although it is feasible to compute exact probabilities for the
! table
! 1 8 5 4 4 2 2
! 5 3 3 4 3 1 0
! 10 1 4 0 0 0 0,
! computing exact probabilities for a similar table which has been
! enlarged by the addition of an extra row (or column) may not be
! feasible.
!-----------------------------------------------------------------------
      subroutine fexact (nrow, ncol, table, ldtabl, expect, percnt,  &
     &                   emin, prt, pre)
! SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, ldtabl
      double precision expect, percnt, emin, prt, pre, table(ldtabl,*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, i1, i10, i2, i3, i3a, i3b, i3c, i4, i5, i6, i7,  &
     &           i8, i9, i9a, iiwk, ireal, irwk, iwkmax, iwkpt,  &
     &           j, k, kk, ldkey, ldstp, mult, nco, nro, & 
     &           ntot, numb
! SPECIFICATIONS FOR INTRINSICS
      intrinsic  max0
      integer    max0
! SPECIFICATIONS FOR SUBROUTINES
      external   prterr, f2xact
! SPECIFICATIONS FOR FUNCTIONS
      external   iwork
      integer    iwork
!***********************************************************************
! To increase workspace, increase the
! size of of rwrk and set the value of 
! IWKMAX to the new dimension
!
! When changing precision, the 
! following declarations should not be
! changed.
!***********************************************************************
      real       rwrk(8000000)
      double precision dwrk(8000000)
      integer    iwrk(8000000)
      equivalence (rwrk(1), iwrk(1)), (rwrk(1),dwrk(1))
! Set workspace size
      iwkmax = 8000000
!***********************************************************************
! To increase the length of the table
! of paste path lengths relative to the
! length of the hash table, increase
! MULT 
!***********************************************************************
      mult   = 40
!***********************************************************************
! Set IREAL = 4 for DOUBLE PRECISION
! Set IREAL = 3 for SINGLE PRECISION
!***********************************************************************
      ireal  = 4
!***********************************************************************
! AMISS is a missing value indicator
! which is returned when the
! probability is not defined.
!***********************************************************************
      amiss = -12345.0d0
!
      iwkpt  = 1
!
      if (nrow .gt. ldtabl) then
         call prterr (1, 'NROW must be less than or equal to '// &
     &               'LDTABL.')
      end if
      ntot = 0
      do 20  i=1, nrow
         do 10  j=1, ncol
            if (table(i,j) .lt. 0) then
               call prterr (2, 'All elements of TABLE must '//  &
     &                     'be positive.')
            end if
            ntot = ntot + table(i,j)
   10    continue
   20 continue
      if (ntot .eq. 0) then
         call prterr (3, 'All elements of TABLE are zero.  '//  &
     &               'PRT and PRE are set to missing values '// &
     &               '(NaN, not a number).')
         prt = amiss
         pre = amiss
         go to 9000
      end if
!
      nco = max0(nrow,ncol)
      nro = nrow + ncol - nco
      k   = nrow + ncol + 1
      kk  = k*max0(nrow,ncol)
!
      i1   = iwork(iwkmax,iwkpt,ntot+1,ireal)
      i2   = iwork(iwkmax,iwkpt,nco,2)
      i3   = iwork(iwkmax,iwkpt,nco,2)
      i3a  = iwork(iwkmax,iwkpt,nco,2)
      i3b  = iwork(iwkmax,iwkpt,nro,2)
      i3c  = iwork(iwkmax,iwkpt,nro,2)
      iiwk = iwork(iwkmax,iwkpt,max0(5*k+2*kk,800+7*max0(nrow,ncol)),2)
      irwk = iwork(iwkmax,iwkpt,max0(400+max0(nrow,ncol)+1,k),ireal)
! Double precision
      if (ireal .eq. 4) then
         numb  = 18 + 10*mult
         ldkey = (iwkmax-iwkpt+1)/numb
      else
! Real workspace
         numb  = 12 + 8*mult
         ldkey = (iwkmax-iwkpt+1)/numb
      end if
!
      ldstp = mult*ldkey
      i4    = iwork(iwkmax,iwkpt,2*ldkey,2)
      i5    = iwork(iwkmax,iwkpt,2*ldkey,2)
      i6    = iwork(iwkmax,iwkpt,2*ldstp,ireal)
      i7    = iwork(iwkmax,iwkpt,6*ldstp,2)
      i8    = iwork(iwkmax,iwkpt,2*ldkey,ireal)
      i9    = iwork(iwkmax,iwkpt,2*ldkey,ireal)
      i9a   = iwork(iwkmax,iwkpt,2*ldkey,ireal)
      i10   = iwork(iwkmax,iwkpt,2*ldkey,2)
!***********************************************************************
! To convert to double precision,
! change RWRK to WWRK in the next CALL
!***********************************************************************
!
      call f2xact (nrow, ncol, table, ldtabl, expect, percnt, emin,     &
     &             prt, pre, dwrk(i1), iwrk(i2), iwrk(i3), iwrk(i3a),   &
     &             iwrk(i3b), iwrk(i3c), iwrk(i4), ldkey, iwrk(i5),     &
     &             dwrk(i6), ldstp, iwrk(i7), dwrk(i8), dwrk(i9),       &
     &             dwrk(i9a), iwrk(i10), iwrk(iiwk), dwrk(irwk))
!
 9000 return
      end
!-----------------------------------------------------------------------
! Name: F2XACT
!
! Purpose: Computes Fisher's exact test for a contingency table,
! routine with workspace variables specified.
!
! Usage: CALL F2XACT (NROW, NCOL, TABLE, LDTABL, EXPECT, PERCNT,
! EMIN, PRT, PRE, FACT, ICO, IRO, KYY, IDIF,
! IRN, KEY, LDKEY, IPOIN, STP, LDSTP, IFRQ,
! DLP, DSP, TM, KEY2, IWK, RWK)
!-----------------------------------------------------------------------
      subroutine f2xact (nrow, ncol, table, ldtabl, expect, percnt,    &
     &                   emin, prt, pre, fact, ico, iro, kyy, idif,    &
     &                   irn, key, ldkey, ipoin, stp, ldstp, ifrq,     &
     &                   dlp, dsp, tm, key2, iwk, rwk)
! SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, ldtabl, ldkey, ldstp, ico(*), iro(*),     &
     &           kyy(*), idif(*), irn(*), key(*), ipoin(*), ifrq(*),   &
     &           key2(*), iwk(*)
      double precision expect, percnt, emin, prt, pre, table(ldtabl,*), &
     &           fact(0:*), stp(*), dlp(*), dsp(*), tm(*), rwk(*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, i31, i310, i311, i32, i33, i34, i35, i36, i37,     &
     &           i38, i39, i41, i42, i43, i44, i45, i46, i47, i48,     &
     &           iflag, ifreq, ii, ikkey, ikstp, ikstp2, ipn, ipo,     &
     &           itmp, itop, itp, j, jkey, jstp, jstp2, jstp3, jstp4,  &
     &           k, k1, kb, kd, kmax, ks, kval, last, n, ncell, nco,   &
     &           nrb, nro, nro2, ntot, ifault, imax
      double precision dd, ddf, df, drn, dro, dspt, emn, obs, obs2,    &
     &           obs3, pastp, pv, tmp, tol
      logical    chisq, ipsh
! SPECIFICATIONS FOR INTRINSICS
      intrinsic  dlog, dmax1, dmin1, dexp, max0, min0, mod, nint, dble
      integer    max0, min0, mod, nint
      double precision dlog, dmax1, dmin1, dexp, dble
! SPECIFICATIONS FOR SUBROUTINES
      external   prterr, f3xact, f4xact, f5xact, f6xact, f7xact, isort
! SPECIFICATIONS FOR FUNCTIONS
      external   f9xact, gammds
      double precision f9xact, gammds
!***********************************************************************
! IMAX is the largest representable
! integer on the machine
!***********************************************************************
      data imax/2147483647/
!***********************************************************************
! AMISS is a missing value indicator
! which is returned when the
! probability is not defined.
!***********************************************************************
      data amiss/-12345.0d0/
!***********************************************************************
! TOL is chosen as the square root of
! the smallest relative spacing
!***********************************************************************
      data tol/3.45254d-07/
!***********************************************************************
! EMX is a large positive value used 
! in comparing expected values
!***********************************************************************
      data emx/1.0d30/
! Initialize KEY array
      do 10  i=1, 2*ldkey
         key(i)  = -9999
         key2(i) = -9999
   10 continue
! Initialize parameters
      pre  = 0.0
      itop = 0
      if (expect .gt. 0.0d0) then
         emn = emin
      else
         emn = emx
      end if
! Initialize pointers for workspace
      k = max0(nrow,ncol)
! f3xact
      i31  = 1
      i32  = i31 + k
      i33  = i32 + k
      i34  = i33 + k
      i35  = i34 + k
      i36  = i35 + k
      i37  = i36 + k
      i38  = i37 + k
      i39  = i38 + 400
      i310 = 1
      i311 = 401
! f4xact
      k   = nrow + ncol + 1
      i41 = 1
      i42 = i41 + k
      i43 = i42 + k
      i44 = i43 + k
      i45 = i44 + k
      i46 = i45 + k
      i47 = i46 + k*max0(nrow,ncol)
      i48 = 1
! Check table dimensions
      if (nrow .gt. ldtabl) then
         call prterr (1, 'NROW must be less than or equal to '// &
     &               'LDTABL.')
      end if
      if (ncol .le. 1) then
         call prterr (4, 'NCOL must be greater than 1.0.')
      end if
! Compute row marginals and total
      ntot = 0
      do 30  i=1, nrow
         iro(i) = 0
         do 20  j=1, ncol
            if (table(i,j) .lt. -0.0001d0) then
               call prterr (2, 'All elements of TABLE must be '//  &
     &                     'positive.')
            end if
            iro(i) = iro(i) + nint(table(i,j))
            ntot   = ntot + nint(table(i,j))
   20    continue
   30 continue
!
      if (ntot .eq. 0) then
         call prterr (3, 'All elements of TABLE are zero.  '//   &
     &               'PRT and PRE are set to missing values '//  &
     &               '(NaN, not a number).')
         prt = amiss
         pre = amiss
         go to 9000
      end if
! Column marginals
      do 50  i=1, ncol
         ico(i) = 0
         do 40  j=1, nrow
            ico(i) = ico(i) + nint(table(j,i))
   40    continue
   50 continue
! sort
      call isort (nrow, iro)
      call isort (ncol, ico)
! Determine row and column marginals
!
      if (nrow .gt. ncol) then
         nro = ncol
         nco = nrow
! Interchange row and column marginals
         do 60  i=1, nrow
            itmp = iro(i)
            if (i .le. ncol) iro(i) = ico(i)
            ico(i) = itmp
   60    continue
      else
         nro = nrow
         nco = ncol
      end if
!
! Get multiplers for stack
      kyy(1) = 1
      do 70  i=2, nro
! Hash table multipliers
         if (iro(i-1)+1 .le. imax/kyy(i-1)) then
            kyy(i) = kyy(i-1)*(iro(i-1)+1)
            j      = j/kyy(i-1)
         else
            call prterr (5, 'The hash table key cannot be computed'//    &
     &                  ' because the largest key is larger than the'//  &
     &                  ' largest representable integer.  The '//        & 
     &                  'algorithm cannot proceed.')
         end if
   70 continue
! Maximum product
      if (iro(nro-1)+1 .le. imax/kyy(nro-1)) then
         kmax = (iro(nro)+1)*kyy(nro-1)
      else
         call prterr (5, 'The hash table key cannot be computed'//      &
     &               ' because the largest key is larger than the'//    &
     &               ' largest representable integer.  The '//          &
     &               'algorithm cannot proceed.')
         go to 9000
      end if
! Compute log factorials
      fact(0) = 0.0d0
      fact(1) = 0.0d0
      fact(2) = dlog(2.0d0)
      do 80  i=3, ntot, 2
         fact(i) = fact(i-1) + dlog(dble(i))
         j       = i + 1
         if (j .le. ntot) fact(j) = fact(i) + fact(2) + fact(j/2) -    &
     &       fact(j/2-1)
   80 continue
! Compute observed path length: OBS
      obs  = tol
      ntot = 0
      do 100  j=1, nco
         dd = 0.0
         do 90  i=1, nro
            if (nrow .le. ncol) then
               dd   = dd + fact(nint(table(i,j)))
               ntot = ntot + nint(table(i,j))
            else
               dd   = dd + fact(nint(table(j,i)))
               ntot = ntot + nint(table(j,i))
            end if
   90    continue
         obs = obs + fact(ico(j)) - dd
  100 continue
! Denominator of observed table: DRO
      dro = f9xact(nro,ntot,iro,fact)
      prt = dexp(obs-dro)
! Initialize pointers
      k        = nco
      last     = ldkey + 1
      jkey     = ldkey + 1
      jstp     = ldstp + 1
      jstp2    = 3*ldstp + 1
      jstp3    = 4*ldstp + 1
      jstp4    = 5*ldstp + 1
      ikkey    = 0
      ikstp    = 0
      ikstp2   = 2*ldstp
      ipo      = 1
      ipoin(1) = 1
      stp(1)   = 0.0
      ifrq(1)  = 1
      ifrq(ikstp2+1) = -1
!
  110 kb = nco - k + 1
      ks   = 0
      n    = ico(kb)
      kd   = nro + 1
      kmax = nro
! IDIF is the difference in going to th
! daughter
      do 120  i=1, nro
         idif(i) = 0
  120 continue
! Generate the first daughter
  130 kd = kd - 1
      ntot     = min0(n,iro(kd))
      idif(kd) = ntot
      if (idif(kmax) .eq. 0) kmax = kmax - 1
      n = n - ntot
      if (n.gt.0 .and. kd.ne.1) go to 130
      if (n .ne. 0) go to 310
!
      k1   = k - 1
      n    = ico(kb)
      ntot = 0
      do 140  i=kb + 1, nco
         ntot = ntot + ico(i)
  140 continue
! Arc to daughter length=ICO(KB)
  150 do 160  i=1, nro
         irn(i) = iro(i) - idif(i)
  160 continue
! Sort irn
      if (k1 .gt. 1) then
         if (nro .eq. 2) then
            if (irn(1) .gt. irn(2)) then
               ii     = irn(1)
               irn(1) = irn(2)
               irn(2) = ii
            end if
         else if (nro .eq. 3) then
            ii = irn(1)
            if (ii .gt. irn(3)) then
               if (ii .gt. irn(2)) then
                  if (irn(2) .gt. irn(3)) then
                     irn(1) = irn(3)
                     irn(3) = ii
                  else
                     irn(1) = irn(2)
                     irn(2) = irn(3)
                     irn(3) = ii
                  end if
               else
                  irn(1) = irn(3)
                  irn(3) = irn(2)
                  irn(2) = ii
               end if
            else if (ii .gt. irn(2)) then
               irn(1) = irn(2)
               irn(2) = ii
            else if (irn(2) .gt. irn(3)) then
               ii     = irn(2)
               irn(2) = irn(3)
               irn(3) = ii
            end if
         else
            do 180  j=2, nro
               i  = j - 1
               ii = irn(j)
  170          if (ii .lt. irn(i)) then
                  irn(i+1) = irn(i)
                  i        = i - 1
                  if (i .gt. 0) go to 170
               end if
               irn(i+1) = ii
  180       continue
         end if
! Adjust start for zero
         do 190  i=1, nro
            if (irn(i) .ne. 0) go to 200
  190    continue
  200    nrb = i
         nro2 = nro - i + 1
      else
         nrb  = 1
         nro2 = nro
      end if
! Some table values
      ddf = f9xact(nro,n,idif,fact)
      drn = f9xact(nro2,ntot,irn(nrb),fact) - dro + ddf
! Get hash value
      if (k1 .gt. 1) then
         kval = irn(1) + irn(2)*kyy(2)
         do 210  i=3, nro
            kval = kval + irn(i)*kyy(i)
  210    continue
! Get hash table entry
         i = mod(kval,2*ldkey) + 1
! Search for unused location
         do 220  itp=i, 2*ldkey
            ii = key2(itp)
            if (ii .eq. kval) then
               go to 240
            else if (ii .lt. 0) then
               key2(itp) = kval
               dlp(itp)  = 1.0d0
               dsp(itp)  = 1.0d0
               go to 240
            end if
  220    continue
!
         do 230  itp=1, i - 1
            ii = key2(itp)
            if (ii .eq. kval) then
               go to 240
            else if (ii .lt. 0) then
               key2(itp) = kval
               dlp(itp)  = 1.0
               go to 240
            end if
  230    continue
!
         call prterr (6, 'LDKEY is too small.  It is not possible to '// &
     &               'give thevalue of LDKEY required, but you could '// &
     &               'try doubling LDKEY (and possibly LDSTP).')
      end if
!
  240 ipsh = .true.
! Recover pastp
      ipn   = ipoin(ipo+ikkey)
      pastp = stp(ipn+ikstp)
      ifreq = ifrq(ipn+ikstp)
! Compute shortest and longest path
      if (k1 .gt. 1) then
         obs2 = obs - fact(ico(kb+1)) - fact(ico(kb+2)) - ddf
         do 250  i=3, k1
            obs2 = obs2 - fact(ico(kb+i))
  250    continue
!
         if (dlp(itp) .gt. 0.0d0) then
            dspt = obs - obs2 - ddf
! Compute longest path
            dlp(itp) = 0.0d0
            call f3xact (nro2, irn(nrb), k1, ico(kb+1), dlp(itp),     &
     &                   ntot, fact, iwk(i31), iwk(i32), iwk(i33),    &
     &                   iwk(i34), iwk(i35), iwk(i36), iwk(i37),      &
     &                   iwk(i38), iwk(i39), rwk(i310), rwk(i311), tol)
            dlp(itp) = dmin1(0.0d0,dlp(itp))
! Compute shortest path
            dsp(itp) = dspt
            call f4xact (nro2, irn(nrb), k1, ico(kb+1), dsp(itp),       &
     &                   fact, iwk(i47), iwk(i41), iwk(i42), iwk(i43),  &
     &                   iwk(i44), iwk(i45), iwk(i46), rwk(i48), tol)
            dsp(itp) = dmin1(0.0d0,dsp(itp)-dspt)
! Use chi-squared approximation?
            if (dble(irn(nrb)*ico(kb+1))/dble(ntot) .gt. emn) then
               ncell = 0.0
               do 270  i=1, nro2
                  do 260  j=1, k1
                     if (irn(nrb+i-1)*ico(kb+j) .ge. ntot*expect) then
                        ncell = ncell + 1
                     end if
  260             continue
  270          continue
               if (ncell*100 .ge. k1*nro2*percnt) then
                  tmp = 0.0
                  do 280  i=1, nro2
                     tmp = tmp + fact(irn(nrb+i-1)) -      &
     &                     fact(irn(nrb+i-1)-1)
  280             continue
                  tmp = tmp*(k1-1)
                  do 290  j=1, k1
                     tmp = tmp + (nro2-1)*(fact(ico(kb+j))-fact(ico(kb+  &
     &                     j)-1))
  290             continue
                  df      = (nro2-1)*(k1-1)
                  tmp     = tmp + df*1.83787706640934548356065947281d0
                  tmp     = tmp - (nro2*k1-1)*(fact(ntot)-fact(ntot-1))
                  tm(itp) = -2.0d0*(obs-dro) - tmp
               else
! tm(itp) set to a flag value
                  tm(itp) = -9876.0d0
               end if
            else
               tm(itp) = -9876.0d0
            end if
         end if
         obs3 = obs2 - dlp(itp)
         obs2 = obs2 - dsp(itp)
         if (tm(itp) .eq. -9876.0d0) then
            chisq = .false.
         else
            chisq = .true.
            tmp   = tm(itp)
         end if
      else
         obs2 = obs - drn - dro
         obs3 = obs2
      end if
! Process node with new PASTP
  300 if (pastp .le. obs3) then
! Update pre
         pre = pre + dble(ifreq)*dexp(pastp+drn)
!
      else if (pastp .lt. obs2) then
         if (chisq) then
            df  = (nro2-1)*(k1-1)
            pv  = 1.0 - gammds(dmax1(0.0d0,tmp+2.0d0*(pastp+drn))/     &
     &            2.0d0,df/2.0d0,ifault)
            pre = pre + dble(ifreq)*dexp(pastp+drn)*pv
         else
! Put daughter on queue
            call f5xact (pastp+ddf, tol, kval, key(jkey), ldkey,        &
     &                   ipoin(jkey), stp(jstp), ldstp, ifrq(jstp),     &
     &                   ifrq(jstp2), ifrq(jstp3), ifrq(jstp4), ifreq,  &
     &                   itop, ipsh)
            ipsh = .false.
         end if
      end if
! Get next PASTP on chain
      ipn = ifrq(ipn+ikstp2)
      if (ipn .gt. 0) then
         pastp = stp(ipn+ikstp)
         ifreq = ifrq(ipn+ikstp)
         go to 300
      end if
! Generate a new daughter node
      call f7xact (kmax, iro, idif, kd, ks, iflag)
      if (iflag .ne. 1) go to 150
! Go get a new mother from stage K
  310 iflag = 1
      call f6xact (nro, iro, iflag, kyy, key(ikkey+1), ldkey, last,    &
     &             ipo)
! Update pointers
      if (iflag .eq. 3) then
         k      = k - 1
         itop   = 0
         ikkey  = jkey - 1
         ikstp  = jstp - 1
         ikstp2 = jstp2 - 1
         jkey   = ldkey - jkey + 2
         jstp   = ldstp - jstp + 2
         jstp2  = 2*ldstp + jstp
         do 320  i=1, 2*ldkey
            key2(i) = -9999
  320    continue
         if (k .ge. 2) go to 310
      else
         go to 110
      end if
!
 9000 return
      end
!-----------------------------------------------------------------------
! Name: F3XACT
!
! Purpose: Computes the shortest path length for a given table.
!
! Usage: CALL F3XACT (NROW, IROW, NCOL, ICOL, DLP, MM, FACT, ICO,
! IRO, IT, LB, NR, NT, NU, ITC, IST, STV, ALEN,
! TOL)
!
! Arguments:
! NROW - The number of rows in the table. (Input)
! IROW - Vector of length NROW containing the row sums for the
! table. (Input)
! NCOL - The number of columns in the table. (Input)
! ICOL - Vector of length K containing the column sums for the
! table. (Input)
! DLP - The longest path for the table. (Output)
! MM - The total count in the table. (Output)
! FACT - Vector containing the logarithms of factorials. (Input)
! ICO - Work vector of length MAX(NROW,NCOL).
! IRO - Work vector of length MAX(NROW,NCOL).
! IT - Work vector of length MAX(NROW,NCOL).
! LB - Work vector of length MAX(NROW,NCOL).
! NR - Work vector of length MAX(NROW,NCOL).
! NT - Work vector of length MAX(NROW,NCOL).
! NU - Work vector of length MAX(NROW,NCOL).
! ITC - Work vector of length 400.
! IST - Work vector of length 400.
! STV - Work vector of length 400.
! ALEN - Work vector of length MAX(NROW,NCOL).
! TOL - Tolerance. (Input)
!-----------------------------------------------------------------------
      subroutine f3xact (nrow, irow, ncol, icol, dlp, mm, fact, ico,      &
     &                   iro, it, lb, nr, nt, nu, itc, ist, stv, alen,    &
     &                   tol)
! SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, mm, irow(*), icol(*), ico(*), iro(*),        &
     &           it(*), lb(*), nr(*), nt(*), nu(*), itc(*), ist(*)
      double precision dlp, tol, fact(0:*), stv(*), alen(0:*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, ic1, ic2, ii, ipn, irl, itp, k, key, ks, kyy, lev,    &
     &           n11, n12, nc1, nc1s, nco, nct, nn, nn1, nr1, nro, nrt
      double precision v, val, vmn
      logical    xmin
! SPECIFICATIONS FOR SAVE VARIABLES
      integer    ldst, nitc, nst
      save       ldst, nitc, nst
! SPECIFICATIONS FOR INTRINSICS
      intrinsic  dmin1, int, mod, dble
      integer    int, mod
      double precision dmin1, dble
! SPECIFICATIONS FOR SUBROUTINES
      external   prterr, f10act, isort
!
      data ldst/200/, nst/0/, nitc/0/
!
      do 10  i=0, ncol
         alen(i) = 0.0
   10 continue
      do 20  i=1, 400
         ist(i) = -1
   20 continue
! nrow is 1
      if (nrow .le. 1) then
         if (nrow .gt. 0) then
            dlp = dlp - fact(icol(1))
            do 30  i=2, ncol
               dlp = dlp - fact(icol(i))
   30       continue
         end if
         go to 9000
      end if
! ncol is 1
      if (ncol .le. 1) then
         if (ncol .gt. 0) then
            dlp = dlp - fact(irow(1)) - fact(irow(2))
            do 40  i=3, nrow
               dlp = dlp - fact(irow(i))
   40       continue
         end if
         go to 9000
      end if
! 2 by 2 table
      if (nrow*ncol .eq. 4) then
         n11 = (irow(1)+1)*(icol(1)+1)/(mm+2)
         n12 = irow(1) - n11
         dlp = dlp - fact(n11) - fact(n12) - fact(icol(1)-n11) -    &
     &         fact(icol(2)-n12)
         go to 9000
      end if
! Test for optimal table
      val  = 0.0
      xmin = .false.
      if (irow(nrow) .le. irow(1)+ncol) then
         call f10act (nrow, irow, ncol, icol, val, xmin, fact, lb, nu, &
     &                nr)
      end if
      if (.not.xmin) then
         if (icol(ncol) .le. icol(1)+nrow) then
            call f10act (ncol, icol, nrow, irow, val, xmin, fact, lb,  &
     &                   nu, nr)
         end if
      end if
!
      if (xmin) then
         dlp = dlp - val
         go to 9000
      end if
! Setup for dynamic programming
      nn = mm
! Minimize ncol
      if (nrow .ge. ncol) then
         nro = nrow
         nco = ncol
!
         do 50  i=1, nrow
            iro(i) = irow(i)
   50    continue
!
         ico(1) = icol(1)
         nt(1)  = nn - ico(1)
         do 60  i=2, ncol
            ico(i) = icol(i)
            nt(i)  = nt(i-1) - ico(i)
   60    continue
      else
         nro = ncol
         nco = nrow
!
         ico(1) = irow(1)
         nt(1)  = nn - ico(1)
         do 70  i=2, nrow
            ico(i) = irow(i)
            nt(i)  = nt(i-1) - ico(i)
   70    continue
!
         do 80  i=1, ncol
            iro(i) = icol(i)
   80    continue
      end if
! Initialize pointers
      vmn  = 1.0d10
      nc1s = nco - 1
      irl  = 1
      ks   = 0
      k    = ldst
      kyy  = ico(nco) + 1
      go to 100
! Test for optimality
   90 xmin = .false.
      if (iro(nro) .le. iro(irl)+nco) then
         call f10act (nro, iro(irl), nco, ico, val, xmin, fact, lb, nu, nr)
      end if
      if (.not.xmin) then
         if (ico(nco) .le. ico(1)+nro) then
            call f10act (nco, ico, nro, iro(irl), val, xmin, fact, lb, nu, nr)
         end if
      end if
!
      if (xmin) then
         if (val .lt. vmn) vmn = val
         go to 200
      end if
! Setup to generate new node
  100 lev = 1
      nr1   = nro - 1
      nrt   = iro(irl)
      nct   = ico(1)
      lb(1) = int(dble((nrt+1)*(nct+1))/dble(nn+nr1*nc1s+1)-tol) - 1
      nu(1) = int(dble((nrt+nc1s)*(nct+nr1))/dble(nn+nr1+nc1s)) - lb(1) + 1
      nr(1) = nrt - lb(1)
! Generate a node
  110 nu(lev) = nu(lev) - 1
      if (nu(lev) .eq. 0) then
         if (lev .eq. 1) go to 200
         lev = lev - 1
         go to 110
      end if
      lb(lev) = lb(lev) + 1
      nr(lev) = nr(lev) - 1
  120 alen(lev) = alen(lev-1) + fact(lb(lev))
      if (lev .lt. nc1s) then
         nn1     = nt(lev)
         nrt     = nr(lev)
         lev     = lev + 1
         nc1     = nco - lev
         nct     = ico(lev)
         lb(lev) = dble((nrt+1)*(nct+1))/dble(nn1+nr1*nc1+1) - tol
         nu(lev) = dble((nrt+nc1)*(nct+nr1))/dble(nn1+nr1+nc1) - lb(lev) + 1
         nr(lev) = nrt - lb(lev)
         go to 120
      end if
      alen(nco) = alen(lev) + fact(nr(lev))
      lb(nco)   = nr(lev)
!
      v = val + alen(nco)
      if (nro .eq. 2) then
! Only 1 row left
         v = v + fact(ico(1)-lb(1)) + fact(ico(2)-lb(2))
         do 130  i=3, nco
            v = v + fact(ico(i)-lb(i))
  130    continue
         if (v .lt. vmn) vmn = v
      else if (nro.eq.3 .and. nco.eq.2) then
! 3 rows and 2 columns
         nn1 = nn - iro(irl) + 2
         ic1 = ico(1) - lb(1)
         ic2 = ico(2) - lb(2)
         n11 = (iro(irl+1)+1)*(ic1+1)/nn1
         n12 = iro(irl+1) - n11
         v   = v + fact(n11) + fact(n12) + fact(ic1-n11) + fact(ic2-n12)
         if (v .lt. vmn) vmn = v
      else
! Column marginals are new node
         do 140  i=1, nco
            it(i) = ico(i) - lb(i)
  140    continue
! Sort column marginals
         if (nco .eq. 2) then
            if (it(1) .gt. it(2)) then
               ii    = it(1)
               it(1) = it(2)
               it(2) = ii
            end if
         else if (nco .eq. 3) then
            ii = it(1)
            if (ii .gt. it(3)) then
               if (ii .gt. it(2)) then
                  if (it(2) .gt. it(3)) then
                     it(1) = it(3)
                     it(3) = ii
                  else
                     it(1) = it(2)
                     it(2) = it(3)
                     it(3) = ii
                  end if
               else
                  it(1) = it(3)
                  it(3) = it(2)
                  it(2) = ii
               end if
            else if (ii .gt. it(2)) then
               it(1) = it(2)
               it(2) = ii
            else if (it(2) .gt. it(3)) then
               ii    = it(2)
               it(2) = it(3)
               it(3) = ii
            end if
         else
            call isort (nco, it)
         end if
! Compute hash value
         key = it(1)*kyy + it(2)
         do 150  i=3, nco
            key = it(i) + key*kyy
  150    continue
! Table index
         ipn = mod(key,ldst) + 1
! Find empty position
         ii = ks + ipn
         do 160  itp=ipn, ldst
            if (ist(ii) .lt. 0) then
               go to 180
            else if (ist(ii) .eq. key) then
               go to 190
            end if
            ii = ii + 1
  160    continue
!
         ii = ks + 1
         do 170  itp=1, ipn - 1
            if (ist(ii) .lt. 0) then
               go to 180
            else if (ist(ii) .eq. key) then
               go to 190
            end if
            ii = ii + 1
  170    continue
!
         call prterr (30, 'Stack length exceeded in f3xact.'// '  This problem should not occur.')
! Push onto stack
  180    ist(ii) = key
         stv(ii) = v
         nst     = nst + 1
         ii      = nst + ks
         itc(ii) = itp
         go to 110
! Marginals already on stack
  190    stv(ii) = dmin1(v,stv(ii))
      end if
      go to 110
! Pop item from stack
  200 if (nitc .gt. 0) then
! Stack index
         itp      = itc(nitc+k) + k
         nitc     = nitc - 1
         val      = stv(itp)
         key      = ist(itp)
         ist(itp) = -1
! Compute marginals
         do 210  i=nco, 2, -1
            ico(i) = mod(key,kyy)
            key    = key/kyy
  210    continue
         ico(1) = key
! Set up nt array
         nt(1) = nn - ico(1)
         do 220  i=2, nco
            nt(i) = nt(i-1) - ico(i)
  220    continue
         go to 90
!
      else if (nro.gt.2 .and. nst.gt.0) then
! Go to next level
         nitc = nst
         nst  = 0
         k    = ks
         ks   = ldst - ks
         nn   = nn - iro(irl)
         irl  = irl + 1
         nro  = nro - 1
         go to 200
      end if
!
      dlp = dlp - vmn
 9000 return
      end
!-----------------------------------------------------------------------
! Name: F4XACT
!
! Purpose: Computes the longest path length for a given table.
!
! Usage: CALL F4XACT (NROW, IROW, NCOL, ICOL, DSP, FACT, ICSTK,
! NCSTK, LSTK, MSTK, NSTK, NRSTK, IRSTK, YSTK,
! TOL)
!
! Arguments:
! NROW - The number of rows in the table. (Input)
! IROW - Vector of length NROW containing the row sums for the
! table. (Input)
! NCOL - The number of columns in the table. (Input)
! ICOL - Vector of length K containing the column sums for the
! table. (Input)
! DSP - The shortest path for the table. (Output)
! FACT - Vector containing the logarithms of factorials. (Input)
! ICSTK - NCOL by NROW+NCOL+1 work array.
! NCSTK - Work vector of length NROW+NCOL+1.
! LSTK - Work vector of length NROW+NCOL+1.
! MSTK - Work vector of length NROW+NCOL+1.
! NSTK - Work vector of length NROW+NCOL+1.
! NRSTK - Work vector of length NROW+NCOL+1.
! IRSTK - NROW by MAX(NROW,NCOL) work array.
! YSTK - Work vector of length NROW+NCOL+1.
! TOL - Tolerance. (Input)
!-----------------------------------------------------------------------
      subroutine f4xact (nrow, irow, ncol, icol, dsp, fact, icstk,         &
     &                   ncstk, lstk, mstk, nstk, nrstk, irstk, ystk,      &
     &                   tol)
! SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, irow(*), icol(*), icstk(ncol,*),              &
     &           ncstk(*), lstk(*), mstk(*), nstk(*), nrstk(*),            &
     &           irstk(nrow,*)
      double precision dsp, tol, fact(0:*), ystk(*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, ic1, ict, ir1, irt, istk, j, k, l, m, mn, n, nco,      &
     &           nro
      double precision amx, y
! SPECIFICATIONS FOR SUBROUTINES
      external   f11act, f8xact
! Take care of the easy cases firstkt
      if (nrow .eq. 1) then
         do 10  i=1, ncol
            dsp = dsp - fact(icol(i))
   10    continue
         go to 9000
      end if
!
      if (ncol .eq. 1) then
         do 20  i=1, nrow
            dsp = dsp - fact(irow(i))
   20    continue
         go to 9000
      end if
!
      if (nrow*ncol .eq. 4) then
         if (irow(2) .le. icol(2)) then
            dsp = dsp - fact(irow(2)) - fact(icol(1)) - fact(icol(2)-irow(2))
         else
            dsp = dsp - fact(icol(2)) - fact(irow(1)) - fact(irow(2)-icol(2))
         end if
         go to 9000
      end if
! initialization before loop
      do 30  i=1, nrow
         irstk(i,1) = irow(nrow-i+1)
   30 continue
!
      do 40  j=1, ncol
         icstk(j,1) = icol(ncol-j+1)
   40 continue
!
      nro      = nrow
      nco      = ncol
      nrstk(1) = nro
      ncstk(1) = nco
      ystk(1)  = 0.0
      y        = 0.0
      istk     = 1
      l        = 1
      amx      = 0.0
!
   50 ir1 = irstk(1,istk)
      ic1 = icstk(1,istk)
      if (ir1 .gt. ic1) then
         if (nro .ge. nco) then
            m = nco - 1
            n = 2
         else
            m = nro
            n = 1
         end if
      else if (ir1 .lt. ic1) then
         if (nro .le. nco) then
            m = nro - 1
            n = 1
         else
            m = nco
            n = 2
         end if
      else
         if (nro .le. nco) then
            m = nro - 1
            n = 1
         else
            m = nco - 1
            n = 2
         end if
      end if
!
   60 if (n .eq. 1) then
         i = l
         j = 1
      else
         i = 1
         j = l
      end if
!
      irt = irstk(i,istk)
      ict = icstk(j,istk)
      mn  = irt
      if (mn .gt. ict) mn = ict
      y = y + fact(mn)
      if (irt .eq. ict) then
         nro = nro - 1
         nco = nco - 1
         call f11act (irstk(1,istk), i, nro, irstk(1,istk+1))
         call f11act (icstk(1,istk), j, nco, icstk(1,istk+1))
      else if (irt .gt. ict) then
         nco = nco - 1
         call f11act (icstk(1,istk), j, nco, icstk(1,istk+1))
         call f8xact (irstk(1,istk), irt-ict, i, nro, irstk(1,istk+1))
      else
         nro = nro - 1
         call f11act (irstk(1,istk), i, nro, irstk(1,istk+1))
         call f8xact (icstk(1,istk), ict-irt, j, nco, icstk(1,istk+1))
      end if
!
      if (nro .eq. 1) then
         do 70  k=1, nco
            y = y + fact(icstk(k,istk+1))
   70    continue
         go to 90
      end if
!
      if (nco .eq. 1) then
         do 80  k=1, nro
            y = y + fact(irstk(k,istk+1))
   80    continue
         go to 90
      end if
!
      lstk(istk)  = l
      mstk(istk)  = m
      nstk(istk)  = n
      istk        = istk + 1
      nrstk(istk) = nro
      ncstk(istk) = nco
      ystk(istk)  = y
      l           = 1
      go to 50
!
   90 if (y .gt. amx) then
         amx = y
         if (dsp-amx .le. tol) then
            dsp = 0.0
            go to 9000
         end if
      end if
!
  100 istk = istk - 1
      if (istk .eq. 0) then
         dsp = dsp - amx
         if (dsp-amx .le. tol) dsp = 0.0
         go to 9000
      end if
      l = lstk(istk) + 1
!
  110 if (l .gt. mstk(istk)) go to 100
      n   = nstk(istk)
      nro = nrstk(istk)
      nco = ncstk(istk)
      y   = ystk(istk)
      if (n .eq. 1) then
         if (irstk(l,istk) .lt. irstk(l-1,istk)) go to 60
      else if (n .eq. 2) then
         if (icstk(l,istk) .lt. icstk(l-1,istk)) go to 60
      end if
!
      l = l + 1
      go to 110
 9000 return
      end
!-----------------------------------------------------------------------
! Name: F5XACT
!
! Purpose: Put node on stack in network algorithm.
!
! Usage: CALL F5XACT (PASTP, TOL, KVAL, KEY, LDKEY, IPOIN, STP,
! LDSTP, IFRQ, NPOIN, NR, NL, IFREQ, ITOP,
! IPSH)
!
! Arguments:
! PASTP - The past path length. (Input)
! TOL - Tolerance for equivalence of past path lengths. (Input)
! KVAL - Key value. (Input)
! KEY - Vector of length LDKEY containing the key values.
! (Input/output)
! LDKEY - Length of vector KEY. (Input)
! IPOIN - Vector of length LDKEY pointing to the linked list
! of past path lengths. (Input/output)
! STP - Vector of length LSDTP containing the linked lists
! of past path lengths. (Input/output)
! LDSTP - Length of vector STP. (Input)
! IFRQ - Vector of length LDSTP containing the past path
! frequencies. (Input/output)
! NPOIN - Vector of length LDSTP containing the pointers to
! the next past path length. (Input/output)
! NR - Vector of length LDSTP containing the right object
! pointers in the tree of past path lengths.
! (Input/output)
! NL - Vector of length LDSTP containing the left object
! pointers in the tree of past path lengths.
! (Input/output)
! IFREQ - Frequency of the current path length. (Input)
! ITOP - Pointer to the top of STP. (Input)
! IPSH - Option parameter. (Input)
! If IPSH is true, the past path length is found in the
! table KEY. Otherwise the location of the past path
! length is assumed known and to have been found in
! a previous call.
!-----------------------------------------------------------------------
      subroutine f5xact (pastp, tol, kval, key, ldkey, ipoin, stp,         &
     &                   ldstp, ifrq, npoin, nr, nl, ifreq, itop, ipsh)
! SPECIFICATIONS FOR ARGUMENTS
      integer    kval, ldkey, ldstp, ifreq, itop, key(*), ipoin(*),        &
     &           ifrq(*), npoin(*), nr(*), nl(*)
      double precision pastp, tol, stp(*)
      logical    ipsh
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    ipn, ird, itmp
      double precision test1, test2
! SPECIFICATIONS FOR SAVE VARIABLES
      integer    itp
      save       itp
! SPECIFICATIONS FOR INTRINSICS
      intrinsic  mod
      integer    mod
! SPECIFICATIONS FOR SUBROUTINES
      external   prterr
!
      if (ipsh) then
! Convert KVAL to integer in range
! 1, ..., LDKEY.
         ird = mod(kval,ldkey) + 1
! Search for an unused location
         do 10  itp=ird, ldkey
            if (key(itp) .eq. kval) go to 40
            if (key(itp) .lt. 0) go to 30
   10    continue
!
         do 20  itp=1, ird - 1
            if (key(itp) .eq. kval) go to 40
            if (key(itp) .lt. 0) go to 30
   20    continue
! Return if KEY array is full
         call prterr(6, 'LDKEY is too small for this problem.  It is '//  &
     &               'not possible to estimate the value of LDKEY '//     &
     &               'required, but twice the current value may be '//    &
     &               'sufficient.')
! Update KEY
   30    key(itp) = kval
         itop       = itop + 1
         ipoin(itp) = itop
! Return if STP array full
         if (itop .gt. ldstp) then
            call prterr(7, 'LDSTP is too small for this problem.  It '//  &
     &                  'is not possible to estimate the value of '//     &
     &                  'LDSTP required, but twice the current value '//  &
     &                  'may be sufficient.')
         end if
! Update STP, etc.
         npoin(itop) = -1
         nr(itop)    = -1
         nl(itop)    = -1
         stp(itop)   = pastp
         ifrq(itop)  = ifreq
         go to 9000
      end if
! Find location, if any, of pastp
   40 ipn = ipoin(itp)
      test1 = pastp - tol
      test2 = pastp + tol
!
   50 if (stp(ipn) .lt. test1) then
         ipn = nl(ipn)
         if (ipn .gt. 0) go to 50
      else if (stp(ipn) .gt. test2) then
         ipn = nr(ipn)
         if (ipn .gt. 0) go to 50
      else
         ifrq(ipn) = ifrq(ipn) + ifreq
         go to 9000
      end if
! Return if STP array full
      itop = itop + 1
      if (itop .gt. ldstp) then
         call prterr(7, 'LDSTP is too small for this problem.  It is '// &
     &               'not possible to estimate the value of LDSTP '//    &
     &               'rerquired, but twice the current value may be '//  &
     &               'sufficient.')
         go to 9000
      end if
! Find location to add value
      ipn  = ipoin(itp)
      itmp = ipn
   60 if (stp(ipn) .lt. test1) then
         itmp = ipn
         ipn  = nl(ipn)
         if (ipn .gt. 0) then
            go to 60
         else
            nl(itmp) = itop
         end if
      else if (stp(ipn) .gt. test2) then
         itmp = ipn
         ipn  = nr(ipn)
         if (ipn .gt. 0) then
            go to 60
         else
            nr(itmp) = itop
         end if
      end if
! Update STP, etc.
      npoin(itop) = npoin(itmp)
      npoin(itmp) = itop
      stp(itop)   = pastp
      ifrq(itop)  = ifreq
      nl(itop)    = -1
      nr(itop)    = -1
!
 9000 return
      end
!-----------------------------------------------------------------------
! Name: F6XACT
!
! Purpose: Pop a node off the stack.
!
! Usage: CALL F6XACT (NROW, IROW, IFLAG, KYY, KEY, LDKEY, LAST,
! IPN)
!
! Arguments:
! NROW - The number of rows in the table. (Input)
! IROW - Vector of length nrow containing the row sums on output.
! (Output)
! IFLAG - Set to 3 if there are no additional nodes to process.
! (Output)
! KYY - Constant mutlipliers used in forming the hash table key.
! (Input)
! KEY - Vector of length LDKEY containing the hash table keys.
! (Input/output)
! LDKEY - Length of vector KEY. (Input)
! LAST - Index of the last key popped off the stack.
! (Input/output)
! IPN - Pointer to the linked list of past path lengths.
! (Output)
!-----------------------------------------------------------------------
      subroutine f6xact (nrow, irow, iflag, kyy, key, ldkey, last, ipn)
! SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, iflag, ldkey, last, ipn, irow(*), kyy(*), key(*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    j, kval
! SPECIFICATIONS FOR SAVE VARIABLES
!
   10 last = last + 1
      if (last .le. ldkey) then
         if (key(last) .lt. 0) go to 10
! Get KVAL from the stack
         kval      = key(last)
         key(last) = -9999
         do 20  j=nrow, 2, -1
            irow(j) = kval/kyy(j)
            kval    = kval - irow(j)*kyy(j)
   20    continue
         irow(1) = kval
         ipn     = last
      else
         last  = 0
         iflag = 3
      end if
      return
      end
!-----------------------------------------------------------------------
! Name: F7XACT
!
! Purpose: Generate the new nodes for given marinal totals.
!
! Usage: CALL F7XACT (NROW, IMAX, IDIF, K, KS, IFLAG)
!
! Arguments:
! NROW - The number of rows in the table. (Input)
! IMAX - The row marginal totals. (Input)
! IDIF - The column counts for the new column. (Input/output)
! K - Indicator for the row to decrement. (Input/output)
! KS - Indicator for the row to increment. (Input/output)
! IFLAG - Status indicator. (Output)
! If IFLAG is zero, a new table was generated. For
! IFLAG = 1, no additional tables could be generated.
!-----------------------------------------------------------------------
      subroutine f7xact (nrow, imax, idif, k, ks, iflag)
! SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, k, ks, iflag, imax(*), idif(*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, k1, m, mm
! SPECIFICATIONS FOR INTRINSICS
      intrinsic  min0
      integer    min0
!
      iflag = 0
! Find node which can be
! incremented, ks
      if (ks .eq. 0) then
   10    ks = ks + 1
         if (idif(ks) .eq. imax(ks)) go to 10
      end if
! Find node to decrement (>ks)
   20 if (idif(k).gt.0 .and. k.gt.ks) then
         idif(k) = idif(k) - 1
   30    k = k - 1
         if (imax(k) .eq. 0) go to 30
         m = k
! Find node to increment (>=ks)
   40    if (idif(m) .ge. imax(m)) then
            m = m - 1
            go to 40
         end if
         idif(m) = idif(m) + 1
! Change ks
         if (m .eq. ks) then
            if (idif(m) .eq. imax(m)) ks = k
         end if
      else
! Check for finish
   50    do 60  k1=k + 1, nrow
            if (idif(k1) .gt. 0) go to 70
   60    continue
         iflag = 1
         go to 9000
! Reallocate counts
   70    mm = 1
         do 80  i=1, k
            mm      = mm + idif(i)
            idif(i) = 0
   80    continue
         k = k1
   90    k = k - 1
         m       = min0(mm,imax(k))
         idif(k) = m
         mm      = mm - m
         if (mm.gt.0 .and. k.ne.1) go to 90
! Check that all counts
! reallocated
         if (mm .gt. 0) then
            if (k1 .ne. nrow) then
               k = k1
               go to 50
            end if
            iflag = 1
            go to 9000
         end if
! Get ks
         idif(k1) = idif(k1) - 1
         ks       = 0
  100    ks = ks + 1
         if (ks .gt. k) go to 9000
         if (idif(ks) .ge. imax(ks)) go to 100
      end if
!
 9000 return
      end
!-----------------------------------------------------------------------
! Name: F8XACT
!
! Purpose: Routine for reducing a vector when there is a zero
! element.
!
! Usage: CALL F8XACT (IROW, IS, I1, IZERO, NEW)
!
! Arguments:
! IROW - Vector containing the row counts. (Input)
! IS - Indicator. (Input)
! I1 - Indicator. (Input)
! IZERO - Position of the zero. (Input)
! NEW - Vector of new row counts. (Output)
!-----------------------------------------------------------------------
      subroutine f8xact (irow, is, i1, izero, new)
! SPECIFICATIONS FOR ARGUMENTS
      integer    is, i1, izero, irow(*), new(*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i
!
      do 10  i=1, i1 - 1
         new(i) = irow(i)
   10 continue
!
      do 20  i=i1, izero - 1
         if (is .ge. irow(i+1)) go to 30
         new(i) = irow(i+1)
   20 continue
!
      i = izero
   30 new(i) = is
   40 i = i + 1
      if (i .gt. izero) return
      new(i) = irow(i)
      go to 40
      end
!-----------------------------------------------------------------------
! Name: F9XACT
!
! Purpose: Computes the log of a multinomial coefficient.
!
! Usage: F9XACT(N, MM, IR, FACT)
!
! Arguments:
! N - Length of IR. (Input)
! MM - Number for factorial in numerator. (Input)
! IR - Vector of length N containing the numebers for the
! denominator of the factorial. (Input)
! FACT - Table of log factorials. (Input)
! F9XACT - The log of the multinomal coefficient. (Output)
!-----------------------------------------------------------------------
      double precision function f9xact (n, mm, ir, fact)
! SPECIFICATIONS FOR ARGUMENTS
      integer    n, mm, ir(*)
      double precision fact(0:*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    k
!
      f9xact = fact(mm)
      do 10  k=1, n
         f9xact = f9xact - fact(ir(k))
   10 continue
!
      return
      end
!-----------------------------------------------------------------------
! Name: F10ACT
!
! Purpose: Computes the shortest path length for special tables.
!
! Usage: CALL F10ACT (NROW, IROW, NCOL, ICOL, VAL, XMIN, FACT, ND,
! NE, M)
!
! Arguments:
! NROW - The number of rows in the table. (Input)
! IROW - Vector of length NROW containing the row totals. (Input)
! NCOL - The number of columns in the table. (Input)
! ICO - Vector of length NCOL containing the column totals.
! (Input)
! VAL - The shortest path. (Output)
! XMIN - Set to true if shortest path obtained. (Output)
! FACT - Vector containing the logarithms of factorials.
! (Input)
! ND - Workspace vector of length NROW.
! NE - Workspace vector of length NCOL.
! M - Workspace vector of length NCOL.
!
! Chapter: STAT/LIBRARY Categorical and Discrete Data Analysis
!-----------------------------------------------------------------------
      subroutine f10act (nrow, irow, ncol, icol, val, xmin, fact, nd, ne, m)
! SPECIFICATIONS FOR ARGUMENTS
      integer    nrow, ncol, irow(*), icol(*), nd(*), ne(*), m(*)
      double precision val, fact(0:*)
      logical    xmin
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, is, ix, nrw1
!
      do 10  i=1, nrow - 1
         nd(i) = 0
   10 continue
!
      is    = icol(1)/nrow
      ne(1) = is
      ix    = icol(1) - nrow*is
      m(1)  = ix
      if (ix .ne. 0) nd(ix) = nd(ix) + 1
!
      do 20  i=2, ncol
         ix    = icol(i)/nrow
         ne(i) = ix
         is    = is + ix
         ix    = icol(i) - nrow*ix
         m(i)  = ix
         if (ix .ne. 0) nd(ix) = nd(ix) + 1
   20 continue
!
      do 30  i=nrow - 2, 1, -1
         nd(i) = nd(i) + nd(i+1)
   30 continue
!
      ix   = 0
      nrw1 = nrow + 1
      do 40  i=nrow, 2, -1
         ix = ix + is + nd(nrw1-i) - irow(i)
         if (ix .lt. 0) return
   40 continue
!
      do 50  i=1, ncol
         ix  = ne(i)
         is  = m(i)
         val = val + is*fact(ix+1) + (nrow-is)*fact(ix)
   50 continue
      xmin = .true.
!
      return
      end
!-----------------------------------------------------------------------
! Name: F11ACT
!
! Purpose: Routine for revising row totals.
!
! Usage: CALL F11ACT (IROW, I1, I2, NEW)
!
! Arguments:
! IROW - Vector containing the row totals. (Input)
! I1 - Indicator. (Input)
! I2 - Indicator. (Input)
! NEW - Vector containing the row totals. (Input)
!-----------------------------------------------------------------------
      subroutine f11act (irow, i1, i2, new)
! SPECIFICATIONS FOR ARGUMENTS
      integer    i1, i2, irow(*), new(*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i
!
      do 10  i=1, i1 - 1
         new(i) = irow(i)
   10 continue
!
      do 20  i=i1, i2
         new(i) = irow(i+1)
   20 continue
!
      return
      end
!-----------------------------------------------------------------------
! Name: ERPRT
!
! Purpose: Print an error message and stop.
!
! Usage: CALL ERPRT (ICODE, MES)
!
! Arguments:
! ICODE - Integer code for the error message. (Input)
! MES - Character string containing the error message. (Input)
!-----------------------------------------------------------------------
      subroutine prterr (icode, mes)
! SPECIFICATIONS FOR ARGUMENTS
      integer    icode
      character  mes*(*)
!
      write (*,*) 'FEXACT ERROR: ', icode, ' ', mes
      stop
      end
!-----------------------------------------------------------------------
! Name: IWORK
!
! Purpose: Routine for allocating workspace.
!
! Usage: IWORK (IWKMAX, IWKPT, NUMBER, ITYPE)
!
! Arguments:
! IWKMAX - Maximum length of workspace. (Input)
! IWKPT - Amount of workspace currently allocated. (Input/output)
! NUMBER - Number of elements of workspace desired. (Input)
! ITYPE - Worspace type. (Input)
! ITYPE TYPE
! 2 Integer
! 3 Real
! 4 Double Precision
! IWORK - Index in RWRK, DWRK, or IWRK of the beginning of the
! first element in the workspace array. (Output)
!-----------------------------------------------------------------------
      integer function iwork (iwkmax, iwkpt, number, itype)
! SPECIFICATIONS FOR ARGUMENTS
      integer    iwkmax, iwkpt, number, itype
! SPECIFICATIONS FOR INTRINSICS
      intrinsic  mod
      integer    mod
! SPECIFICATIONS FOR SUBROUTINES
      external   prterr
!
      iwork = iwkpt
      if (itype.eq.2 .or. itype.eq.3) then
         iwkpt = iwkpt + number
      else
         if (mod(iwork,2) .ne. 0) iwork = iwork + 1
         iwkpt = iwkpt + 2*number
         iwork = iwork/2
      end if
      if (iwkpt .gt. iwkmax+1) then
         call prterr (40, 'Out of workspace.')
      end if
      return
      end
!-----------------------------------------------------------------------
! Name: ISORT
!
! Purpose: Shell sort for an integer vector.
!
! Usage: CALL ISORT (N, IX)
!
! Arguments:
! N - Lenth of vector IX. (Input)
! IX - Vector to be sorted. (Input/output)
!-----------------------------------------------------------------------
      subroutine isort (n, ix)
! SPECIFICATIONS FOR ARGUMENTS
      integer    n, ix(*)
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    i, ikey, il(10), it, iu(10), j, kl, ku, m
! SPECIFICATIONS FOR SUBROUTINES
      external   prterr
! Sort IX
      m = 1
      i = 1
      j = n
   10 if (i .ge. j) go to 40
      kl   = i
      ku   = j
      ikey = i
      j    = j + 1
! Find element in first half
   20 i = i + 1
      if (i .lt. j) then
         if (ix(ikey) .gt. ix(i)) go to 20
      end if
! Find element in second half
   30 j = j - 1
      if (ix(j) .gt. ix(ikey)) go to 30
! Interchange
      if (i .lt. j) then
         it    = ix(i)
         ix(i) = ix(j)
         ix(j) = it
         go to 20
      end if
      it       = ix(ikey)
      ix(ikey) = ix(j)
      ix(j)    = it
! Save upper and lower subscripts of
! the array yet to be sorted
      if (m .lt. 11) then
         if (j-kl .lt. ku-j) then
            il(m) = j + 1
            iu(m) = ku
            i     = kl
            j     = j - 1
         else
            il(m) = kl
            iu(m) = j - 1
            i     = j + 1
            j     = ku
         end if
         m = m + 1
         go to 10
      else
         call prterr (20, 'This should never occur.')
      end if
! Use another segment
   40 m = m - 1
      if (m .eq. 0) go to 9000
      i = il(m)
      j = iu(m)
      go to 10
!
 9000 return
      end
!-----------------------------------------------------------------------
! Name: GAMMDS
!
! Purpose: Cumulative distribution for the gamma distribution.
!
! Usage: PGAMMA (Q, ALPHA,IFAULT)
!
! Arguments:
! Q - Value at which the distribution is desired. (Input)
! ALPHA - Parameter in the gamma distribution. (Input)
! IFAULT - Error indicator. (Output)
! IFAULT DEFINITION
! 0 No error
! 1 An argument is misspecified.
! 2 A numerical error has occurred.
! PGAMMA - The cdf for the gamma distribution with parameter alpha
! evaluated at Q. (Output)
!-----------------------------------------------------------------------
!
! Algorithm AS 147 APPL. Statist. (1980) VOL. 29, P. 113
!
! Computes the incomplete gamma integral for positive
! parameters Y, P using and infinite series.
!
      double precision function gammds (y, p, ifault)
! SPECIFICATIONS FOR ARGUMENTS
      integer    ifault
      double precision y, p
! SPECIFICATIONS FOR LOCAL VARIABLES
      integer    ifail
      double precision a, c, f
! SPECIFICATIONS FOR SAVE VARIABLES
      double precision e, one, zero
      save       e, one, zero
! SPECIFICATIONS FOR INTRINSICS
      intrinsic  dlog, dexp
      double precision dlog, dexp
! SPECIFICATIONS FOR FUNCTIONS
      external   alogam
      double precision alogam
      double precision zexp, zlog
!
      data e, zero, one/1.0d-6, 0.0d0, 1.0d0/
!
      zexp(a) = dexp(a)
      zlog(a) = dlog(a)
!
! Checks for the admissibility of arguments and value of F
!
      ifault = 1
      gammds = zero
      if (y.le.zero .or. p.le.zero) return
      ifault = 2
!
! ALOGAM is natural log of gamma function
! no need to test ifail as an error is impossible
!
      f = zexp(p*zlog(y)-alogam(p+one,ifail)-y)
      if (f .eq. zero) return
      ifault = 0
!
! Series begins
!
! = one
      gammds = one
      a      = p
   10 a = a + one
! = c*y/a
      gammds = gammds + c
      if (c/gammds .gt. e) go to 10
      gammds = gammds*f
      return
      end
!-----------------------------------------------------------------------
! Name: ALOGAM
!
! Purpose: Value of the log-gamma function.
!
! Usage: ALOGAM (X, IFAULT)
!
! Arguments:
! X - Value at which the log-gamma function is to be evaluated.
! (Input)
! IFAULT - Error indicator. (Output)
! IFAULT DEFINITION
! 0 No error
! 1 X .LT. 0
! ALGAMA - The value of the log-gamma function at XX. (Output)
!-----------------------------------------------------------------------
!
! Algorithm ACM 291, Comm. ACM. (1966) Vol. 9, P. 684
!
! Evaluates natural logarithm of gamma(x)
! for X greater than zero.
!
      double precision function alogam (x, ifault)
! SPECIFICATIONS FOR ARGUMENTS
      integer    ifault
      double precision x
! SPECIFICATIONS FOR LOCAL VARIABLES
      double precision f, y, z
! SPECIFICATIONS FOR SAVE VARIABLES
      double precision a1, a2, a3, a4, a5, half, one, seven, zero
      save       a1, a2, a3, a4, a5, half, one, seven, zero
! SPECIFICATIONS FOR INTRINSICS
      intrinsic  dlog
      double precision dlog
      double precision zlog
!
! The following constants are dlog(2PI)/2,
! half, zero, one, seven
!
      data a1, a2, a3, a4, a5/0.918938533204673d0, 0.000595238095238d0,     &
     &     0.000793650793651d0, 0.002777777777778d0,                        &
     &     0.083333333333333d0/
      data half, zero, one, seven/0.5d0, 0.0d0, 1.0d0, 7.0d0/
!
      zlog(f) = dlog(f)
!
      alogam = zero
      ifault = 1
      if (x .lt. zero) return
      ifault = 0
      y      = x
      f      = zero
      if (y .ge. seven) go to 30
      f = y
   10 y = y + one
      if (y .ge. seven) go to 20
      f = f*y
      go to 10
   20 f = -zlog(f)
   30 z = one/(y*y)
      alogam = f + (y-half)*zlog(y) - y + a1 + (((-a2*z+a3)*z-a4)*z+a5) &
    &          /y
      return
      end
!
!! END MODULE FISHER_EXACT
