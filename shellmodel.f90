! Henry Jared Doster
! PHY 981: Nuclear Structure
!
! Configuration-interaction method (shell-model diagonalization)
!
! The purpose of this code is to construct a many-body Hamiltonian matrix in a basis of Slater determinants 
! and then to find the low-lying energy eigenstates.


program diagonalize

  use f90library

  implicit none

  integer :: h,i,j,k, count                       ! Loop variables
  integer :: total                                ! count number of unbroken pairs of particles
  integer :: p,q,r,s                              ! Hamiltonian loop variables
  integer :: pos1, pos2, phase, coef, matelem     ! ladder operator variables
  integer :: NSP                                  ! Number of single-particle states (read from text file)
  integer :: NSD                                  ! Maximum number of slater determinants
  integer :: NUM                                  ! Number of particles
  integer :: mproj
  real(8) :: g                                    ! pairing contribution strengh
  LOGICAL :: test                                 ! test for goto statement (goto can't be in a subroutine)


  real(8), allocatable :: d(:), e(:), z(:,:)      ! arguments for diagonalization functions in f90lib module
  integer :: dim

  integer, allocatable :: ind(:), n(:), l(:), j2(:), jz2(:), T(:) ! Quantum numbers read from text file
  integer, allocatable :: M(:)
  integer,allocatable :: sd(:,:), occ(:)          ! Slater determinant (num of slater dets, num of particles/occupation numbers)
  integer,allocatable :: msd(:,:), psd(:,:), store(:,:)
  real(8), allocatable :: Ham(:,:), Ham1(:,:), Ham2(:,:), Hamm(:,:)   ! The Hamiltonians
  integer, allocatable :: bra(:,:), ket(:,:)      ! bras and kets
  real,allocatable :: energy(:)                   ! Energy values read from text file

  character(len=20) :: filename


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Read data from text files into allocatable arrays !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print*,
  print*, 'Enter name of file with valence space data (including the extension).'
  read*, filename
  print*,

  call opentextfiles(filename)

  read(10,*) NSP

  allocate( ind(NSP), n(NSP), l(NSP), j2(NSP), jz2(NSP), T(NSP), energy(NSP) )

  do i=1,NSP
     read(10,*) ind(i), n(i), l(i), j2(i), jz2(i), T(i), energy(i)
  enddo

!! Print contents of text file
!  do i=1,NSP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     print*, ind(i), n(i), l(i), j2(i), jz2(i), T(i), energy(i)
!  enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Get information from user and allocate remaining arrays !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print*, "There are", NSP, "possible single-particle state in the given valence space."
  print*,
  print*, "How many particles are in the valence space?"
  read(*,*) NUM
  print*,

!! Calculate maximum number of slater determinants and allocate arrays
  NSD = Factorial(NSP) / (Factorial(NUM)*Factorial(NSP-NUM))
  print*, NSP, "choose", NUM, "is", NSD
  print*,

  allocate( sd(NSD,NUM), occ(NUM), M(NSD) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Enumerate all possible Slater determinants !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Initialize the slater determinant array to be in lowest energy 1,2,3,4...
  do i=1,NUM
     occ(i) = i
  enddo

  sd(1,:) = occ


!! Run odometer 
  do i=1,NSD
     do j = NUM,1,-1
        if(occ(j) < NSP-NUM+j )then
           h = occ(j)
           do k = j,NUM
              occ(k) = h+1+k-j ! sequential order                                                                                                                       
           end do
           sd(i+1,:) = occ                                                                                                                                       
           exit
        end if
     end do
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Select Slater determinants of a particular total projection M !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Calculate M for each state

  do i=1,NSD
    count = 0
    do j=1,NUM
        count = count + jz2(sd(i,j))
     enddo
     M(i) = count
  enddo


!  print*, 'All possible projections of Slater determinants'
!  print*,
!  print*, 'Slater determinants', 'Projections'
!  do i=1, NSD
!     print*, sd(i,:), M(i) 
!  enddo
!  print*,

  print*, "Please specify the desired total projection (2*M)."
  read(*,*) mproj
  print*,

!! Count number of states with "mproj"

  count=0
  do i=1,NSD
     if ( M(i) == mproj ) then
        count = count +1
     endif
  enddo


!! Create new array and reassign desired Slater determinants

  allocate( msd(count, NUM) )

  j=1
  do i=1,NSD
     if ( M(i) == mproj ) then
        msd(j,:) = sd(i,:)
        j = j+1
     endif
  enddo
 
  print*, 'There are', count, 'Slater determinants with the total projection', mproj
  print*,

!  do i=1,count
!     print*, msd(i,:)
!  enddo
!  print*,


!! If mproj=0, find the Slater determinants with unbroken pairs
!! REQUIRES: "count" from previous step

  if (mproj==0 .AND. MOD(NUM,2)==0) then ! Determine if projection=0 and there are an even number of particles

     !! Allocate array that is definitely big enough to store all possible states with unbroken pairs
     allocate ( store(count,NUM) )


     !! Store states with unbroken pairs
     total=0
     do i=1,count                        ! Choose a slater determinant with M=0

        k=0
        do j=1, (NUM/2+1), 2             ! Choose a pair of particles
           if ( MOD( msd(i,j),2) /= 0 .AND. msd(i,j)+1==msd(i,j+1)) then   ! check if particles are paired
              k=k+1                                                        ! count the number of unbroken pairs
           endif
        enddo
        
        if (k == NUM/2) then             ! Determine if all pairs in this state are unbroken
           total = total +1              ! Count the number of states with all pair unbroken
           store(total,:) = msd(i,:)
        endif
     enddo


     !! Allocate array of correct size and move states to it
     allocate( psd(total,NUM) )

     do i=1,total
        psd(i,:) = store(i,:)
     enddo


     !! Print results
     print*, 'There are', total, 'Slater determinants with unbroken pairs'
     print*,

     do i=1,total
        print*, psd(i,:)
     enddo
     print*,

  endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Calculate Hamiltonian !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Allocate Hamiltonian and initialize to zero
! For this specific case, the Hamiltonian has dimension equal to the number of SD with unbroken pairs

  allocate(Ham(total,total), Ham1(total,total), Ham2(total,total), Hamm(total,total), bra(total,NUM), ket(total,NUM))
  Ham = 0
  Ham1 = 0
  Ham2 = 0


! We are calculating matrix elements for states with unbroken pairs only
! assigned SD with unbroken pairs to the bras and kets 

  do i=1,total  
     bra(i,:) = psd(i,:)
     ket(i,:) = bra(i,:)
  enddo

! In the following do loops, the loop is skipped if the ladder operator returns 0
! becuase skipping the loop is effectively the same as making it zero.


!!!!!!!!!!!!!!!!!!!!!!! One-Body Hamiltonian !!!!!!!!!!!!!!!!!!!!!!!!

! Each single-particle index "k" from the text file has a unique set of
! single-particle quantum numbers


  do j=1,total          ! Choose a column (ket state) of the Hamiltonian matrix

     do k=1,NSP       ! Sum over all possible many body indicies

        ! Initializations
        ket(j,:) = bra(j,:)  ! re-initilize ket state
        phase = 0
        count = 0       ! Phase counter
        pos1=0          ! position of first lowering operator
        test = .TRUE.
                 
        ! Calculate one substate of the ket state
        call lower(ket(j,:), k, pos1, count, test)
        if (test .eqv. .FALSE.) then
!!           print*,'lower cycle'
!!           print*,'-----------------------------------------------------------'
           cycle
        endif

        call raise(ket(j,:), k, pos1)

        call normal_order(ket(j,:), count, phase)
        
        do i=1,total        ! Sum over rows (bra states)
           call orthonormality(bra(j,:), ket(i,:), coef)
           Ham1(i,j) = Ham1(i,j) + (n(k)-1)*coef*phase      ! add results to previous substate calculation
        enddo
                 
     
     enddo

!     print*, 'final ket =', ket(j,:)

!     print*, 'end loop'
!     print*,'--------------------------------------------------'
!     print*,

  enddo

  

! Print final Hamiltonian

  print*, 'The One-Body Hamiltonian H_0 is'
  print*,
  do i=1,total
     print'(10(F10.2))', Ham1(i,:)
  enddo
  print*,

!!!!!!!!!!!!!!!!!!!!!!! Two-Body Hamiltonian !!!!!!!!!!!!!!!!!!!!!!!!

! Potential V is zero for all states except those where pairing is unbroken
! This makes calculating the matrix elements an easier task.

! ususally, we would have to read all possible combinations of slater determinants
! into a function that calculates the matrix elements.
! But in this case, it's always "g" or "0".  So, just check the matrix elements for
! the constraints to see if it is 0.

! FORTRAN is a row-major order programming language (i=row, j=column)


  do j=1,total          ! Choose a column (ket state) of the Hamiltonian matrix
     
     do p=1,NSP         ! Sum over ALL possible ket substates
        do q=p+1,NSP                 
           do r=1,NSP
              do s=r+1,NSP
                 
!!                 read*,
!!                 print*,'p=', p
!!                 print*,'q=', q
!!                 print*,'r=', r
!!                 print*,'s=', s
!!                 print*,
                 

!!! UNBROKEN PAIRS MATRIX ELEMENT CONSTRAINTS: principle quantum number, spin opposites

                 ! Check that bra state satisfies matrix element constraints FOR UNBROKEN PAIRS
                 if ( n(p) /= n(q) ) then
!!                    print*,'first matrix element cycle'
!!                    print*,'-----------------------------------------------------------'
                    cycle
                 elseif ( jz2(p) /= -jz2(q) ) then
!!                    print*,'second matrix element cycle'
!!                    print*,'-----------------------------------------------------------'
                    cycle
                 endif
                 
                ! Check that ket substate satisfies matrix element constraints FOR UNBROKEN PAIRS
                 if ( n(r) /= n(s) ) then
!!                    print*,'third matrix element cycle'
!!                    print*,'-----------------------------------------------------------'
                    cycle
                 elseif ( jz2(r) /= -jz2(s) ) then
!!                    print*,'fourth matrix element cycle'
!!                    print*,'-----------------------------------------------------------'
                    cycle
                 endif
                 
                 ! Check spin contraint in matrix elements between bra and ket state       
                 if ( jz2(p) /= jz2(r) ) then
!!                    print*,'fifth matrix element cycle'
!!                    print*,'-----------------------------------------------------------'
                    cycle
                 elseif ( jz2(q) /= jz2(s) ) then
!!                    print*,'six matrix element cycle'
!!                    print*,'-----------------------------------------------------------'
                    cycle
                 endif
              
              
              
!!! CALCULATE LADDER OPERATORS
                 
                 ! Initializations
                 ket(j,:) = bra(j,:)  ! re-initilize ket state
                 phase=0
                 count = 0       ! Phase counter
                 pos1=0          ! position of first lowering operator
                 pos2=0          ! position of second lowering operator
                 test = .TRUE.
                 
                 ! Calculate one substate of the ket state
                 call lower(ket(j,:), s, pos1, count, test)
                 if (test .eqv. .FALSE.) then
!!                    print*,'first lower cycle'
!!                    print*,'-----------------------------------------------------------'
                    cycle
                 endif
                 
                 call lower(ket(j,:), r, pos2, count, test)
                 if (test .eqv. .FALSE.) then
!!                    print*,'second lower cycle'
!!                    print*,'-----------------------------------------------------------'
                    cycle
                 endif
                 
                 call raise(ket(j,:), q, pos2)
                 
                 call raise(ket(j,:), p, pos1)
                 
                 call normal_order(ket(j,:), count, phase)
                 
!!                 print*, 'initial bra = ', bra(j,:)
!!                 print*, 'final ket =', ket(j,:)
!!                 print*, 'phase =', phase
!!                 print*,
!!                 read*,
!!                 print*, 'start orthogonality'
                 
                 do i=1,total        ! Sum over rows (bra states)
                    call orthonormality(bra(i,:), ket(j,:), coef)
                    Ham2(i,j) = Ham2(i,j) + coef*phase        ! add results to previous substate calculation
                 enddo

!!                 print*, 'end loop'
!!                 print*,'--------------------------------------------------'
!!                 print*,
                 
              enddo
           enddo
        enddo
     enddo
     
  enddo
  

! Divide by 2 (there is double counting somewhere that I can't find)

  Ham2 = Ham2/2


! Print Hamiltonian

  print*, 'The two-Body Hamiltonian H_1 is g*h where h is'
  print*,
  do i=1,total
     print'(10(F10.2))', Ham2(i,:)
  enddo
  print*,

  

  print*, "What is the value of the pairing contribution g strength?"
  read*, g
  print*, g
  print*,

  Hamm = Ham2*g

  Ham = Ham1 + Hamm

  print*, "So, the total Hamiltonian H_0 + H_1 is"
  print*,
  do i=1,total
     print'(10(F10.2))', Ham(i,:)
  enddo
  print*,


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Diagonalize Hamiltonian !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dim = total
  allocate( d(dim), e(dim), z(dim,dim) )


  call tred2(Ham, dim, d, e)

  call tqli(d,e,dim,z)    ! d contains the eigenvalues

  print*, 'The eigenvalues of the Hamiltonian are'
  print '(10(F10.5))', d
  print*,


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Plot different values of g !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i = -50,50
     g = i/10

     Hamm = Ham2*g
     Ham = Ham1 + Hamm

     call tred2(Ham, dim, d, e)
     call tqli(d,e,dim,z)  

     WRITE(11,'(10(F10.2))') g, d(1), d(2), d(3), d(4), d(5), d(6)
     print '(10(F10.2))', g, d(1), d(2), d(3), d(4), d(5), d(6)
  enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Clean up !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call closetextfiles

  deallocate (ind,n,l,j2,jz2,T,energy,msd,store,psd)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Subroutines  !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  subroutine opentextfiles(filename)
    integer :: OPEN_STATUS
    character(len=20),intent(in) :: filename

    OPEN(UNIT=10,FILE=filename,STATUS="OLD",ACTION="READ",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, file not opened properly------------"
    endif

    OPEN(UNIT=11,FILE='data_for_plot.txt',STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, lambda1 not opened properly------------"
    endif

  end subroutine opentextfiles

  subroutine closetextfiles
    CLOSE(UNIT=10)

    CLOSE(UNIT=11)
  end subroutine closetextfiles


  INTEGER FUNCTION  Factorial(n)
    implicit none
    integer, intent(in) :: n
    integer :: i, Ans

    Ans = 1
    do i=1,n
       Ans = Ans*i
    enddo
    Factorial = Ans
  END FUNCTION Factorial



  SUBROUTINE bubble_sort(a, count)
    integer, intent(inout) :: count
    integer, intent(inout), dimension(:) :: a
    integer :: temp
    INTEGER :: i, j
    LOGICAL :: swapped = .TRUE.    !!!!!!!!! does this work??

    DO j = size(a)-1, 1, -1
       swapped = .FALSE.
       DO i = 1, j
          IF (a(i) > a(i+1)) THEN
             temp = a(i)
             a(i) = a(i+1)
             a(i+1) = temp
             swapped = .TRUE.
             count = count + 1    ! count all swops                                                                                                                              
          endif 
       END DO
       IF (.NOT. swapped) then   !!!!!!!!!!!! What does this do???
          return  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IS THIS RIGHT???
       endif
    END DO
  END SUBROUTINE bubble_sort


  subroutine normal_order(ket, count, phase)
    integer, intent(inout) :: count, phase
    integer, intent(inout), dimension(:) :: ket
    
    ! Perform normal ordering
    call bubble_sort(ket, count)

    if ( mod(count,2) == 0 ) then
       phase = 1
    elseif (mod(count,2) /=0 ) then
       phase = -1
    else
       print*, 'A mistake has been made.'
    endif
  endsubroutine normal_order




! Position is recorded by lower and used by raise

  subroutine lower(ket, index, position, count, test)
    integer, intent(in) :: index
    integer, intent(inout), dimension(:) :: ket
    integer, intent(inout) :: count, position
    logical, intent(inout) :: test
    integer :: check, location

    check=0

    ! Check for existance
    do i=1,size(ket)                  ! Run through all indicies of the SD
       if ( ket(i) == index ) then    ! If the index exists in the SD, then lower to zero
          ket(i) = 0
          position = i
          check = check + 1
       endif
    enddo

    if (check==0) then                  ! If the index does not exist, exit
       test = .FALSE.
!!!       print*, 'This index does not exist in the many-body wavefunction'
    endif
  endsubroutine lower


! Raising operator is defined such that the lowering operators must occur first
! So, operators must be in reduced, normal form for this function to work

  subroutine raise(ket, index, position)
    integer, intent(in) :: index
    integer, intent(inout), dimension(:) :: ket
    integer, intent(inout) :: position

    ! Check for non-zero
    if (ket(position) /= 0) then
       print*, 'This index cannot be raised.'
    endif

    ! Perform raising operation
    ket(position) = index
  end subroutine raise


  subroutine orthonormality(bra, ket, coef)
    integer, intent(in), dimension(:) :: bra,ket
    integer, intent(inout) :: coef
    integer :: i
    integer :: test

    test=0
    
    do i=1,size(bra)
       if (bra(i) == ket(i)) then
          test = test + 1
       endif
    enddo

    if( test == size(ket) ) then
       coef=1
    else
       coef=0
    endif

!!!    print*, '(', bra, ')', '      ', '(', ket,')', '=', coef

  end subroutine orthonormality


end program diagonalize
