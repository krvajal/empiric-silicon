program psudo
    use types
    use constants
    use utils
    implicit none

    integer,parameter :: BasisDimension= 113, NumKPoints = 100
    integer :: GPoints(BasisDimension,3)
    real(dp),parameter :: LatticeConstant = 10.320606615 ! 10.26121_dp ! in  bhor radius

    real(dp) :: KPoints(NumKPoints,3)
    real(dp) :: Hamiltonian(BasisDimension, BasisDimension)
    real(dp) :: EValues(BasisDimension)
    real(dp) :: offset = 0
    integer :: BANDS_UNIT
    real(dp) ::  RhoPoint(3), XPoint(3), WPoint(3), KPoint(3), LPoint(3)


    open(newunit=BANDS_UNIT, file="bands.txt", status="replace")

    !points in the reciprocal unit cell
    RhoPoint  = (/ 0.0_dp, 0.0_dp, 0.0_dp/); 
    XPoint = (/1.0_dp, 0.0_dp, 0.0_dp/);
    WPoint = (/ 1.0_dp, 0.5_dp, 0.0_dp/);
    LPoint  = 0.5_dp;
    KPoint  =   3.0_dp/4; KPoint(3) = 0

    call GenerateGValues

    call ComputeBands(RhoPoint,XPoint, offset)
    
    call ComputeBands(XPoint,WPoint, offset)
    ! print *, "offset", offset
    call ComputeBands(WPoint,LPoint, offset)
    ! print *, "offset", offset
    call ComputeBands(LPoint,RhoPoint, offset)
    ! print *, "offset", offset
    call ComputeBands(RhoPoint,KPoint, offset)
    print *, "offset", offset



    contains

    subroutine GenerateHamiltonian(k)
        implicit none
        integer :: i, j, l
        real(dp) :: k(3)
        integer :: DeltaG(3), ModDeltaGSquared
        real(dp) :: coeff  =  (2 * pi) /LatticeConstant

        do i  = 1, BasisDimension
            do j=1, BasisDimension
                forall(l=1:3)DeltaG(l) = GPoints(i,l) - GPoints(j,l)
                
                ModDeltaGSquared = sum(DeltaG**2)

                Hamiltonian(i,j) = cos(real(sum(DeltaG)) * pi/4.0_dp ) * SiliconPsudopot(ModDeltaGSquared)
                if(i==j)then
                    Hamiltonian(i,i) =  0.5*(norm2(k + GPoints(i,:))**2) *coeff*coeff
                endif
            enddo

        enddo


    end subroutine GenerateHamiltonian

    !=====================================
    ! Return the values of the silicon
    ! psudopotential
    !=====================================
    real(dp) function SiliconPsudopot(ModGSquared)
        implicit none
        integer :: ModGSquared
    
        if (  ModGSquared == 3) then
            SiliconPsudopot = -0.2242_dp/2
        else if ( ModGSquared == 8)then
            SiliconPsudopot =  0.0552_dp/2
        else if ( ModGSquared == 11) then
            SiliconPsudopot = 0.0724_dp/2
        else
            SiliconPsudopot = 0.0
            
        endif
    end function SiliconPsudopot

    subroutine GenerateGValues()
        implicit none
        integer :: n1,n2,n3
        real(dp) :: ModGSquared, GvalueCutoff
        integer :: count
        count = 0
        GvalueCutoff = 22
        do n1 = -10,10
            do n2 = -10,10
                do n3 = -10,10
                    ModGSquared = 3 * n1**2  + 3* n2**2 + 3 * n3**2 - 2*n1*n2 - 2 * n1*n3 - 2 *n2*n3
                    if (ModGSquared < GvalueCutoff) then
                        count = count + 1
                        GPoints(count,:) = [ (-n1 + n2 + n3), (n1 - n2 + n3),( n1 + n2 - n3) ]
                    endif
                enddo
            enddo
        enddo
        print *,"Basis Size:",count


    end subroutine GenerateGValues


    subroutine GenerateKPoints(N,p,pmin,pmax)
        implicit none
        integer, intent(in)     :: N
        real(dp),intent(out) :: p(N,3)
        integer:: i
        real(dp) :: spacing
        real(dp),intent(in) :: pmax(3)
        real(dp),intent(in) :: pmin(3)
        real(dp) :: len
        real(dp) :: direct_vec(3)
        direct_vec  = pmax - pmin
        len = sqrt( sum( (direct_vec)**2))

        direct_vec = direct_vec/len !normalize vector

        spacing = len/(N-1)

        do i = 0,N-1
            p(i+1,:)= pmin + i * spacing*direct_vec
        enddo

    end subroutine GenerateKPoints

    ! ============================================
    ! Solves the eigenvaule equation Ax = E x
    ! Input:
    !            Dim: positive integer with the dimension of A
    !            A: real symetric matrix to be diagonalized
    ! Ouput:
    !            Values: vector with the eigenvalues of A sorted asc
    ! ============================================
    subroutine SymetricEigenvalues(Dim,A, Values)
        integer :: Dim
        real(dp),intent(in) :: A(Dim,Dim)
        real(dp) :: A_(Dim,Dim)
        real(dp),intent(out) :: Values(Dim)
        integer,parameter :: LWORK_MAX = 1000
        real(dp) :: Work(LWORK_MAX)
        integer :: lWork, info
        A_ = A
        !query optimal lwork
        lWork = -1
        
        call dsyev("N", "L",Dim,A_,Dim,Values,Work,lWork,info)
        lWork = min(LWORK_MAX, int(Work(1)))
        
        !solve eigenproblem
        call dsyev("N", "U", Dim, A_, Dim, Values,Work,lWork,info)
        
        call assert(info == 0)


    end subroutine SymetricEigenvalues


    subroutine ComputeBands(StartPoint,EndPoint, offset)
        implicit none 
        real(dp),intent(in) :: StartPoint(3), EndPoint(3)
        real(dp), intent(inout) :: offset
        real(dp) ::KPoint(3)
        integer :: i

        

        call GenerateKPoints(NumKPoints, KPoints, StartPoint, EndPoint )

        do i = 1, NumKPoints
            KPoint  = KPoints(i,:)
            call GenerateHamiltonian(KPoint)
            call symetric_eigenvalues(BasisDimension, Hamiltonian, EValues)
            print *, KPoint
            write (BANDS_UNIT,*) norm2(KPoint - StartPoint) + offset, EValues
        enddo
        offset = offset + norm2(EndPoint- StartPoint,1)

    end subroutine ComputeBands

end program psudo