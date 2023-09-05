! Shell Lab code. Written by Jacob Monroe I believe.
! Called CosAngle3, but really returns angle given three coordinate positions
real*8 function CosAngle3(Pos1, Pos2, Pos3)
  implicit none
  real*8, parameter :: pi = 3.1415926535897931D0
  real*8, parameter :: DegPerRad = 180.D0/pi
  real*8, dimension(3), intent(in) :: Pos1, Pos2, Pos3
  real*8, dimension(3) :: Vec21, Vec23
  real*8 :: Norm, Phi
  if (all(Pos1 == Pos2) .or. all(Pos2 == Pos3)) then
    CosAngle3 = 0.
    return
  endif
  Vec21 = Pos1 - Pos2
  Vec23 = Pos3 - Pos2
  Norm = sqrt(sum(Vec21*Vec21)*sum(Vec23*Vec23))
  !CosAngle3 = dot_product(Vec21, Vec23) / Norm
  Phi = min(1.D0, max(-1.D0, dot_product(Vec21, Vec23) / Norm))
  Phi = acos(Phi)
  CosAngle3 = mod(Phi + pi, pi*2.D0) - pi
  if (CosAngle3 < -pi) CosAngle3 = CosAngle3 + pi*2.D0
  CosAngle3 = CosAngle3 * DegPerRad
end function

! Finds nearest neighbors in list Pos for all positions in list subPos
! Returns NsubPos by NPos matrix of true and false
! NOT EFFICIENT IF subPos = Pos (can do much better in that case)
! but algorithm here works well if subPos << Pos
! If have subPos = Pos instead use allNearNeighbors
subroutine nearNeighbors(subPos, Pos, BoxL, lowCut, highCut, NNeighbors, NsubPos, NPos)
    implicit none
    real(8), dimension(NsubPos, 3), intent(in) :: subPos
    real(8), dimension(NPos, 3), intent(in) :: Pos
    real(8), dimension(3), intent(in) :: BoxL
    real(8), intent(in) :: lowCut, highCut
    integer, intent(in) :: NsubPos
    integer, intent(in) :: NPos
    logical, dimension(NsubPos, NPos), intent(out) :: NNeighbors
    integer :: i, j
    real(8), dimension(3) :: iBoxL
    real(8), dimension(3) :: pos1, pos2, distvec
    real(8) :: lowsq, highsq, distvecSq
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0) 
    lowsq = lowCut*lowCut
    highsq = highCut*highCut
    NNeighbors = .false.
    ! Loop over all positions in subPos (may be subset of Pos)
    do i = 1, NsubPos
        pos1 = subPos(i,:)
        ! Loop over all positions in Pos
        do j = 1, NPos
            pos2 = Pos(j,:)
            distvec = pos2 - pos1
            ! MUST re-image individually to account for periodicity to be accurate
            distvec  = distvec - BoxL * anint(distvec * iBoxL)
            distvecSq = sum(distvec**2)
            if ( (distvecSq > lowsq) .and. (distvecSq <= highsq) ) then
                ! lower cut-off not included... lazy, really, but prevents including same atom
                NNeighbors(i,j) = .true.
            endif
        enddo
    enddo
end subroutine

! Same as neighbor searching above, but much more efficient if don't have 
! a subPos, i.e. want all neighbors for all coordinates in Pos
subroutine allNearNeighbors(Pos, BoxL, lowCut, highCut, NNeighbors, NPos)
    implicit none
    real(8), dimension(NPos, 3), intent(in) :: Pos
    real(8), dimension(3), intent(in) :: BoxL
    real(8), intent(in) :: lowCut, highCut
    integer, intent(in) :: NPos
    logical, dimension(NPos, NPos), intent(out) :: NNeighbors
    integer :: i, j
    real(8), dimension(3) :: iBoxL
    real(8), dimension(3) :: pos1, pos2, distvec
    real(8) :: lowsq, highsq, distvecSq
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0) 
    lowsq = lowCut*lowCut
    highsq = highCut*highCut
    NNeighbors = .false.
    ! Loop over all positions in Pos
    do i = 1, NPos
        pos1 = Pos(i,:)
        ! Loop over all other positions
        do j = i+1, NPos
            pos2 = Pos(j,:)
            distvec = pos2 - pos1
            ! MUST re-image individually to account for periodicity to be accurate
            distvec  = distvec - BoxL * anint(distvec * iBoxL)
            distvecSq = sum(distvec**2)
            if ( (distvecSq > lowsq) .and. (distvecSq <= highsq) ) then
                ! lower cut-off not included... lazy, really, but prevents including same atom
                NNeighbors(i,j) = .true.
                NNeighbors(j,i) = .true.
            endif
        enddo
    enddo
end subroutine

! Assesses tetrahedrality of specified set 
! Returns all three-body nearest neighbor angles for all combinations of nearest neighbors
! Note that return is just NneighPos x NneighPos symmetric array
subroutine tetraCosAng(refPos, neighPos, BoxL, allAngs, NneighPos)
    implicit none
    real(8), dimension(3), intent(in) :: refPos
    real(8), dimension(NneighPos, 3), intent(in) :: neighPos
    real(8), dimension(3), intent(in) :: BoxL
    integer, intent(in) :: NneighPos
    real(8), dimension(NneighPos,NneighPos), intent(out) :: allAngs
    integer :: i, j
    real(8) :: tempAng
    real(8), dimension(3) :: distvec1, distvec2, iBoxL
    real*8, external :: CosAngle3
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    ! Loop over nearest neighbor positions
    do i = 1, NneighPos
        distvec1 = neighPos(i,:) - refPos
        distvec1 = distvec1 - BoxL * anint(distvec1 * iBoxL)
        distvec1 = refPos + distvec1
        ! Loop over other neighbor positions (want combinations)
        do j = i+1, NneighPos
            distvec2 = neighPos(j,:) - refPos
            distvec2 = distvec2 - BoxL * anint(distvec2 * iBoxL)
            distvec2 = refPos + distvec2
            ! Find cosine of the angle - MAKE SURE REFPOS IN CENTER
            tempAng = CosAngle3(distvec1, refPos, distvec2)
            allAngs(i,j) = tempAng
            allAngs(j,i) = tempAng
        enddo
    enddo
end subroutine
