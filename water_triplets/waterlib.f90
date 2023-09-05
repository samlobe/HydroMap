! Fortran library for investigating structure and properties of liquid water
! Works for both bulk water and nearby solutes
! SpherePoints and SphereSurfaceAreas are copied from sim.geom.py
! with minor modifications made.


! Same as in geom.py
! Finds centroid of given atom set
subroutine Centroid(Pos, Ret, N, Dim)
    implicit none
    integer, intent(in) :: N, Dim
    real(8), dimension(N,Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(out) :: Ret
    Ret = sum(Pos, 1) / real(N)
end subroutine

! Taken from geom.py
subroutine crossProd3(r1, r2, Ret, Dim)
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: r1, r2
    real(8), dimension(Dim), intent(out) :: Ret
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    Ret(1) = r1(2)*r2(3) - r1(3)*r2(2)
    Ret(2) = r1(3)*r2(1) - r1(1)*r2(3)
    Ret(3) = r1(1)*r2(2) - r1(2)*r2(1)
end subroutine

! Same as in geom.py
subroutine reimage(Pos, RefPos, BoxL, ReimagedPos, NPos, Dim)
    implicit none
    integer, intent(in) :: NPos, Dim
    real(8), dimension(NPos, Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(in) :: RefPos
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NPos, Dim), intent(out) :: ReimagedPos
    integer :: i
    real(8), dimension(Dim) :: distvec, iBoxL
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    do i = 1, NPos
        distvec = Pos(i,:) - RefPos
        distvec = distvec - BoxL * anint(distvec * iBoxL)
        ReimagedPos(i,:) = RefPos + distvec
    enddo            
end subroutine

! Same as in proteinlib.f90
real(8) function RgWeights(Pos, Weights, NAtom, Dim)
  implicit none
  integer, intent(in) :: NAtom, Dim
  real(8), dimension(NAtom,Dim), intent(in) :: Pos
  real(8), dimension(NAtom), intent(in) :: Weights
  real(8), dimension(Dim) :: Center
  integer :: i
  Center = sum(Pos, 1) / real(NAtom)
  RgWeights = 0.
  do i = 1, NAtom
    RgWeights = RgWeights + Weights(i) * sum((Pos(i,:) - Center)**2)
  enddo
  RgWeights = RgWeights / sum(Weights)
  RgWeights = sqrt(RgWeights)
end function

! Same as in geom.py
! Places points on a sphere to later define SASA
subroutine SpherePoints(N, Points)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N,3), intent(out) :: Points
    real(8) :: off, y, phi, r
    real(8) :: inc 
    real(8), parameter :: pi = 3.1415926535897931D0 
    integer :: k
    inc = pi * (3. - sqrt(5.))
    Points = 0.
    off = 2. / real(N)
    do k = 1, N
        y = real(k-1) * off - 1. + (off * 0.5)
        r = sqrt(max(1. - y*y, 0.))
        phi = real(k-1) * inc
        Points(k,1) = cos(phi)*r
        Points(k,2) = y
        Points(k,3) = sin(phi)*r
    enddo
end subroutine

! Modified to also return Exposed, so can find which points (then atoms) are exposed
subroutine SphereSurfaceAreas(Pos, Radii, Points, nExp, BoxL, Areas, Exposed, NSphere, NPoints, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(NSphere, Dim), intent(in) :: Pos
    real(8), dimension(NSphere), intent(in) :: Radii
    real(8), dimension(NPoints, Dim), intent(in) :: Points
    integer, intent(in) :: nExp
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NSphere), intent(out) :: Areas
    logical, dimension(NSphere), intent(out) :: Exposed
    real(8), parameter :: pi = 3.141592653589D0 
    integer, intent(in) :: NSphere, NPoints
    integer :: i, j, k
    real(8), dimension(NPoints,Dim) :: ThisPoints
    real(8) :: AreaPerPoint
    real(8), dimension(NSphere) :: RadiiSq
    real(8), dimension(Dim) :: iPos, jPos
    real(8), dimension(Dim) :: distvec, iBoxL
    logical, dimension(NPoints) :: tempExposed
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0) 
    Areas = 0.
    Exposed = .false.
    RadiiSq = Radii*Radii
    do i = 1, NSphere
        iPos = Pos(i,:)
        AreaPerPoint = 4.*pi*Radii(i)**2 / real(NPoints)
        tempExposed = .true.
        do k = 1, NPoints
            ThisPoints(k,:) = Points(k,:) * Radii(i) + iPos
        enddo
        do j = 1, NSphere
            if (i == j) cycle
            jPos = Pos(j,:)
            distvec = jPos - iPos 
            distvec = distvec - BoxL * anint(distvec * iBoxL)
            jPos = iPos + distvec
            if (.not. any(tempExposed)) exit
            !first check if spheres are far from each other
            if (sum((jPos-iPos)**2) > (Radii(i) + Radii(j))**2) cycle
            do k = 1, NPoints
                if (.not. tempExposed(k)) cycle
                if (sum((ThisPoints(k,:) - jPos)**2) < RadiiSq(j)) tempExposed(k) = .false.
            enddo
        enddo
        Areas(i) = AreaPerPoint * real(count(tempExposed))
        if (count(tempExposed) >= nExp) Exposed(i) = .true.
    enddo
end subroutine

! Same as in geom.py
subroutine SphereVolumes(Pos, Radii, dx, Volumes, NSphere, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(NSphere, Dim), intent(in) :: Pos
    real(8), dimension(NSphere), intent(in) :: Radii
    real(8),  intent(in) :: dx
    real(8), dimension(NSphere), intent(out) :: Volumes
    integer, intent(in) :: NSphere
    real(8), dimension(NSphere) :: RadiiSq
    real(8) :: minDistSq, DistSq, dV
    integer :: i,j
    real(8), dimension(Dim) :: Pos2, minPos, maxPos
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    RadiiSq = Radii*Radii
    Volumes = 0.
    dV = dx*dx*dx
    minPos = (/minval(Pos(:,1) - Radii), minval(Pos(:,2) - Radii), minval(Pos(:,3) - Radii)/)
    maxPos = (/maxval(Pos(:,1) + Radii), maxval(Pos(:,2) + Radii), maxval(Pos(:,3) + Radii)/)
    maxPos = maxPos + dx * 0.5    
    !first do a coarse grid check to see which spheres are where
    Pos2 = minPos
    do while (all(Pos2 < maxPos))
        j = 0
        minDistSq = huge(1.d0)
        do i = 1, NSphere
            DistSq = sum((Pos(i,:) - Pos2)**2)
            if (DistSq < minDistSq .and. DistSq < RadiiSq(i)) then
                minDistSq = DistSq
                j = i
            endif
        enddo
        if (j > 0) Volumes(j) = Volumes(j) + dV
        Pos2(1) = Pos2(1) + dx
        do i = 1, 2
            if (Pos2(i) >= maxPos(i)) then
                Pos2(i) = minPos(i)
                Pos2(i+1) = Pos2(i+1) + dx
            endif
        enddo
    enddo   
end subroutine

! Calculates average g(r) for atoms in Pos2 to atoms in Pos1
! Will work best for protein if use SURFACE atoms as Pos1, not all
! Assumes a bulk density of atoms in Pos2 provided... i.e. system can be inhomogeneous
! rather than assuming homogeneously distributed atoms in Pos2 over whole volume of box
subroutine RadialDist(Pos1, Pos2, binwidth, totbins, BulkDens, BoxL, rdf, NPos1, NPos2)
    implicit none
    real(8), dimension(NPos1, 3), intent(in) :: Pos1
    real(8), dimension(NPos2, 3), intent(in) :: Pos2
    real(8), intent(in) :: binwidth
    integer, intent(in) :: totbins
    ! Total number of bins, so max bin value is totbins*binwidth; if larger ignored
    real(8), intent(in) :: BulkDens
    real(8), dimension(3), intent(in) :: BoxL
    real(8), dimension(totbins), intent(out) :: rdf
    integer, intent(in) :: NPos1, NPos2
    real(8), parameter :: pi = 3.141592653589D0
    real(8), dimension(totbins) :: counts
    real(8), dimension(3) :: iPos, jPos, iBoxL, distVec
    real(8) :: dist
    integer :: nbin, i, j, k
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    counts = 0.
    do i = 1, NPos2
        iPos = Pos2(i,:)
        do j = 1, NPos1
            jPos = Pos1(j,:)
            distVec = jPos - iPos
            distVec = distVec - BoxL * anint(distVec * iBoxL)
            dist = sqrt(sum(distVec**2))
            ! Place in correct bin
            nbin = ceiling(dist / binwidth) ! Bin i has all distances less than i but greater than i - 1
            if (nbin <= totbins) then
                counts(nbin) = counts(nbin) + 1.
            endif
            ! If beyond max of totbins*binwidth, just don't include
        enddo
    enddo
    ! Normalize counts correctly
    do k = 1, totbins
        rdf(k) = counts(k) / (NPos1 * BulkDens * (4./3.) * pi * (binwidth**3) * ((k**3) - ((k-1)**3)))
    enddo
    ! Gives g(r) for single snapshot... for whole trajectory, do for each frame and average
end subroutine

subroutine RadialDistSame(Pos, binwidth, totbins, BulkDens, BoxL, rdf, NPos)
    implicit none
    real(8), dimension(NPos, 3), intent(in) :: Pos
    real(8), intent(in) :: binwidth
    integer, intent(in) :: totbins
    ! Total number of bins, so max bin value is totbins*binwidth; if larger ignored
    real(8), intent(in) :: BulkDens
    real(8), dimension(3), intent(in) :: BoxL
    real(8), dimension(totbins), intent(out) :: rdf
    integer, intent(in) :: NPos
    real(8), parameter :: pi = 3.141592653589D0
    real(8), dimension(totbins) :: counts
    real(8), dimension(3) :: iPos, jPos, iBoxL, distVec
    real(8) :: dist
    integer :: nbin, i, j, k
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    counts = 0.
    do i = 1, NPos
        iPos = Pos(i,:)
        do j = i+1, NPos
            jPos = Pos(j,:)
            distVec = jPos - iPos
            distVec = distVec - BoxL * anint(distVec * iBoxL)
            dist = sqrt(sum(distVec**2))
            ! Place in correct bin
            nbin = ceiling(dist / binwidth) ! Bin i has all distances less than i but greater than i - 1
            if (nbin <= totbins) then
                counts(nbin) = counts(nbin) + 1.
            endif
            ! If beyond max of totbins*binwidth, just don't include
        enddo
    enddo
    ! Normalize counts correctly
    do k = 1, totbins
        rdf(k) = counts(k) / (NPos * BulkDens * (4./3.) * pi * (binwidth**3) * ((k**3) - ((k-1)**3)))
    enddo
    ! Gives g(r) for single snapshot... for whole trajectory, do for each frame and average
end subroutine

! Write more general function for finding probability distribution of distances in any dimension, with reimaging
! If Pos1 and Pos2 are the same, then divide the number of counts by 2
! If the distance is zero, it's not counted
subroutine PairDistanceHistogram(Pos1, Pos2, binwidth, totbins, BoxL, hist, NPos1, NPos2, Dim)
    implicit none
    integer, intent(in) :: NPos1, NPos2, Dim
    real(8), dimension(NPos1, Dim), intent(in) :: Pos1
    real(8), dimension(NPos2, Dim), intent(in) :: Pos2
    real(8), intent(in) :: binwidth
    integer, intent(in) :: totbins
    ! Total number of bins, so max bin value is totbins*binwidth; if larger ignored
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(totbins), intent(out) :: hist
    real(8), dimension(Dim) :: iPos, jPos, iBoxL, distVec
    real(8) :: dist
    integer :: nbin, i, j
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    hist = 0.
    do i = 1, NPos1
        iPos = Pos1(i,:)
        do j = 1, NPos2
            jPos = Pos2(j,:)
            distVec = jPos - iPos
            distVec = distVec - BoxL * anint(distVec * iBoxL)
            dist = sqrt(sum(distVec**2))
            if (dist == 0.d0) cycle ! Skip it if it's the same position
            ! Place in correct bin
            nbin = ceiling(dist / binwidth) ! Bin i has all distances less than i but greater than i - 1
            if (nbin <= totbins) then
                hist(nbin) = hist(nbin) + 1.
            endif
            ! If beyond max of totbins*binwidth, just don't include
        enddo
    enddo
end subroutine

! Finds which atoms in first and second solvation shells
! Should give EXTERIOR/SURFACE atoms (best efficiency) to Pos1 and solvent atoms to Pos2
subroutine SolventShells(Pos1, Pos2, dist1, dist2, shell1, shell2, NPos1, NPos2)
    implicit none
    real(8), dimension(NPos1, 3), intent(in) :: Pos1
    real(8), dimension(NPos2, 3), intent(in) :: Pos2
    real(8), intent(in) :: dist1, dist2
    logical, dimension(NPos2), intent(out) :: shell1, shell2
    integer, intent(in) :: NPos1, NPos2
    real(8), dimension(NPos1) :: tempdists
    real(8), dimension(3) :: iPos, jPos
    integer :: i, j
    real(8) :: dist1Sq, dist2Sq
    shell1 = .false.
    shell2 = .false.
    dist1Sq = dist1*dist1
    dist2Sq = dist2*dist2
    do i = 1, NPos2
        iPos = Pos2(i,:)
        tempdists = 0.
        do j = 1, NPos1
            jPos = Pos1(j,:)
            tempdists(j) = sum((jPos - iPos)**2)
        enddo
        if (minval(tempdists) <= dist1Sq) then
            shell1(i) = .true.
        else if (minval(tempdists) <= dist2Sq) then
            shell2(i) = .true.
        endif
    enddo
end subroutine

! Finds H-bonds SPECIFICALLY between a peptide and nearby waters, positions provided
! Weird implementation - also, H-bonds found geometrically not that useful  
! Used to analyze simulations of peptide pulled of CH3 SAM, so keeping
! To really analyze H-bonds, should really use updated version, generalHbonds
subroutine FindHbonds(pepAcc, pepDon, watPos, distCut, angCut, &
                      NBonds, watAccOut, watDonOut, pepAccOut, pepDonOut, NpepAcc, NpepDon, NwatPos)
    implicit none
    real(8), dimension(NpepAcc, 3), intent(in) :: pepAcc
    real(8), dimension(NpepDon, 3), intent(in) :: pepDon
    real(8), dimension(NwatPos, 3), intent(in) :: watPos
    real(8), intent(in) :: distCut, angCut
    integer, intent(out) :: NBonds
    integer, dimension(NwatPos), intent(out) :: watAccOut, watDonOut 
    integer, dimension(NpepAcc), intent(out) :: pepAccOut
    integer, dimension(NpepDon), intent(out) :: pepDonOut
    integer, intent(in) :: NpepAcc, NpepDon, NwatPos
    real(8) :: distCutSq, cosAngCut
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), parameter :: RadPerDeg = pi / 180.D0
    integer j, i
    real(8), dimension(3) :: oxPos, h1Pos, h2Pos, vecWat1, vecWat2
    real(8), dimension(3) :: apos, bondVec1, bondVec2, heavyPos, HPos, vecPep
    real(8) :: bondDist1Sq, bondDist2Sq, cosAng1, cosAng2
    NBonds = 0
    watDonOut = 0
    watAccOut = 0
    pepAccOut = 0
    pepDonOut = 0
    distCutSq = distCut * distCut
    cosAngCut = cos(angCut * RadPerDeg)
    ! Check length of water position array
    if (mod(NwatPos, 3) /= 0) then
        print *, "Water position array does not have length of a factor of 3."
        stop
    endif
    if (mod(NpepDon, 2) /= 0) then
        print *, "Peptide donor position array does not have length of a factor of 2."
        stop
    endif
    ! Loop over all water molecules
    do j = 1, NwatPos, 3
        oxPos = watPos(j,:)
        h1Pos = watPos(j+1,:)
        h2Pos = watPos(j+2,:)
        vecWat1 = h1Pos - oxPos
        vecWat2 = h2Pos - oxPos
        ! Loop over peptide acceptors
        do i = 1, NpepAcc
            apos = pepAcc(i,:)
            bondVec1 = apos - h1Pos
            bondDist1Sq = sum(bondVec1**2)
            if (bondDist1Sq < distCutSq) then
                cosAng1 = dot_product(bondVec1, vecWat1) / sqrt(bondDist1Sq * sum(vecWat1**2))
                if (cosAng1 > cosAngCut) then
                    NBonds = NBonds + 1
                    pepAccOut(i) = pepAccOut(i) + 1
                    watDonOut(j+1) = watDonOut(j+1) + 1
                    cycle
                    ! If already have H-bond, don't check for second!
                endif
            endif
            bondVec2 = apos - h2Pos
            bondDist2Sq = sum(bondVec2**2)
            if (bondDist2Sq < distCutSq) then
                cosAng2 = dot_product(bondVec2, vecWat2) / sqrt(bondDist2Sq * sum(vecWat2**2))
                if (cosAng2 > cosAngCut) then
                    NBonds = NBonds + 1
                    pepAccOut(i) = pepAccOut(i) + 1
                    watDonOut(j+2) = watDonOut(j+2) + 1
                endif
            endif
        enddo
        ! Loop over peptide donors (hydrogens with attached heavy atom 1 index behind!)
        do i = 1, NpepDon, 2
            heavyPos = pepDon(i,:)
            HPos = pepDon(i+1,:)
            bondVec1 = oxPos - HPos
            bondDist1Sq = sum(bondVec1**2)
            if (bondDist1Sq < distCutSq) then
                vecPep = HPos - heavyPos
                cosAng1 = dot_product(bondVec1, vecPep) / sqrt(bondDist1Sq * sum(vecPep**2))
                if (cosAng1 > cosAngCut) then
                    NBonds = NBonds + 1
                    pepDonOut(i+1) = pepDonOut(i+1) + 1
                    watAccOut(j) = watAccOut(j) + 1
                endif
            endif
        enddo
    enddo
end subroutine

! Finds H-Bonds for peptide backbone only 
! Used to analyze simulations of peptide pulled of CH3 SAM, so keeping
! To really analyze H-bonds, should really use updated version, generalHbonds
subroutine BBHbonds(pepAcc, pepDon, distCut, angCut, NBonds, pepAccOut, pepDonOut, NpepAcc, NpepDon)
    implicit none
    real(8), dimension(NpepAcc, 3), intent(in) :: pepAcc
    real(8), dimension(NpepDon, 3), intent(in) :: pepDon
    real(8), intent(in) :: distCut, angCut
    integer, intent(out) :: NBonds
    integer, dimension(NpepAcc), intent(out) :: pepAccOut
    integer, dimension(NpepDon), intent(out) :: pepDonOut
    integer, intent(in) :: NpepAcc, NpepDon
    real(8) :: distCutSq, cosAngCut
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), parameter :: RadPerDeg = pi / 180.D0
    integer j, i
    real(8), dimension(3) :: apos, bondVec1, heavyPos, HPos, vecPep
    real(8) :: bondDist1Sq, cosAng1
    NBonds = 0
    pepAccOut = 0
    pepDonOut = 0
    distCutSq = distCut * distCut
    cosAngCut = cos(angCut * RadPerDeg)
    if (mod(NpepDon, 2) /= 0) then
        print *, "Peptide donor position array does not have length of a factor of 2."
        stop
    endif
    ! Loop over peptide donors (hydrogens with attached heavy atom 1 index behind!)
    do i = 1, NpepDon, 2
        heavyPos = pepDon(i,:)
        HPos = pepDon(i+1,:)
        ! Loop over peptide acceptors
        ! Could be more efficient by excluding H attached to same heavy atom...
        ! Simpler this way, and ok because will register cosAng1 of -1 and not be considered
        do j = 1, NpepAcc
            apos = pepAcc(j,:)
            bondVec1 = apos - HPos
            bondDist1Sq = sum(bondVec1**2)
            if (bondDist1Sq < distCutSq) then
                vecPep = HPos - heavyPos
                cosAng1 = dot_product(bondVec1, vecPep) / sqrt(bondDist1Sq * sum(vecPep**2))
                if (cosAng1 > cosAngCut) then
                    NBonds = NBonds + 1
                    pepDonOut(i+1) = pepDonOut(i+1) + 1
                    pepAccOut(j) = pepAccOut(j) + 1
                endif
            endif
        enddo
    enddo
end subroutine

! Finds H-Bonds between waters in a group, given full water coordinates
! Need to re-image EACH computed distance for accuracy!
! Not like above where only want H-bonds close to peptide, which is re-image reference in python code
! Used to analyze simulations of peptide pulled of CH3 SAM, so keeping
! To really analyze H-bonds, should really use updated version, generalHbonds
subroutine WatHbonds(watPos, allWatPos, BoxL, distCut, angCut, NBonds, watAccOut, watDonOut, NwatPos, NallWatPos)
    implicit none
    real(8), dimension(NwatPos, 3), intent(in) :: watPos
    real(8), dimension(NallWatPos, 3), intent(in) :: allWatPos
    real(8), dimension(3), intent(in) :: BoxL
    real(8), intent(in) :: distCut, angCut
    integer, intent(out) :: NBonds
    integer, dimension(NwatPos), intent(out) :: watAccOut
    integer, dimension(NwatPos), intent(out) :: watDonOut
    integer, intent(in) :: NwatPos
    integer, intent(in) :: NallWatPos
    real(8) :: distCutSq, cosAngCut
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), parameter :: RadPerDeg = pi / 180.D0
    integer j, i
    real(8), dimension(3) :: bondVec1, bondVec2, heavyPos, H1Pos, H2Pos, avec1, avec2
    real(8), dimension(3) :: oxpos, watH1pos, watH2pos
    real(8) :: bondDist1Sq, bondDist2Sq, cosAng1, cosAng2
    real(8), dimension(3) :: iBoxL
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    NBonds = 0
    watAccOut = 0
    watDonOut = 0
    distCutSq = distCut * distCut
    cosAngCut = cos(angCut * RadPerDeg)
    if (mod(NwatPos, 3) /= 0) then
        print *, "Waters don't all have 3 atoms!"
        stop
    endif
    ! Loop over water molecules in set of interest, indexing off H1, with O at i-1, H2 at i+1
    do i = 2, NwatPos, 3
        heavyPos = watPos(i-1,:)
        H1Pos = watPos(i,:)
        H2Pos = watPos(i+1,:)
        ! Loop over ALL water molecules (including subset)
        ! Less efficient, but ok because will give cosAng1 of -1 if donor and acceptor in same molecule
        do j = 1, NallWatPos, 3
            oxpos = allWatPos(j,:)
            ! First look for H-bonds involving in-set hydrogens (donors) with all oxygens (acceptors)
            bondVec1 = oxpos - H1Pos 
            bondVec1 = bondVec1 - BoxL * anint(bondVec1 * iBoxL)
            bondDist1Sq = sum(bondVec1**2)
            if (bondDist1Sq < distCutSq) then
                avec1 = H1Pos - heavyPos
                ! Don't re-image because water molecules kept together (by AMBER)
                cosAng1 = dot_product(bondVec1, avec1) / sqrt(bondDist1Sq * sum(avec1**2))
                if (cosAng1 > cosAngCut) then
                    NBonds = NBonds + 1
                    watDonOut(i) = watDonOut(i) + 1
                    if (j <= NwatPos) then
                        watAccOut(j) = watAccOut(j) + 1
                    endif
                    ! Don't bother checking other hydrogen if already one bonded to this O
                    cycle
                endif
            endif
            bondVec2 = oxpos - H2Pos
            bondVec2 = bondVec2 - BoxL * anint(bondVec2 * iBoxL)
            bondDist2Sq = sum(bondVec2**2)
            if (bondDist2Sq < distCutSq) then
                avec2 = H2Pos - heavyPos
                cosAng2 = dot_product(bondVec2, avec2) / sqrt(bondDist2Sq * sum(avec2**2))
                if (cosAng2 > cosAngCut) then
                    NBonds = NBonds + 1
                    watDonOut(i+1) = watDonOut(i+1) + 1
                    if (j <= NwatPos) then
                        watAccOut(j) = watAccOut(j) + 1
                    endif
                endif
            endif
        enddo
        do j = 2, NallWatPos, 3
            oxpos = allWatPos(j-1,:)
            watH1pos = allWatPos(j,:)
            watH2pos = allWatPos(j+1,:)
            ! And repeat, but for H-bonds involving in-set oxygens and out of set hydrogens
            bondVec1 = heavyPos - watH1pos
            bondVec1 = bondVec1 - BoxL * anint(bondVec1 * iBoxL)
            bondDist1Sq = sum(bondVec1**2)
            if (bondDist1Sq < distCutSq) then
                avec1 = watH1pos - oxpos
                ! Don't re-image because water molecules kept together (by AMBER)
                cosAng1 = dot_product(bondVec1, avec1) / sqrt(bondDist1Sq * sum(avec1**2))
                if (cosAng1 > cosAngCut) then
                    NBonds = NBonds + 1
                    if (j <= NwatPos) then
                        watDonOut(j) = watDonOut(j) + 1
                    endif
                    watAccOut(i-1) = watAccOut(i-1) + 1
                    ! Don't bother checking other hydrogen if already one bonded to this O
                    cycle
                endif
            endif
            bondVec2 = heavyPos - watH2pos
            bondVec2 = bondVec2 - BoxL * anint(bondVec2 * iBoxL)
            bondDist2Sq = sum(bondVec2**2)
            if (bondDist2Sq < distCutSq) then
                avec2 = watH2pos - oxpos
                cosAng2 = dot_product(bondVec2, avec2) / sqrt(bondDist2Sq * sum(avec2**2))
                if (cosAng2 > cosAngCut) then
                    NBonds = NBonds + 1
                    if (j <= NwatPos) then
                        watDonOut(j+1) = watDonOut(j+1) + 1
                    endif
                    watAccOut(i-1) = watAccOut(i-1) + 1
                endif
            endif
        enddo
    enddo
end subroutine

! Modified from geometrylib.py to return angle instead of cosine of angle
! So called CosAngle3, but really returns angle given three coordinate positions
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

! Calculates squared displacements of list of atoms from themselves, applying unwrapping 
! as necessary - hinges on fact that ANINT returns nearest integer and that position 
! cannot physically change by more than half box length between pos and prevPos
subroutine calcSD(pos, prevPos, refPos, BoxL, allSD, newPos, Npos, NprevPos, NrefPos)
    implicit none
    real(8), dimension(Npos, 3), intent(in) :: pos
    real(8), dimension(NprevPos, 3), intent(in) :: prevPos
    real(8), dimension(NrefPos, 3), intent(in) :: refPos
    real(8), dimension(3), intent(in) :: BoxL
    integer, intent(in) :: Npos, NprevPos, NrefPos
    real(8), dimension(Npos, 3), intent(out) :: allSD
    real(8), dimension(Npos, 3), intent(out) :: newPos
    real(8), dimension(3) :: iBoxL, avec, distvec
    integer :: i
    ! First check if Npos matches NrefPos, if not quit
    if (Npos /= NrefPos .or. NprevPos /= Npos) then
        print *, 'Number of current and reference coordinates must match.', Npos, NprevPos, NrefPos
        stop
    endif
    ! Set up imaging/unwrapping
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    ! Loop over coordinates, calc squared displacement
    do i = 1, Npos
        ! First unwrap to find new coordinate without PBCs
        avec = pos(i,:) - prevPos(i,:)
        avec = avec - BoxL * anint(avec * iBoxL)
        newPos(i,:) = prevPos(i,:) + avec
        ! Now calculate squared displacement from reference
        distvec = newPos(i,:) - refPos(i,:)
        allSD(i,:) = distvec**2
    enddo
end subroutine 

! Define function to just return angle given two normalized vectors
real*8 function AngBetween(Vec1, Vec2)
  implicit none
  real*8, parameter :: pi = 3.1415926535897931D0
  real*8, parameter :: DegPerRad = 180.D0/pi
  real*8, dimension(3), intent(in) :: Vec1, Vec2
  real*8 :: Phi
  Phi = min(1.D0, max(-1.D0, dot_product(Vec1, Vec2)))
  Phi = acos(Phi)
  AngBetween = mod(Phi + pi, pi*2.D0) - pi
  if (AngBetween < -pi) AngBetween = AngBetween + pi*2.D0
  AngBetween = AngBetween * DegPerRad
end function

! Given oxygen and hydrogen positions (arrays must be in same order, so first two H
! correspond to the first O) and a reference vector, computes the angle between
! the reference and the water dipoles, and the vector normal to the plane of the
! oxygen and hydrogen atoms in the water molecules. This is much more informative
! than reporting both the dipole and both the OH bond vectors, I think.
! Note that the water dipole direction is just the sum of the OH vectors.
subroutine watOrient(opos, hpos, refvec, BoxL, angDip, angPlane, Nopos, Nhpos)
    implicit none
    real(8), dimension(Nopos, 3), intent(in) :: opos
    real(8), dimension(Nhpos, 3), intent(in) :: hpos
    real(8), dimension(3), intent(in) :: refvec
    real(8), dimension(3), intent(in) :: BoxL
    integer, intent(in) :: Nopos, Nhpos
    real(8), dimension(Nopos), intent(out) :: angDip
    real(8), dimension(Nopos), intent(out) :: angPlane
    external :: crossProd3
    real*8, external :: AngBetween
    real(8), dimension(3) :: vecoh1, vecoh2, vecdip, vecplane, refvecnorm, iBoxL
    real(8) :: normdip, normplane
    integer :: i
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    ! Check to make sure number of hydrogens matches oxygens
    if (2*Nopos /= Nhpos) then
        print *, 'Number of hydrogens must be two times number of oxygens.'
        stop
    endif
    ! Make sure refvec normalized
    refvecnorm = refvec/sqrt(sum(refvec*refvec))
    ! Now loop over oxygens
    do i = 1, Nopos
        vecoh1 = hpos(2*i - 1,:) - opos(i,:)
        ! Image just in case
        vecoh1 = vecoh1 - BoxL * anint(vecoh1 * iBoxL)
        vecoh2 = hpos(2*i,:) - opos(i,:)
        vecoh2 = vecoh2 - BoxL * anint(vecoh2 * iBoxL)
        vecdip = vecoh1 + vecoh2
        vecdip = vecdip - BoxL * anint(vecdip * iBoxL)
        normdip = sqrt(sum(vecdip*vecdip))
        angDip(i) = AngBetween(vecdip/normdip, refvecnorm)
        call crossProd3(vecoh1, vecoh2, vecplane, 3)
        normplane = sqrt(sum(vecplane*vecplane))
        angPlane(i) = AngBetween(vecplane/normplane, refvecnorm)
    enddo
end subroutine

! Given oxygen and hydrogen positions (arrays must be in same order, so first two H
! correspond to the first O) and a reference vector, computes the angle between
! the reference and the water dipoles, and the vector normal to the plane of the
! oxygen and hydrogen atoms in the water molecules. This is much more informative
! than reporting both the dipole and both the OH bond vectors, I think.
! Note that the water dipole direction is just the sum of the OH vectors.
subroutine watOHvec(opos, hpos, BoxL, normdip, Nopos, Nhpos)
    implicit none
    real(8), dimension(Nopos, 3), intent(in) :: opos
    real(8), dimension(Nhpos, 3), intent(in) :: hpos
    real(8), dimension(3), intent(in) :: BoxL
    integer, intent(in) :: Nopos, Nhpos
    real(8), dimension(Nopos,3), intent(out) :: normdip
    real(8), dimension(3) :: vecoh1, vecoh2, vecdip, iBoxL
    integer :: i
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    ! Check to make sure number of hydrogens matches oxygens
    if (2*Nopos /= Nhpos) then
        print *, 'Number of hydrogens must be two times number of oxygens.'
        stop
    endif
    ! Now loop over oxygens
    do i = 1, Nopos
        vecoh1 = hpos(2*i - 1,:) - opos(i,:)
        ! Image just in case
        vecoh1 = vecoh1 - BoxL * anint(vecoh1 * iBoxL)
        vecoh2 = hpos(2*i,:) - opos(i,:)
        vecoh2 = vecoh2 - BoxL * anint(vecoh2 * iBoxL)
        vecdip = vecoh1 + vecoh2
        vecdip = vecdip - BoxL * anint(vecdip * iBoxL)
        normdip(i,:) = vecdip/sqrt(sum(vecdip*vecdip))
    enddo
end subroutine

subroutine watHHvec(h1pos, h2pos, BoxL, normhh, Nh1pos, Nh2pos)
    implicit none
    real(8), dimension(Nh1pos, 3), intent(in) :: h1pos
    real(8), dimension(Nh2pos, 3), intent(in) :: h2pos
    real(8), dimension(3), intent(in) :: BoxL
    integer, intent(in) :: Nh1pos, Nh2pos
    real(8), dimension(Nh1pos,3), intent(out) :: normhh
    real(8), dimension(3) :: vechh, iBoxL
    integer :: i
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    ! Check to make sure number of hydrogens matches oxygens
    if (Nh1pos /= Nh2pos) then
        print *, 'Number of hydrogens must be the same.'
        stop
    endif
    ! Now loop over oxygens
    do i = 1, Nh1pos
        vechh = h1pos(i,:) - h2pos(i,:)
        ! Image just in case
        vechh = vechh - BoxL * anint(vechh * iBoxL)
        normhh(i,:) = vechh/sqrt(sum(vechh*vechh))
    enddo
end subroutine

! Bins x,y,z coordinates of oxygens based on provided x, y, z bins
subroutine binOnGrid(opos, xbins, ybins, zbins, outhist, Nopos, Nx, Ny, Nz)
    implicit none
    real(8), dimension(Nopos, 3), intent(in) :: opos
    real(8), dimension(Nx), intent(in) :: xbins
    real(8), dimension(Ny), intent(in) :: ybins
    real(8), dimension(Nz), intent(in) :: zbins
    integer, intent(in) :: Nopos, Nx, Ny, Nz
    integer, dimension(Nx-1, Ny-1, Nz-1), intent(out) :: outhist
    real(8) :: binwidth, radsq, thisdistsq
    integer :: thisxbin, thisybin, thiszbin
    real(8), dimension(3) :: thisvec
    integer :: i
    ! Find all bin widths - ASSUMES UNIFORM BINS!
    binwidth = xbins(2) - xbins(1)
    ! Also assumes all binwidths are the same in each dimension, i.e. have cubes
    if (ybins(2) - ybins(1) /= binwidth .or. zbins(2) - zbins(1) /= binwidth) then
        print *, 'Must break volume into CUBES. Currently, bin-widths do not match.'
        stop 
    endif
    radsq = binwidth*binwidth/4.0d0
    ! Initiate outhist
    outhist = 0
    do i = 1, Nopos
        ! Find x, y, and z bins - inclusive on left edge, exclusive on right
        ! Not considered if out of bin range
        thisxbin = floor((opos(i,1) - xbins(1)) / binwidth) + 1
        ! Make sure actually in range of bins
        if (thisxbin < 1 .or. thisxbin > (Nx-1)) then
            cycle
        else
            thisybin = floor((opos(i,2) - ybins(1)) / binwidth) + 1
            if (thisybin < 1 .or. thisybin > (Ny-1)) then
                cycle
            else
                thiszbin = floor((opos(i,3) - zbins(1)) / binwidth) + 1
                if (thiszbin < 1 .or. thiszbin > (Nz-1)) then
                    cycle
                else
                    ! Add one to outhist if also in sphere within cube
                    thisvec(1) = opos(i,1) - (xbins(thisxbin) + binwidth*0.5d0)
                    thisvec(2) = opos(i,2) - (ybins(thisybin) + binwidth*0.5d0)
                    thisvec(3) = opos(i,3) - (zbins(thiszbin) + binwidth*0.5d0)
                    thisdistsq = sum(thisvec*thisvec)
                    if (thisdistsq <= radsq) then
                        outhist(thisxbin, thisybin, thiszbin) = outhist(thisxbin, thisybin, thiszbin) + 1
                    else
                        cycle
                    endif
                endif
            endif
        endif
    enddo
end subroutine

! Also want to be able to use a finely spaced grid, but with arbitrarily sized spherical probe volumes
! with centers at the grid points. So given an array of 3D positions, reimages around
! each grid point and selects all positions within probeRadius. Then returns an array
! with the number of positions from Pos in each grid position in gridPos. Not as efficient
! as binOnGrid, but allows more samples per simulation frame.
subroutine probeGrid(Pos, gridPos, probeRadius, BoxL, numGrid, NPos, NgridPos)
    implicit none
    real(8), dimension(NPos, 3), intent(in) :: Pos
    real(8), dimension(NgridPos, 3), intent(in) :: gridPos
    real(8), intent(in) :: probeRadius
    real(8), dimension(3), intent(in) :: BoxL
    integer, intent(in) :: NPos, NgridPos
    integer, dimension(NgridPos), intent(out) :: numGrid
    real(8), dimension(3) :: iBoxL, refPos, thisPos, distVec
    real(8) :: radiusSq, distSq
    integer :: i, j
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    radiusSq = probeRadius*probeRadius
    ! Loop over grid points
    do i = 1, NgridPos
        refPos = gridPos(i,:)
        ! Loop over other positions
        do j = 1, NPos
            thisPos = Pos(j,:)
            ! Check if squared reimaged distance is below cut-off
            distVec = thisPos - refPos
            distVec = distVec - BoxL * anint(distVec * iBoxL)
            distSq = sum(distVec*distVec)
            if (distSq <= radiusSq) then
                numGrid(i) = numGrid(i) + 1
            endif
        enddo
    enddo
end subroutine

! Below is a very general algorithm for computing H-bonds
! This should be used instead of the other functions for finding H-bonds
! found in this library. Those other functions have been left because they
! were used in the peptide-surface pulling simulation analysis.
! This code will not be efficient in every situation, but I tried to make it
! as general, simple to understand, and reliable as possible.
! acceptorPos - the acceptor heavy-atom positions
! donorPos - the donor heavy-atom positions
! donorHPos - the donor hydrogen positions; note that for something like water
!             with two hydrogens per heavy atom oxygen donor, each oxygen will
!             need to be listed TWICE in donorPos such that donorPos and donorHPos
!             are the same length with entries with indices that match for the same
!             water molecule.
! BoxL - box dimensions for wrapping
! distCut - distance cut-off
! angCut - angle cut-off, in degrees, between acceptor-H vector and donor-H vector
!          In other words, 180 is a straight line in same direction and 0 is opposite directions
! bondBool - output in the form of a boolean array that is of size Nacc x Ndon
!            Entry of 1 or True means that the two members are H-bonded
!            The boolean array is NOT symmetric 
subroutine generalHbonds(acceptorPos, donorPos, donorHPos, BoxL, distCut, angCut, bondBool, Nacc, Ndon, Ndonh)
    implicit none
    real(8), dimension(Nacc, 3), intent(in) :: acceptorPos
    real(8), dimension(Ndon, 3), intent(in) :: donorPos
    real(8), dimension(Ndonh, 3), intent(in) :: donorHPos
    real(8), dimension(3), intent(in) :: BoxL
    real(8), intent(in) :: distCut, angCut
    integer, intent(in) :: Nacc, Ndon, Ndonh
    logical, dimension(Nacc, Ndon), intent(out) :: bondBool
    real*8, external :: AngBetween
    real(8), dimension(3) :: iBoxL, posAcc, posDon, distVec, posH, accVec, donVec
    real(8) :: distSq, distCutSq, ang
    integer :: i, j
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    ! Check to make sure that every donor H has a heavy atom associated with it
    if (Ndonh /= Ndon) then
        print *, 'Number of donor hydrogens and heavy-atoms do not match.'
        stop
    endif
    distCutSq = distCut*distCut
    bondBool = .false.
    ! Loop over acceptors
    do i = 1, Nacc
        posAcc = acceptorPos(i,:)
        ! Loop over donor heavy-atoms
        do j = 1, Ndon
            posDon = donorPos(j,:)
            ! Image the vector between to get distance
            distVec = posDon - posAcc
            distVec = distVec - BoxL * anint(distVec * iBoxL)
            distSq = sum(distVec*distVec)
            if (distSq > distCutSq .OR. distSq <= 1.0E-2) then
                cycle ! Above also excludes if distance is zero... i.e. same atom
            else
                ! Now check angle cut-off, making sure to normalize vectors
                posH = donorHPos(j,:)
                accVec = posAcc - posH
                ! Make sure to image
                accVec = accVec - BoxL * anint(accVec * iBoxL)
                accVec = accVec / sqrt(sum(accVec*accVec))
                donVec = posDon - posH
                donVec = donVec - BoxL * anint(donVec * iBoxL)
                donVec = donVec / sqrt(sum(donVec*donVec))
                ang = AngBetween(accVec, donVec)
                ! To be H-bonded, angle must be ABOVE cut-off for D-H----A angle
                ! Note that straight line is 180 degrees since defined H as origin for two vectors
                if (ang < angCut) then
                    cycle
                else
                    bondBool(i, j) = .true.
                endif
            endif
        enddo 
    enddo
end subroutine


! Would like to create a Willard-Chandler interface from instantaneous solvent positions
! However, this requires an algorithm called "marching cubes" to find isosurfaces
! I have no desire to code this algorithm myself, and it is available from scikit-image,
! both in its classic version, and with improvements
! Below generates the density field (density values and unit normal vectors) for a given
! set of water oxygen locations and x, y, and z grid locations
! If want to generalize, it just creates a Gaussian density field from points, with each
! Gaussian truncated and shifted to zero beyond 3*smoothlen
! pos - water oxygen positions
! gridx - x-dimension grid locations
! gridy - y-dimension grid locations
! gridz - z-dimension grid locations
! BoxL - the box dimensions (for wrapping)
! smoothlen - the smoothing len (sigma)
! densvals - all density values at the grid points
! densnorms - all normal vectors at the grid points
subroutine WillardDensityField(pos, gridx, gridy, gridz, BoxL, smoothlen, densvals, densnorms, Npos, Nx, Ny, Nz)
    implicit none
    real(8), dimension(Npos, 3), intent(in) :: pos
    real(8), dimension(Nx), intent(in) :: gridx
    real(8), dimension(Ny), intent(in) :: gridy
    real(8), dimension(Nz), intent(in) :: gridz
    real(8), dimension(3), intent(in) :: BoxL
    real(8), intent(in) :: smoothlen
    integer, intent(in) :: Npos, Nx, Ny, Nz
    real(8), dimension(Nx, Ny, Nz), intent(out) :: densvals
    real(8), dimension(Nx, Ny, Nz, 3), intent(out) :: densnorms
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), dimension(3) :: iBoxL, apos, watpos, normvec, thisvec, normfunc
    real(8) :: shiftterm, densval, rvalsq, expterm, densfunc
    integer :: i, j, k, l
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    shiftterm = exp(-9.0D0/2.0D0) / ((2.0D0*pi*(smoothlen*smoothlen))**(1.5D0))
    ! Loop over grid in each dimension to get all points
    do i = 1, Nx
        do j = 1, Ny
            do k = 1, Nz
                apos(1) = gridx(i)
                apos(2) = gridy(j)
                apos(3) = gridz(k)
                ! Set density function and gradient of density to zero
                ! These are both functions of all the solvent coordinates, so will sum as go
                densval = 0.0D0
                normvec = 0.0D0
                do l = 1, Npos
                    watpos = pos(l,:)
                    thisvec = apos - watpos
                    thisvec = thisvec - BoxL * anint(thisvec * iBoxL)
                    rvalsq = sum(thisvec*thisvec)
                    ! Check if beyond cut-off for this Gaussian
                    if (rvalsq >= (9.0D0*(smoothlen**2))) then
                        densfunc = 0.0D0
                        normfunc = 0.0D0
                    ! If not beyond cut-off, add to density field at this point
                    else
                        expterm = -rvalsq / (2.0D0*(smoothlen*smoothlen))
                        densfunc = exp(expterm) / ((2.0D0*pi*(smoothlen*smoothlen))**(1.5D0)) - shiftterm
                        ! Normal to interface is in direction of gradient (increasing density)
                        ! For minimization, usually go in direction of negative gradient - don't get confused
                        normfunc = -thisvec * (densfunc + shiftterm) / (smoothlen*smoothlen)
                    endif
                    densval = densval + densfunc
                    normvec = normvec + normfunc
                enddo
                ! Record the density value and normal at this grid point
                densvals(i, j, k) = densval
                ! Normalize the density gradient (surface normal vector)
                densnorms(i, j, k, :) = normvec / sqrt(sum(normvec*normvec))
            enddo
        enddo
    enddo
end subroutine


! Also want to be able to just get the value of the density field at specific points, not on a grid
! pos - water oxygen positions
! denspts - x,y,z positions at which to compute the density field and normal vectors
! BoxL - the box dimensions (for wrapping)
! smoothlen - the smoothing len (sigma)
! densvals - all density values at the grid points
! densnorms - all normal vectors at the grid points
subroutine WillardDensityPoints(pos, denspts, BoxL, smoothlen, densvals, densnorms, Npos, Npts)
    implicit none
    real(8), dimension(Npos, 3), intent(in) :: pos
    real(8), dimension(Npts, 3), intent(in) :: denspts
    real(8), dimension(3), intent(in) :: BoxL
    real(8), intent(in) :: smoothlen
    integer, intent(in) :: Npos, Npts
    real(8), dimension(Npts), intent(out) :: densvals
    real(8), dimension(Npts, 3), intent(out) :: densnorms
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), dimension(3) :: iBoxL, apos, watpos, normvec, thisvec, normfunc
    real(8) :: shiftterm, densval, rvalsq, expterm, densfunc
    integer :: i, j
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    shiftterm = exp(-9.0D0/2.0D0) / ((2.0D0*pi*(smoothlen*smoothlen))**(1.5D0))
    ! Loop over grid in each dimension to get all points
    do i = 1, Npts
        apos = denspts(i,:)
        ! Set density function and gradient of density to zero
        ! These are both functions of all the solvent coordinates, so will sum as go
        densval = 0.0D0
        normvec = 0.0D0
        do j = 1, Npos
            watpos = pos(j,:)
            thisvec = apos - watpos
            thisvec = thisvec - BoxL * anint(thisvec * iBoxL)
            rvalsq = sum(thisvec*thisvec)
            ! Check if beyond cut-off for this Gaussian
            if (rvalsq >= (9.0D0*(smoothlen**2))) then
                densfunc = 0.0D0
                normfunc = 0.0D0
            ! If not beyond cut-off, add to density field at this point
            else
                expterm = -rvalsq / (2.0D0*(smoothlen*smoothlen))
                densfunc = exp(expterm) / ((2.0D0*pi*(smoothlen*smoothlen))**(1.5D0)) - shiftterm
                ! Normal to interface is in direction of gradient (increasing density)
                ! For minimization, usually go in direction of negative gradient - don't get confused
                normfunc = -thisvec * (densfunc + shiftterm) / (smoothlen*smoothlen)
            endif
            densval = densval + densfunc
            normvec = normvec + normfunc
        enddo
        ! Record the density value and normal at this grid point
        densvals(i) = densval
        ! Normalize the density gradient (surface normal vector)
        densnorms(i, :) = normvec / sqrt(sum(normvec*normvec))
    enddo
end subroutine


! Given a set of points defining an interface and normal vectors at each point, along with solvent
! positions, the below finds which grid points are closest to each water, which water is closest
! to each grid point, and the number of waters within the magnitude of the normal vector from the
! interface
! pos - the positions of solvent
! gridpos - the position of points defining an interface
! gridnorm - the normal vectors at each interface point (should be of norm 1)
! cutoff - the cutoff for the water distance projected onto the nearest surface normal vector
! BoxL - the box dimensions
! watclose - the index of the surface point in gridpos closest to each water
! surfclose - the index of the water in pos closest to each gridpoint in gridpos
! numwater - the number of waters within the cutoff of the interface along the surface normals
! allwatdists - the distance of each water to the interface (useful for building density profiles)
subroutine InterfaceWater(pos, gridpos, gridnorm, cutoff, BoxL, watclose, surfclose, numwater, allwatdists, Npos, Ngridpos)
    implicit none
    real(8), dimension(Npos, 3), intent(in) :: pos
    real(8), dimension(Ngridpos, 3), intent(in) :: gridpos
    real(8), dimension(Ngridpos, 3), intent(in) :: gridnorm
    real(8), intent(in) :: cutoff
    real(8), dimension(3), intent(in) :: BoxL
    integer, intent(in) :: Npos, NgridPos
    integer, dimension(Npos), intent(out) :: watclose
    integer, dimension(Ngridpos), intent(out) :: surfclose
    integer, intent(out) :: numwater
    real(8), dimension(Npos), intent(out) :: allwatdists
    real(8), dimension(Ngridpos) :: griddists
    real(8), dimension(3) :: iBoxL, watpos, gpos, distvec, closenorm, normvec
    real(8) :: distsq, watdist, projectdist
    integer :: i, j
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    ! Assign large number to smallest distance for finding waters close to interface points
    griddists = 1000.0D0
    ! Loop over all solvent
    do i = 1, Npos
        watpos = pos(i,:)
        ! Reset smallest distance to large number and closenorm vector
        watdist = 1000.0D0
        closenorm = 0.0D0
        ! Loop over all interface points
        do j = 1, NgridPos
            gpos = gridpos(j,:)
            distvec = watpos - gpos
            distvec = distvec - BoxL * anint(distvec * iBoxL)
            distsq = sum(distvec*distvec)
            ! If this distance is smaller than that stored, update with closer point
            if (distsq < watdist) then
                watclose(i) = j
                watdist = distsq
                closenorm = gridnorm(j,:)
            endif
            if (distsq < griddists(j)) then
                surfclose(j) = i
                griddists(j) = distsq
            endif
        enddo
        ! Now that tried all interface points for this water, pick normal vector and see if water is counted
        ! If the normal vector, closenorm is normalized so that it's magnitude is 1, then the dot product
        ! of the water-surface vector with closenorm should be the projection we're looking for...
        ! i.e. it's some fraction of the total magnitude of the water-surface vector, or this magnitude
        ! multiplied by the cosine of the angle between this and the surface normal
        normvec = watpos - gridpos(watclose(i),:)
        normvec = normvec - BoxL * anint(normvec * iBoxL)
        projectdist = sum(normvec*closenorm)
        allwatdists(i) = projectdist
        if (projectdist <= cutoff) then
            numwater = numwater + 1
        endif
    enddo
end subroutine


!Given two vectors, computes the distance via a pseudo-spherical distance metric
!Have functions for 1-, 2-, and 3-body degrees of freedom
subroutine distanceMetric1B(vec1, vec2, rsq, sintw, dist)
    implicit none
    real(8), dimension(6), intent(in) :: vec1
    real(8), dimension(6), intent(in) :: vec2
    real(8), intent(in) :: rsq, sintw
    real(8), intent(out) :: dist
    real(8), dimension(6) :: diffs, dsqvec
    diffs = (vec2 - vec1)*(vec2 - vec1)
    dsqvec(1) = diffs(1)
    dsqvec(2) = diffs(2)
    dsqvec(3) = diffs(3)
    dsqvec(4) = rsq*diffs(4)
    dsqvec(5) = rsq*sin(vec2(4))*sin(vec1(4))*diffs(5)
    dsqvec(6) = rsq*sintw*diffs(6)
    dist = sqrt(sum(dsqvec))
end subroutine

    
subroutine distanceMetric2B(vec1, vec2, rsq, sintw, dist)
    implicit none
    real(8), dimension(12), intent(in) :: vec1
    real(8), dimension(12), intent(in) :: vec2
    real(8), intent(in) :: rsq, sintw
    real(8), intent(out) :: dist
    real(8), dimension(12) :: diffs, dsqvec
    diffs = (vec2 - vec1)*(vec2 - vec1)
    dsqvec(1) = diffs(1)
    dsqvec(2) = diffs(2)
    dsqvec(3) = diffs(3)
    dsqvec(4) = rsq*diffs(4)
    dsqvec(5) = rsq*sin(vec2(4))*sin(vec1(4))*diffs(5)
    dsqvec(6) = rsq*sintw*diffs(6)
    dsqvec(7) = diffs(7)
    dsqvec(8) = rsq*diffs(8)
    dsqvec(9) = rsq*diffs(9)
    dsqvec(10) = rsq*sin(vec2(9))*sin(vec1(9))*diffs(10)
    dsqvec(11) = rsq*sintw*diffs(11)
    dsqvec(12) = rsq*sintw*diffs(12)
    dist = sqrt(sum(dsqvec))
end subroutine


subroutine distanceMetric3B(vec1, vec2, rsq, sintw, dist)
    implicit none
    real(8), dimension(18), intent(in) :: vec1
    real(8), dimension(18), intent(in) :: vec2
    real(8), intent(in) :: rsq, sintw
    real(8), intent(out) :: dist
    real(8), dimension(18) :: diffs, dsqvec
    diffs = (vec2 - vec1)*(vec2 - vec1)
    dsqvec(1) = diffs(1)
    dsqvec(2) = diffs(2)
    dsqvec(3) = diffs(3)
    dsqvec(4) = rsq*diffs(4)
    dsqvec(5) = rsq*sin(vec2(4))*sin(vec1(4))*diffs(5)
    dsqvec(6) = rsq*sintw*diffs(6)
    dsqvec(7) = diffs(7)
    dsqvec(8) = rsq*diffs(8)
    dsqvec(9) = rsq*diffs(9)
    dsqvec(10) = rsq*sin(vec2(9))*sin(vec1(9))*diffs(10)
    dsqvec(11) = rsq*sintw*diffs(11)
    dsqvec(12) = rsq*sintw*diffs(12)
    dsqvec(13) = diffs(13)
    dsqvec(14) = vec2(13)*vec1(13)*diffs(14)
    dsqvec(15) = vec2(13)*vec1(13)*sin(vec2(14))*sin(vec2(14))*diffs(15)
    dsqvec(16) = rsq*diffs(16)
    dsqvec(17) = rsq*sin(vec2(16))*sin(vec1(16))*diffs(17)
    dsqvec(18) = rsq*sintw*diffs(18)
    dist = sqrt(sum(dsqvec))
end subroutine


! Calculates 3-dimensional histogram for tiplets of particles in Pos
! The bins are for distance 1-2, distance 1-3, and the 3-body angle
! All binning will start at zero for both distances and angles!
! Could change, but seems tricky to do so... just did simplest histogramming procedure possible
subroutine histrr3b(Pos, BoxL, distWidth, dNum, angWidth, aNum, histOut, NPos)
    implicit none
    real(8), dimension(NPos, 3), intent(in) :: Pos
    real(8), dimension(3), intent(in) :: BoxL
    real(8), intent(in) :: distWidth, angWidth
    integer, intent(in) :: dNum, aNum
    integer, intent(in) :: NPos
    real(8), dimension(dNum, dNum, aNum), intent(out) :: histOut
    integer :: i, j, k, dbin1, dbin2, abin
    real(8) :: dist1, dist2, ang
    real(8), dimension(3) :: refzero, refpos, pos1, pos2, distvec1, distvec2, iBoxL
    real*8, external :: CosAngle3
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    refzero = 0.
    histOut = 0.
    ! Loop over all triples and place in bins
    do i = 1, NPos
        refpos = Pos(i,:)
        do j = 1, NPos
            if (j == i) cycle
            pos1 = Pos(j,:)
            distvec1 = pos1 - refpos
            distvec1 = distvec1 - BoxL * anint(distvec1 * iBoxL)
            dist1 = sqrt(sum(distvec1*distvec1))
            ! Binning rule is to exclude lower edges, include upper edges
            ! In practice, should not make a big difference
            dbin1 = ceiling(dist1/distWidth)
            if (dbin1 > dNum) cycle
            do k = j+1, NPos
                if (k == i) cycle
                pos2 = Pos(k,:)
                distvec2 = pos2 - refpos
                distvec2 = distvec2 - BoxL * anint(distvec2 * iBoxL)
                dist2 = sqrt(sum(distvec2*distvec2))
                dbin2 = ceiling(dist2/distWidth)
                if (dbin2 > dNum) cycle
                ang = CosAngle3(distvec1, refzero, distvec2)
                abin = ceiling(ang/angWidth)
                if (abin > aNum) cycle
                histOut(dbin1, dbin2, abin) = histOut(dbin1, dbin2, abin) + 1
            enddo
        enddo
    enddo
end subroutine
                
            

