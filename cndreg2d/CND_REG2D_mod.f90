!This is essentially CND_REG.c in 2D rewritten for healpix maps.
!Provides cluster analysis of the excursion sets, Euler number
!and boundary length computation for each isolated cluster.
!l_map can be output for full cluster analysis.
!
!June 15th, 2011, D. Pogosyan    -    fix memory management
!June 10th, 2011, D. Pogosyan    -    allows choice what statistics to do
!May 25th, 2011, D. Pogosyan     -    new code for boundary length
!June 9th, 2009, D. Pogosyan     -    arrange as a single module
!May  4th, 2009, D. Pogosyan	    -    cosmetic changes
!July 4th, 2001, D. Pogosyan	    -    organized as a subroutine
!June, 2001, Dmitri Pogosyan     -    input-output changes
!April, 2001, Dmitri Pogosyan    -    first implementation for Healpix in F90
!
!Input:
!tmap(:)    - Healpix map
!nside1     - Healpix NSIDE
!ordering   - Healpix map ordering. Must be correct.
!level      - dT/T threshold to analyse
!
!Output:
!stats%n()  - sizes in pixels of isolated clusters above the level
!%g()       - 2 x euler number of each isolated cluster
!%s()       - number of each cluster side elements
!%sf()      - length of each cluster boundary
!nvoids     - number of isolated clusters (why do I call it so ??)
!
!Optional  CND_CNTRL(type CND_CNTRL_TYPE) structure has control variables
!CONNECT    - statistical (=CND_CNTRL%STAT or 0, this is default)
!             or exact (=CND_CNTRL%EXACT or 1)
!EXCURSION  - above (=CND_CNTRL%ABOVE or  1, this is default)
!             or below (=CND_CNTRL%BELOW or -1) the threshold
!
!Note:
!for 'statistical' cluster connectivity,
!only sum(genus)/4 over all clusters makes real sense.
!please divide by 4 only after sum

module CND_REG2D_mod
use healpix_types
use pix_tools

implicit none
private

type EXCRS_STAT
     integer(I4B),pointer    :: n(:) !sizes in pixels of isolated clusters above the level
     integer(I4B),pointer    :: s(:) !number of each cluster side elements
     integer(I4B),pointer    :: sm(:)
     real(DP),pointer        :: sf(:) !length of each cluster boundary
     integer(I4B),pointer    :: g(:) !2 x euler number of each isolated cluster
     logical                 :: do_genus=.true., do_sides=.false., do_length=.true.
end type EXCRS_STAT

type CND_CNTRL_TYPE
     integer                 :: RING=1, NESTED=2 !HEALPix map ordering
     integer                 :: STAT=0, EXACT=1 !'statistical' or 'exact' cluster connectivity
     integer                 :: ABOVE=1, BELOW=-1 !excursion logic
     integer                 :: CONNECT=0, EXCURSION=1
     real(DP)                :: MASK_VALUE=-1.d29
end type CND_CNTRL_TYPE

integer,parameter            :: UNCHECKED=-1, OUTSIDE=0, MASKED=-2, BORDER=-4, MASKED_BORDER=-6 !INSIDE > 0
integer(I4B)                 :: NSIDE, NSIDESQ, NPIX, ORDERING, MAXSTACK=0
integer(I4B),allocatable     :: l_map(:)
integer(I4B),allocatable     :: pixstack(:)
integer(I1B),dimension(0:23) :: vertex,vertex_exct,vertex_stat

! The indexing corresponds to circular labelling of pixels around the vertice
! choice of -2 for diagonal links is 'exact'. Choice of 0 is for 'statistics'
data vertex_exct / 0,1,1,0,1,-2,0,-1,1,0,-2,-1,0,-1,-1,0,0,1,1,0,1,0,0,1 /
data vertex_stat / 0,1,1,0,1, 0,0,-1,1,0, 0,-1,0,-1,-1,0,0,1,1,0,1,0,0,1 /

public                       :: CND_REG, CND_CNTRL_TYPE, EXCRS_STAT

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  subroutine CND_REG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CND_REG (tmap,nside1,ordering1,level,stats,nvoids,CND_CNTRL_IN)

implicit none
real(DP),intent(in),dimension(0:)     :: tmap
real(DP),intent(in)                   :: level
integer(I4B),intent(in)               :: nside1, ordering1
integer(I4B),intent(out)              :: nvoids
type(EXCRS_STAT),intent(inout)        :: stats
type(CND_CNTRL_TYPE),optional         :: CND_CNTRL_IN
integer(I4B)                          :: i, i_nest, nabove, k
integer(I4B)                          :: NV, NB
type(CND_CNTRL_TYPE)                  :: CND_CNTRL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Data read-in and setup

if (.not.stats%do_genus .and. .not.stats%do_sides .and. .not.stats%do_length) then
  write(0,*)'No job was requested, return'
  return
endif

if (PRESENT(CND_CNTRL_IN)) CND_CNTRL = CND_CNTRL_IN

if (ordering1/=CND_CNTRL%NESTED .and. ordering1/=CND_CNTRL%RING) then
  write(0,*)'HEALPix ordering',ordering1,'unrecognized'
endif

NSIDE = nside1
NSIDESQ = nside1**2
NPIX = nside2npix(NSIDE)
ORDERING = ordering1

!Allocate and initialized the index map
allocate(l_map(0:NPIX-1))
l_map = OUTSIDE

nabove=0
do i=0,NPIX-1
i_nest=i
if (tmap(i) < CND_CNTRL%MASK_VALUE) then
  if (ordering1 == CND_CNTRL%RING) call ring2nest(NSIDE,i,i_nest)
  l_map(i_nest) = MASKED
  elseif (((CND_CNTRL%EXCURSION==CND_CNTRL%ABOVE).and.(tmap(i)>=level)).or. &
          ((CND_CNTRL%EXCURSION==CND_CNTRL%BELOW).and.(tmap(i)<level))) then
    if (ordering1 == CND_CNTRL%RING) call ring2nest(NSIDE,i,i_nest)
    l_map(i_nest) = UNCHECKED
    nabove = nabove+1
  endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Cluster Analysis
  
NV=min(nabove+1,NPIX/4)             ! max number of the connected regions
allocate(stats%n(1:NV))           ! will contain volume of clusters

NB=nabove                   ! Stack size NB should be at least equal number
allocate(pixstack(0:NB))  ! of border points for largest cluster.
                              ! Which can be as large as all points nabove
stats%n=0

k=0
do i=0,NPIX-1                     ! Main loop begins
  if (l_map(i) == UNCHECKED) then ! we hit cluster interior
    k=k+1                         ! This cluster will be  k+1
    l_map(i) = k                  ! Mark pixel i as belonging to cluster k
    pixstack(0) = i               ! Put pixel i on stack
    if(k.gt.NV) then
      print*, 'k>NV!!!, should never happen, stop'
      stop
    else
      stats%n(k) = regions(k)  ! Search inside cluster k and get its size
    endif
  endif
enddo

nvoids=min(k,NV)
deallocate(pixstack)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Genus determination

if (stats%do_genus) then
  allocate(stats%g(1:nvoids)) ! will contain genus of clusters
  if (CND_CNTRL%CONNECT == CND_CNTRL%STAT) then
    vertex = vertex_stat
  else if (CND_CNTRL%CONNECT == CND_CNTRL%EXACT) then
    vertex = vertex_exct
  else
    write(0,*)'requested connectivity type ',CND_CNTRL%CONNECT ,'unknown'
    write(0,*)'default to STATISTICAL'
    vertex = vertex_stat
  endif
  stats%g=0
  do i=0,NPIX-1 ! Loop over vertices.
                ! need to check INSIDE or BORDER cells only
                ! to get all the vertices
    if (l_map(i) > OUTSIDE .or. l_map(i) <= BORDER) then
      call update_genus(i,stats%g)
    endif
  enddo
  ! Two vertices - North and South poles - are extra
  call update_genus_pole(0,stats%g) ! North pole in ring notation
  call update_genus_pole(NPIX-4,stats%g) ! South pole in ring notation
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Number of border segments determination

if (stats%do_sides) then
  allocate(stats%s(1:nvoids)) ! contains number of cluster sides
  allocate(stats%sm(1:nvoids)) ! contains count of masked boudary
  stats%s=0
  stats%sm=0
  do i=0,NPIX-1 ! Loop over pixels
                ! To get all the border edges
    if (l_map(i) == BORDER) then ! need to check BORDER pixels
      call update_sides(i,stats%s)
    else if (l_map(i) == MASKED_BORDER) then
      call update_sides(i,stats%sm) ! or BORDER+MASKED only
    endif
  enddo
endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Border length determination

if (stats%do_length) then
  allocate(stats%sf(1:nvoids)) ! contains surface length of clusters
  stats%sf=0.d0
  do i=0,NPIX-1 ! Loop over vertices
                ! MASKED BORDER is not included
    if (l_map(i) > OUTSIDE .or. l_map(i) == BORDER) then
      call update_surf(tmap,level,i,stats%sf)
    endif
  enddo
endif
deallocate(l_map)
  
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  auxiliary functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function regions (k)
integer k,new_on_stack,istack
istack = 0
regions = 1
do while (istack >= 0)
new_on_stack = neibr(k,istack)
regions = regions + new_on_stack
istack = istack + new_on_stack - 1
enddo
return
end function


integer function neibr(k,istack)
integer k,istack
integer i,new_on_stack,ipix
integer(I4B) center_pix,neighbour_list(8),n_neighbours
center_pix = pixstack(istack)
call neighbours_nest(NSIDE,center_pix,neighbour_list,n_neighbours)
neibr = 0
do i = 1,n_neighbours
ipix=neighbour_list(i)
if (l_map(ipix) == UNCHECKED) then
pixstack(istack+neibr) = ipix
if (istack+neibr > MAXSTACK) MAXSTACK=istack+neibr
l_map(ipix) = k
neibr = neibr+1
else if (l_map(ipix) == OUTSIDE) then
l_map(ipix) = BORDER
else if (l_map(ipix) == MASKED)  then
l_map(ipix) = MASKED_BORDER
endif
enddo
return
end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This is main function for genus analysis. Sets code of the vertex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE update_genus(vert,genus)
INTEGER(I4B), INTENT(IN)    :: vert
INTEGER(I4B), INTENT(INOUT) :: genus(:)
INTEGER(I4B) neighbour_list(8),nneigh
INTEGER(I4B) bp1,bp2,bp3,set_code
INTEGER(I1B) code

set_code = OUTSIDE
code = 0
! We select the vertice as the one produced by: :
! vert + three last neighbours (for 4 edge vertice)
! vert + two last neighbour    (for 3 edge vertice)
! Order is: vert, last, last-1, last-2 (if needed), corresponding
! to first 4 bits in 'code', 5th bit is 1 for 3 edge vertice
!
! This emuneration points to each pixel's eastmost vertex
! (Which is much easier than sourthern or northern !!!!
! Eastmost vertex is 3 edge if vert, last, last-1 belong
! to three different basic pixels

call neighbours_nest(NSIDE,vert,neighbour_list,nneigh)

if (l_map(vert)  > OUTSIDE) then
code = IBSET(code,0)
set_code = l_map(vert)
endif

if (l_map(neighbour_list(nneigh)) > OUTSIDE)  then
code = IBSET(code,1)
set_code = l_map(neighbour_list(nneigh))
endif

if (l_map(neighbour_list(nneigh-1)) > OUTSIDE) then
code = IBSET(code,2)
set_code = l_map(neighbour_list(nneigh-1))
endif

if (nneigh <= 7) then        ! Special case, our vertex can be 3 edge
bp1=vert/NSIDESQ
bp2=neighbour_list(nneigh)/NSIDESQ

if (bp1 /= bp2) then       ! 3 edge is still possible
bp3=neighbour_list(nneigh-1)/NSIDESQ
if ((bp1 /= bp3).and.(bp2 /= bp3)) then
code = IBSET(code,4)              ! yes, it is
code = IBCLR(code,3)              ! 4th beighbour is not used
endif
endif
endif

if (.not.BTEST(code,4)) then         !4th neighbour checked for 4 edge
if (l_map(neighbour_list(nneigh-2)) > OUTSIDE) then
code = IBSET(code,3)
set_code = l_map(neighbour_list(nneigh-2))
endif
endif

if (set_code > OUTSIDE) genus(set_code) = genus(set_code)+vertex(code)
RETURN
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Special version of set code for South and North poles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE update_genus_pole (pole,genus)
INTEGER(I4B), INTENT(IN)    :: pole
INTEGER(I4B), INTENT(INOUT) :: genus(:)

INTEGER(I1B) code
INTEGER(I4B) set_code_pole,inest,i

set_code_pole = OUTSIDE
code = 0

do i=0,3
call ring2nest(NSIDE,pole+i,inest)
if (l_map(inest) > OUTSIDE) then
code = IBSET(code,i)
set_code_pole = l_map(inest)
endif
enddo

if (set_code_pole > OUTSIDE) then
genus(set_code_pole)=genus(set_code_pole)+vertex(code)
endif
RETURN
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This is function for border analysis. Counts border segments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE update_sides (vert,s)
INTEGER(I4B), INTENT(IN)    :: vert
INTEGER(I4B), INTENT(INOUT) :: s(:)
INTEGER(I4B) neighbour_list(8),nneigh
INTEGER(I4B) i,k,bp1,bp2,bp3

! We are at a border pixel and need to count boundary edges,
! but not boundary vertices. In case of a regular 8 neighbour pixel,
! the neighbours with common edge ! are number 2,4,6,8

! If nneigh=7, pixel has 3-edge vertex. (NSIDE=1 can have nneigh < 7 !)
! The vertex is 3 pixel, if its neighbours belong to different
! 12 basis pixels

call neighbours_nest(NSIDE,vert,neighbour_list,nneigh)

if (nneigh == 8) then ! Straighforward loop, separate for speed

do i=2,nneigh,2
k=l_map(neighbour_list(i))
if (k > OUTSIDE) s(k)=s(k)+1
enddo

else                    ! Special case, our cell has 3 edge vertex
! this code is general, but needed for 8 cells only
bp1=vert/NSIDESQ

i=nneigh
do while (i >= 2)
k=l_map(neighbour_list(i))
if (k > OUTSIDE) s(k)=s(k)+1
! if three edge vertex, reset counter
bp3=neighbour_list(i-1)/NSIDESQ
if (bp1 /= bp3) then
bp2=neighbour_list(i)/NSIDESQ
if ((bp2 /= bp3).and.(bp1 /= bp2)) i=i+1
endif
i=i-2
enddo

endif

RETURN
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This is new function for surface analysis. Sets surface length in radians
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE update_surf(tmap,level,vert,s)

REAL(DP),     intent(in), dimension(0:) :: tmap
REAL(DP),     intent(in)    :: level
INTEGER(I4B), intent(in)    :: vert
REAL(DP),     intent(inout) :: s(:)

TYPE(CND_CNTRL_TYPE)        :: CND_CNTRL
INTEGER(I4B) neighbour_list(8),best_pix(4),nneigh,vneigh
INTEGER(I4B) bp1,bp2,bp3,k,i,i_ring,cross_count
INTEGER(I1B) code

REAL(DP), DIMENSION(3,4)    :: vertex_xyz,pix_xyz
REAL(DP), DIMENSION(3,2)    :: cross_xyz
REAL(DP), DIMENSION(4)      :: field
REAL(DP), DIMENSION(3)      :: bfit
REAL(DP), DIMENSION(4,3)    :: AA
REAL(DP)                    :: ds,ds1,DDOT

! We use the same code to classify vertex as for genus
k = OUTSIDE
code = 0
! We select the vertice as the one produced by: :
! vert + three last neighbours (for 4 edge vertice)
! vert + two last neighbour    (for 3 edge vertice)
! Order is: vert, last, last-1, last-2 (if needed), corresponding
! to first 4 bits in 'code', 5th bit is 1 for 3 edge vertice
!
! This emuneration points to each pixel's eastmost vertex
! (Which is much easier than sourthern or northern !!!!
! Eastmost vertex is 3 edge if vert, last, last-1 belong
! to three different basic pixels
!
! If any of the neighbors were masked, the vertice doesn't contribute
! we know in advance that 'vert' was not masked

call neighbours_nest(NSIDE,vert,neighbour_list,nneigh)

if (l_map(vert)  > OUTSIDE) then
code = IBSET(code,0)
k = l_map(vert)
endif

if (l_map(neighbour_list(nneigh)) > OUTSIDE)  then
code = IBSET(code,1)
k = l_map(neighbour_list(nneigh))
elseif (l_map(neighbour_list(nneigh)) == MASKED_BORDER) then
return
endif

if (l_map(neighbour_list(nneigh-1)) > OUTSIDE) then
code = IBSET(code,2)
k = l_map(neighbour_list(nneigh-1))
elseif (l_map(neighbour_list(nneigh-1)) == MASKED_BORDER) then
return
endif

if (nneigh <= 7) then        ! Special case, our vertex can be 3 edge
bp1=vert/NSIDESQ
bp2=neighbour_list(nneigh)/NSIDESQ

if (bp1 /= bp2) then       ! 3 edge is still possible
bp3=neighbour_list(nneigh-1)/NSIDESQ
if ((bp1 /= bp3).and.(bp2 /= bp3)) then
code = IBSET(code,4)
code = IBCLR(code,3)              ! 4th beighbour is not used
endif
endif
endif

if (.not.BTEST(code,4)) then         !4th neighbour checked for 4 edge
if (l_map(neighbour_list(nneigh-2)) > OUTSIDE) then
code = IBSET(code,3)
k = l_map(neighbour_list(nneigh-2))
elseif (l_map(neighbour_list(nneigh-2)) == MASKED_BORDER) then
return
endif
endif

! -------------------------------------------------------------------
! Choose what pixels to use to approximate the field by the plane
! Based on "code" we choose 3 or 4, treating "saddle point" case specially

call select_pixels_for_fit(code,vert,neighbour_list,nneigh,best_pix,vneigh)

if (vneigh > 0) then  ! vneigh=0 means we are not at the boundary

! a) Get temperature values in neghbouring pixels and shift by level.
!    tmap could have been in ring notation, so we'll account for it
if (ORDERING == CND_CNTRL%RING) then
do i=1,vneigh
call nest2ring(NSIDE,best_pix(i),i_ring)
field(i) = tmap(i_ring) - level
enddo
else
do i=1,vneigh
field(i) = tmap(best_pix(i)) - level
enddo
endif

! b) first just get coordinates of the vertex and vert pixel
!    Vertex is the EAST(4th) one of the orig pixel "vert"
call pix2vec_nest(NSIDE,vert,pix_xyz(:,1),vertex_xyz)

!    vert is first pixel if vneigh=4, and who knows where if vneigh=3
!    so now get coordinates of three pixels (1:3 or 2:4)
do i=vneigh-2,vneigh
call pix2vec_nest(NSIDE,best_pix(i),pix_xyz(:,i))
enddo

! c) Project all pixel positions onto the plane tangent at vertex
!    and move the coordinate origin to the vertex
do i=1,vneigh
pix_xyz(:,i)=pix_xyz(:,i)/DDOT(3,pix_xyz(:,i),1,vertex_xyz(:,4),1) - vertex_xyz(:,4)
enddo

! d) Linearly approximate the field in the tangent plane
!    f = f_0 + a (x_1 . x) + b (x_2 . x) where I'll take vert and
!    it's last neighbour positions for basis vectors x_1 and x_2
!    The following has a structure to include 4 pixel noise covariance
!    which is currently set by hand to unit matrix

if (code == 5 .or. code == 10) vneigh = 3 ! for saddle, first triplet
call planar_fit(vneigh,field,pix_xyz,AA,bfit)

! e) Solve for intersection of the surface line with boundaries
!    around the vertex (that connect nominal pixel centers)
!    There should be no more (barring degeneracy) and
!    certainly no less than 2 of those

call intersections(vneigh,pix_xyz,AA,bfit,cross_xyz)

! Last step, get the length. We should map intersection points back onto
! sphere and calculate angular distance between them, but
! explicit mapping onto sphere can be avoided by using Healpix angdist

cross_xyz(:,1) = cross_xyz(:,1) + vertex_xyz(:,4)
cross_xyz(:,2) = cross_xyz(:,2) + vertex_xyz(:,4)
call angdist(cross_xyz(:,1),cross_xyz(:,2),ds)

if (code == 5 .or. code == 10) then ! saddle points has second segment
! get it replacing pix 1 with 4
field(1)     = field(4)
pix_xyz(:,1) = pix_xyz(:,4)
call planar_fit(vneigh,field,pix_xyz,AA,bfit)
call intersections(vneigh,pix_xyz,AA,bfit,cross_xyz)
cross_xyz(:,1) = cross_xyz(:,1) + vertex_xyz(:,4)
cross_xyz(:,2) = cross_xyz(:,2) + vertex_xyz(:,4)
call angdist(cross_xyz(:,1),cross_xyz(:,2),ds1)
ds = ds + ds1
endif

! Time to increment the lengths of the "k" region
s(k) = s(k) + ds
endif

RETURN
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine select_pixels_for_fit(code,vert,neighbour_list,nneigh,best_pix,vneigh)

INTEGER(I1B), intent(in)    :: code
INTEGER(I4B), intent(in)    :: vert,nneigh,neighbour_list(8)
INTEGER(I4B), intent(out)   :: vneigh, best_pix(4)

INTEGER(I4B)                :: i

if (code==0 .or. code==15 .or. code==16 .or. code==23) then
vneigh=0    ! not a boundary, nothing to be done
return
endif

if (code > 16) then                   ! 3 sided vertex, use all
best_pix = (/ 1, 2, 3, -1 /)
vneigh = 3
else if (code==1 .or. code==14) then  !One above, three below, or reverse
best_pix = (/ 1, 2, 4, -1 /)
vneigh = 3
else if (code==2 .or. code==13) then
best_pix = (/ 2, 3, 1, -1 /)
vneigh = 3
else if (code==4 .or. code==11) then
best_pix = (/ 3, 4, 2, -1 /)
vneigh = 3
else if (code==8 .or. code==7) then
best_pix = (/ 4, 1, 3, -1 /)       ! First two pixels never coplanar
vneigh = 3
else if (code==5 .or. code==10) then  ! saddle point, use all 4 points
best_pix = (/ 1, 2, 4, 3 /)        ! but in two runs of 3
vneigh = 4
else                           ! Use all 4 pixels, applies to code 3,6,9,12
best_pix = (/ 1, 2, 3, 4 /)
vneigh = 4
endif

do i=1,4
if (best_pix(i) /= -1) then
if (best_pix(i) == 1) then
best_pix(i) = vert
else
best_pix(i) = neighbour_list(nneigh-best_pix(i)+2)
endif
endif
enddo

return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine planar_fit(vneigh,field,pix_xyz,AA,bfit)

INTEGER(I4B), intent(in)                :: vneigh
REAL(DP), intent(in), DIMENSION(3,4)    :: pix_xyz
REAL(DP), intent(in), DIMENSION(4)      :: field
REAL(DP), intent(out),DIMENSION(3)      :: bfit
REAL(DP), intent(out),DIMENSION(4,3)    :: AA

REAL(DP), DIMENSION(4,4)    :: CNpp
REAL(DP), DIMENSION(4,3)    :: CNA
REAL(DP), DIMENSION(3,3)    :: AtCNA
REAL(DP), DIMENSION(3)      :: WORK
INTEGER(I4B), DIMENSION(3)  :: IPIV

INTEGER(I4B)                :: i,i1,INFO
REAL(DP)                    :: DDOT

DATA CNpp /1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.0d0/

do i=1,vneigh
AA(i,1) = 1.d0
AA(i,2) = DDOT(3,pix_xyz(:,i),1,pix_xyz(:,1),1)
AA(i,3) = DDOT(3,pix_xyz(:,i),1,pix_xyz(:,2),1)
enddo
call DSYMM('L','L',vneigh,3,1.d0,CNpp,4,AA,4,0.d0,CNA,4)
call DGEMM('T','N',3,3,vneigh,1.d0,AA,4,CNA,4,0.d0,AtCNA,3)
call DGEMV('T',vneigh,3,1.d0,CNA,4,field,1,0.d0,bfit,1)
call DSYSV('L',3,1,AtCNA,3,IPIV,bfit,3,WORK,3,INFO)

return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine intersections(vneigh,pix_xyz,AA,bfit,cross_xyz)

INTEGER(I4B), intent(in)                :: vneigh
REAL(DP), intent(in), DIMENSION(3,4)    :: pix_xyz
REAL(DP), intent(in), DIMENSION(4,3)    :: AA
REAL(DP), intent(in), DIMENSION(3)      :: bfit
REAL(DP), intent(out),DIMENSION(3,2)    :: cross_xyz

INTEGER(I4B)                :: i,i1,cross_count
REAL(DP)                    :: bi,bi1,alpha

cross_count = 0

do i=1,vneigh
if (i == vneigh) then
i1=1
else
i1=i+1
endif
bi = bfit(2)*AA(i,2)+bfit(3)*AA(i,3)
bi1= bfit(2)*AA(i1,2)+bfit(3)*AA(i1,3)
alpha = -(bfit(1) + bi)/(bi1 - bi)

if (alpha >= 0.d0 .and. alpha <= 1.d0) then
cross_count = cross_count + 1
cross_xyz(:,cross_count) = (1.d0-alpha)*pix_xyz(:,i) + alpha*pix_xyz(:,i1)
if (cross_count == 2) exit
endif
enddo

if (cross_count < 2) then
write(0,*)'Did not find intersection where expected, stop'
stop
endif

return
end subroutine

end module CND_REG2D_mod
