!Minkowski functionals for single scalar map
!it reads only the first map from a given HEALPix fits file

!authored by Anne Ducout and Dmitri Pogosyan
!adjusted by Jiaxin Wang

!input arguments of the main routine:
!
! 1st, path to the HEALPix signal map
!
! 2nd, path to the HEALPix mask map
!
! 3rd, HEALPix Nside for input maps
!
! 4th, number of MFs thresholds


program main
  
use healpix_types
use fitstools
use CND_REG2D_mod

implicit none
!warning, I4B may not be sufficient
character*80                       :: p_chr !number of thresholds char
integer(I4B)                       :: p !number of thresholds
integer(I4B),parameter             :: nmap = 1 !single map to be read
character*80                       :: nside_chr !HEALPix Nside char
integer(I4B)                       :: nside, npix !HEALPix Nside and Npix
integer(I4B)                       :: umpix !number of un-masked pixels
real(DP)                           :: mean, variance, rms, minval, maxval !basic stastics of map
real(DP),parameter                 :: maxnu = 3.5d0 !mamimum MFs threshold value, in units of RMS
real(DP),dimension(:),allocatable  :: mapdata !memory handle for single HEALPix map
real(DP),dimension(:,:),allocatable:: mapmask, mapfits
character*80,parameter             :: nameout = 'mfs.dat' !MFs output file path
character*80                       :: mapname !mapname, path to the target scalar map
character*80                       :: maskname !maskname, path to the mask map
logical                            :: anynull !nuisance parameter for HEALPix read map function
integer(I4B)                       :: i
real(DP)                           :: nullval !nuisane parameter for HEALPix read map function
real(DP),parameter                 :: badvalue = -1.6375000d30 !HEALPix bad value

write(*,*)'info: in main routine'
  
!===============================================================================
!            Parsing parameters (customized zone)
!===============================================================================
  
!make sure the right number of inputs have been provided
if(command_argument_count().ne.4) then
write(*,*)'4 command-line arugments required: map path, mask path, HEALPix Nside, number of thresholds'
stop
endif

call get_command_argument(1,mapname) !path to target scalar map
call get_command_argument(2,maskname) !path to mask map
call get_command_argument(3,nside_chr)
call get_command_argument(3,p_chr)

read(nside_chr,*) nside !convert input char to integer
npix = 12*nside*nside !calculate Npix from Nside
read(p_chr,*) p !number of thresholds

!===============================================================================
!            Reading map
!===============================================================================

allocate(mapdata(0:npix-1))
allocate(mapfits(0:npix-1,1:nmap)) !follows HEALPix FITS map handler definition
!HEALPix function
!https://healpix.sourceforge.io/html/sub_read_bintab.htm
!input arguments: #1 FITS-file name, #3 total pixel size
!output arguments: #2 map handler, #4 value of missing pix, #5 if missing pix
call read_bintab(mapname,mapfits,npix,nmap,nullval,anynull)
!map takes over mapfits
mapdata = mapfits(:,1)
!delete mapfits
deallocate(mapfits)
write(*,*)'info: target map read'

!===============================================================================
!  Importation of the galactic mask, ring
!===============================================================================

allocate(mapmask(0:npix-1,1:nmap))
call read_bintab(maskname,mapmask,npix,nmap,nullval,anynull) !HEALPix function
where(mapmask(:,1)==0.) !set masked area by HEALPix.bad_value
mapdata = badvalue
end where
umpix=count(mapmask(:,1) > 0) !counting unmasked pixels
write(*,*)'info: mask map applied'
deallocate(mapmask)

!===============================================================================
!       Properties and normalisation of the maps (excluding masked pixels)
!===============================================================================

mean = sum(mapdata(0:npix-1), mask = mapdata>badvalue)/umpix
write(*,*)'info: partial sky mean =',mean
variance = 0d0
minval = mapdata(0)
maxval = mapdata(0)

do i = 0,npix-1
  if (mapdata(i)>badvalue) then !un-masked area
    variance = variance+(mapdata(i)-mean)**2
    !update min and max values
    if (mapdata(i)<minval) minval = mapdata(i)
    if (mapdata(i)>maxval) maxval = mapdata(i)
  endif
enddo
write(*,*)'info: partial sky min =',minval,', max =',maxval

variance = variance/(umpix-1)
write(*,*)'info: partial sky var =',variance

rms = sqrt(variance)
write(*,*)'info: partial skay std =',rms

!===============================================================================
!            Minkowski functionals of the map
!===============================================================================

!call the subroutine defined in the following
call mink(mapdata,nside,npix,nameout,p,rms,maxnu,minval,mean,umpix)
deallocate(mapdata)

end program main


subroutine mink(tmap,nside,npix,nameout,p,rms,maxnu,minval,meanmap,umpix)

!Minkowski functionals of the map. Loop on threshold values, functionals are
!calculated at each threshold and written on an ASCII file.
!arguments:
!input
!#1 tmap, target map (masked)
!#2 nside, target map HEALPix Nside
!#2 nameout, output plain file path
!#3 p, number of thresholds
!#4 rms, target map std (partial sky)
!#5 maxnu, threshold upper limit, in std unit
!#6 minval, target map min (parital sky)
!#7 meanmap, target map mean value (normalized)
!#8 umpix, un-masked pixel number
!#9 nside, input map and mask nside

use healpix_types
use fitstools
use CND_REG2D_mod
use pix_tools
 
implicit none
integer(I4B),intent(in)          :: p, umpix
integer(I4B),intent(in)          :: nside, npix
character*80,intent(in)          :: nameout
real(DP), dimension(0:npix-1)    :: tmap
real(DP)                         :: level, rms, maxnu, meanmap, minval
integer(I4B)                     :: j, k
!TYPE EXCRS_STST defined in CND_REG2D_mod.f90
type(EXCRS_STAT)                 :: stats
integer(I4B)                     :: nvoids, ordering
integer(I4B)                     :: result_volume, result_genus
integer(I4B)                     :: result_sides, result_maskedsides
real(DP)                         :: result_length
real(DP)                         :: v0, v1, v2
real(DP)                         :: fsky !sky coverage
!TYPE CND_CNTRL_TYPE defined in CND_REG2D_mod.f90
TYPE(CND_CNTRL_TYPE)             :: CND_CNTRL

!===============================================================================
!            Initialisation
!===============================================================================

write(*,*)'info: in subroutine mink'

ordering = CND_CNTRL%RING !ordering =  CND_CNTRL%NESTED
CND_CNTRL%CONNECT = CND_CNTRL%STAT !CND_CNTRL%CONNECT = CND_CNTRL%EXACT

fsky = umpix
fsky = fsky/npix

open(1,file=nameout)
write(1,*)'# info: sky fraction =',fsky
write(1,*)'# info: initial RMS of the map =',rms
write(*,*)'info: mink initialized'

!===============================================================================
!            Loop on thresholds
!===============================================================================

stats%do_sides =.false.

!loop through threshold levels
do j=0,p-1
  ! level in signal units, maxnu is the maximum sigma span of std
  level = meanmap-maxnu*rms+j*2.d0*maxnu*rms/(p-1)
  if (CND_CNTRL%CONNECT == CND_CNTRL%STAT) then
    if (level < meanmap) then
      CND_CNTRL%EXCURSION = CND_CNTRL%BELOW
    else
      CND_CNTRL%EXCURSION = CND_CNTRL%ABOVE
    endif
  endif

  !call CND_REG subroutine defined in CND_REG2D_mod.f90
  call CND_REG(tmap,nside,ordering,level,stats,nvoids,CND_CNTRL)

  result_genus = 0
  result_sides = 0
  result_maskedsides = 0
  result_volume = 0
  result_length = 0.d0

  do k = 1,nvoids
    result_volume = result_volume+stats%n(k)
    if (stats%do_genus) result_genus = result_genus+stats%g(k)
    if (stats%do_sides) then
      result_sides = result_sides+stats%s(k)
      result_maskedsides = result_maskedsides+stats%sm(k)
    endif
    if (stats%do_length) result_length = result_length+stats%sf(k)
  enddo

  if ( CND_CNTRL%CONNECT == CND_CNTRL%STAT ) then
    if ( level < meanmap ) then
      result_volume = npix - result_volume
      result_genus  = 8 - result_genus
    endif
  endif

  !v0: first Minkowski funtional, volume (area)
  if ( level<meanmap ) then
    v0 = result_volume-(npix-umpix)
  else
    v0 = result_volume
  endif
  v0 = v0/umpix
  !v1: second Minkowski funtional, perimeter
  v1=result_length/(4.d0*Pi*fsky)/4.d0
  !v2: third Minkowski funtional, genus
  v2 = result_genus
  v2 = v2/(4.d0*Pi*fsky)/4.d0

  write(1,'(e10.3,2X,e10.3,2X,e10.3,2X,e10.3,2X,i7)')level,v0,v1,v2,nvoids

  ! clear stats arrays
  deallocate(stats%n)
  if (stats%do_sides) then
    deallocate(stats%s)
    deallocate(stats%sm)
  endif
  if (stats%do_length) deallocate(stats%sf)
  if (stats%do_genus) deallocate(stats%g)

enddo
close(1)

! Get the properties of the manifold. Can be used for correction
stats%do_genus=.true.
stats%do_sides=.true.
stats%do_length=.false.

CND_CNTRL%EXCURSION = CND_CNTRL%ABOVE

call CND_REG(tmap,nside,ordering,minval-rms,stats,nvoids,CND_CNTRL)

result_genus=0
result_sides=0
result_maskedsides=0
result_volume=0

do k=1,nvoids
  result_volume =result_volume+stats%n(k)
  result_genus = result_genus+stats%g(k)
  result_sides = result_sides+stats%s(k)
  result_maskedsides = result_maskedsides+stats%sm(k)
enddo

write(*,*)'Properties of the manifold',result_genus/4.,result_volume,result_sides,result_maskedsides,nvoids

! clear stats arrays
deallocate(stats%n)
if (stats%do_sides) then
  deallocate(stats%s)
  deallocate(stats%sm)
endif
if (stats%do_length) deallocate(stats%sf)
if (stats%do_genus) deallocate(stats%g)

end subroutine mink
