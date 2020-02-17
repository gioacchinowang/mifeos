!Minkowski functionals for single scalar map

!shared by Dr. Anne Ducout
!modified by Dr. Jiaxin Wang


program main
  
  USE healpix_types
  USE fitstools
  USE CND_REG2D_mod
  
  implicit none
  !warning, I4B may not be sufficient
  integer(I4B)                       :: p !number of thresholds
  integer(I4B),parameter             :: nmap=1 !single map to be read
  character*80                       :: nside_c !HEALPix Nside char
  integer(I4B)                       :: nside, npix !HEALPix Nside and Npix
  integer(I4B)                       :: i,umpix
  real(DP)                           :: mean,variance,minval,maxval
  real(DP)                           :: ecartt,nullval
  real(DP)                           :: badvalue
  !maxnu, maximum threshold (symetric)
  real(DP)                           :: maxnu,meanmap
  real(DP),dimension(:),allocatable  :: map
  real(DP),dimension(:,:),allocatable:: mapmask,mapfits
  character*80                       :: nameout
  !mapname, path to the target scalar map
  !maskname, path to the mask map
  character*80                       :: mapname,maskname
  logical                            :: anynull

  write(*,*)'info: in main routine'
  !===========================================================================
  !            Parameters (customized zone)
  !=========================================================================== 
  
  !make sure the right number of inputs have been provided
  if(COMMAND_ARGUMENT_COUNT().ne.3) then
    write(*,*)'error, three command-line arugments required: map name, mask name, HEALPix Nside'
    stop
  endif

  call GET_COMMAND_ARGUMENT(1,mapname)
  call GET_COMMAND_ARGUMENT(2,maskname)
  call GET_COMMAND_ARGUMENT(3,nside_c)

  !mapname='map_cmb_ns64.fits' !path to target scalar map
  !maskname='mask_gal_fsky0.80_ns64.fits' !path to mask map
  read(nside_c,*) nside
  npix = 12*nside*nside

  p=100 !number of thresholds
  maxnu=3.5d0 !maxinum threshold val, in units of sigma(map)

  badvalue=-1.6375000d30 !HEALPix.bad_value

  !===========================================================================
  !            Reading map
  !===========================================================================

  allocate(map(0:npix-1))
  !HEALPix FITS map handler definition
  allocate(mapfits(0:npix-1,1:nmap))

  !HEALPix function
  !https://healpix.sourceforge.io/html/sub_read_bintab.htm
  !input arguments: #1 FITS-file name, #3 total pixel size
  !output arguments: #2 map handler, #4 value of missing pix, #5 if missing pix
  call read_bintab(mapname,mapfits,npix,nmap,nullval,anynull)

  !map takes over mapfits
  map=mapfits(:,1)
  !delete mapfits
  deallocate(mapfits)
  write(*,*)'info: target map read'

  !==========================================================================
  !  Importation of the galactic mask, ring
  !========================================================================== 

  allocate(mapmask(0:npix-1,1:nmap))

  !HEALPix function
  call read_bintab(maskname,mapmask,npix,nmap,nullval,anynull)

  !set masked area by HEALPix.bad_value
  where(mapmask(:,1)==0.)
    map=badvalue
  end where

  umpix=count(mapmask(:,1) > 0) !unmasked pixels
  write(*,*)'info: mask map applied'
  deallocate(mapmask)

  !==========================================================================
  !       Properties and normalisation of the maps (excluding masked pixels)
  !========================================================================== 

  mean=sum(map(0:npix-1),mask=map>badvalue)/umpix
  write(*,*)'info: partial sky mean =',mean
  !where (map>badvalue) map=map-mean !remove the monopole of partial sky
  !mean=sum(map(0:npix-1),mask=map>badvalue)/umpix
  !write(*,*)'info: regulated partial sky mean =',mean

  variance=0d0
  minval=map(0)
  maxval=map(0)

  do i=0,npix-1
    if (map(i)>badvalue) then !un-masked area
      variance=variance+(map(i)-mean)**2
      !update min and max values
      if (map(i)<minval) minval=map(i)
      if (map(i)>maxval) maxval=map(i)
    endif
  enddo
  write(*,*)'info: partial sky min =',minval,', max =',maxval

  variance=variance/(umpix-1)
  write(*,*)'info: partial sky var =',variance

  ecartt=sqrt(variance)
  write(*,*)'info: partial skay std =',ecartt

  !normalize map by its std
  !where (map>badvalue) map=map/ecartt
  !meanmap=mean/ecartt
  !minval=minval/ecartt
  !maxval=maxval/ecartt
  !write(*,*)'info: target map normalized'

  !==========================================================================
  !            Minkowski functionals of the map
  !========================================================================== 
   
  nameout='mf.dat'

  !call the subroutine defined in the following
  call mink(map,nameout,p,ecartt,maxnu,minval,mean,umpix,nside,npix)
  
  deallocate(map)

end program main


subroutine mink(dt,nameout,p,ecartt,maxnu,minval,meanmap,umpix,nside,npix)

  !Minkowski functionals of the map. Loop on threshold values, functionals are
  !calculated at each threshold and written on an ASCII file.
  !arguments:
  !input
  !#1 dt, target map (masked and normalized)
  !#2 nameout, output plain file path
  !#3 p, number of thresholds
  !#4 ecartt, target map std (partial sky)
  !#5 maxnu, threshold upper limit, in std unit
  !#6 minval, target map min (parital sky)
  !#7 meanmap, target map mean value (normalized)
  !#8 umpix, un-masked pixel number

  USE healpix_types
  USE fitstools
  USE CND_REG2D_mod
  USE pix_tools
 
  implicit none
  integer(I4B),intent(in)          :: p, umpix
  integer(I4B),intent(in)          :: nside, npix
  real(DP), dimension(0:npix-1)    :: dt
  real(DP)                         :: level,ecartt,maxnu,meanmap,minval
  integer(I4B)                     :: j,k
  !TYPE EXCRS_STST defined in CND_REG2D_mod.f90
  type(EXCRS_STAT)                 :: stats
  integer(I4B)                     :: nvoids,ordering
  integer(I4B)                     :: result_volume,result_genus
  integer(I4B)                     :: result_sides, result_maskedsides
  real(DP)                         :: result_length
  real(DP)                         :: v0,v1,v2,fsky
  character*80,intent(in)          :: nameout
  !TYPE CND_CNTRL_TYPE defined in CND_REG2D_mod.f90
  TYPE(CND_CNTRL_TYPE)             :: CND_CNTRL

  !=========================================================================
  !            Initialisation
  !=========================================================================

  write(*,*)'info: in subroutine mink'

  ordering =  CND_CNTRL%RING !  ordering =  CND_CNTRL%NESTED
  CND_CNTRL%CONNECT = CND_CNTRL%STAT !  CND_CNTRL%CONNECT = CND_CNTRL%EXACT 

  fsky=umpix
  fsky=fsky/npix 

  open(1,file=nameout)
  write(1,*)'# info: sky fraction =',fsky
  write(1,*)'# info: initial RMS of the map =',ecartt
  write(*,*)'info: mink initialized'

  !=========================================================================
  !            Loop on thresholds
  !========================================================================= 

  stats%do_sides =.false.

  !loop through threshold levels
  do j=0,p-1
    ! level in signal units, maxnu is the maximum sigma span of std
    level=meanmap-maxnu*ecartt+j*2.d0*maxnu*ecartt/(p-1)
    if ( CND_CNTRL%CONNECT == CND_CNTRL%STAT ) then
      if ( level < meanmap ) then
        CND_CNTRL%EXCURSION = CND_CNTRL%BELOW
      else
        CND_CNTRL%EXCURSION = CND_CNTRL%ABOVE
      endif
    endif

    !call CND_REG subroutine defined in CND_REG2D_mod.f90
    call CND_REG(dt,nside,ordering,level,stats,nvoids,CND_CNTRL)

    result_genus=0
    result_sides=0
    result_maskedsides=0
    result_volume=0
    result_length=0.d0

    do k=1,nvoids
      result_volume =result_volume+stats%n(k)
      if (stats%do_genus) result_genus = result_genus+stats%g(k)
      if (stats%do_sides) then
        result_sides = result_sides+stats%s(k)
        result_maskedsides = result_maskedsides+stats%sm(k)
      endif
      if (stats%do_length)result_length= result_length+stats%sf(k)
    enddo

    if ( CND_CNTRL%CONNECT == CND_CNTRL%STAT ) then
      if ( level < meanmap ) then
        result_volume = npix - result_volume
        result_genus  = 8 - result_genus
      endif
    endif

    !v0: first Minkowski funtional, volume (area)

    if ( level<meanmap ) then 
      v0=result_volume-(npix-umpix)
    else 
      v0=result_volume 
    endif
    v0=v0/umpix

    !v1: second Minkowski funtional, perimeter
    v1=result_length/(4.d0*Pi*fsky)/4.d0

    !v2: third Minkowski funtional, genus
    v2=result_genus
    v2=v2/(4.d0*Pi*fsky)/4.d0

    write(1,'(e10.3,2X,e10.3,2X,e10.3,2X,e10.3,2X,i7)')level,v0,v1,v2,nvoids

    if ( associated(stats%n) ) deallocate(stats%n)  ! clear stats arrays
    !if ( associated(stats%s) ) deallocate(stats%s)
    !if ( associated(stats%sm) ) deallocate(stats%sm)
    if ( associated(stats%sf)) deallocate(stats%sf)
    if ( associated(stats%g) ) deallocate(stats%g)

  enddo
  close(1)

  ! Get the properties of the manifold. Can be used for correction

  stats%do_genus=.true.
  stats%do_sides=.true.
  stats%do_length=.false.

  CND_CNTRL%EXCURSION = CND_CNTRL%ABOVE 

  call CND_REG(dt,nside,ordering,minval-ecartt,stats,nvoids,CND_CNTRL)

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

  write(*,*) 'Properties of the manifold',result_genus/4.,result_volume,result_sides,result_maskedsides,nvoids

  if ( associated(stats%n) ) deallocate(stats%n)  ! clear stats arrays
  if ( associated(stats%s) ) deallocate(stats%s)
  if ( associated(stats%sm) ) deallocate(stats%sm)
  !if ( associated(stats%sf)) deallocate(stats%sf)
  if ( associated(stats%g) ) deallocate(stats%g)
  
end subroutine mink
