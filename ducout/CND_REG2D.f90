!Minkowski functionals for one temperature map

!v1, September 2012


!p=number of thresholds
!maxnu=maximum threshold (symetric)
!map='map_cmb_ns64.fits'
!mask='mask_gal_fsky0.80_ns64.fits'

!results in text file "mf.dat"


program main

  USE healpix_types
  USE fitstools
  USE CND_REG2D_mod


  implicit none
  integer(I4B)                       :: p
  integer(I4B),parameter             :: nmap=1
  integer(I4B),parameter             :: npix=49152
  integer(I4B)                       :: i,umpix
  real(DP)                           :: mean,variance,minval,maxval  
  real(DP)                           :: ecartt,nullval
  real(DP)                           :: galvalue
  real(DP)                           :: maxnu,meanmap
  real(DP),dimension(:),allocatable  :: map
  real(DP),dimension(:,:),allocatable:: mapmask,mapfits
  character*80                       :: nameout  
  character*80                       :: mapname,maskname
  logical                            :: anynull





!===========================================================================
!            Parameters
!=========================================================================== 

mapname='map_cmb_ns64.fits'
maskname='mask_gal_fsky0.80_ns64.fits'

p=26 !p=number of thresholds
maxnu=3.5d0 !in units of sigma(map)


galvalue=-1.6375000d30 !healpix.bad_value

!===========================================================================
!            Reading map
!=========================================================================== 


allocate(map(0:npix-1))

allocate(mapfits(0:npix-1,1:nmap))
  
call read_bintab(mapname,mapfits,npix,nmap,nullval,anynull)
      
map=mapfits(:,1)
      
deallocate(mapfits)
      

!==========================================================================
!  Importation of the galactic mask, ring, nside=64
!========================================================================== 
   
allocate(mapmask(0:npix-1,1:nmap))
   
call read_bintab(maskname,mapmask,npix,nmap,nullval,anynull)
   
where(mapmask(:,1)==0.)
   map=galvalue
end where

umpix=count(mapmask(:,1) > 0) !unmasked pixels
write(*,*)'umpix=',umpix
  
deallocate(mapmask)


!==========================================================================
!       Properties and normalisation of the maps (excluding masked pixels)
!========================================================================== 
   
mean=sum(map(0:npix-1),mask = map > -1d15)/umpix

where (map>-1d15) map=map-mean
mean=sum(map(0:npix-1),mask = map > -1d15)/umpix
  
variance=0d0
minval=10d0
maxval=-10d0

do i=0,npix-1
   if (map(i) > -1d15) then
      variance=variance+(map(i)-mean)**2
      if (map(i) < minval) minval=map(i)
      if (map(i) > maxval) maxval=map(i)
   endif
enddo

variance=variance/(umpix-1)

ecartt=sqrt(variance)
write(*,*)'standard deviation=',ecartt
    
!    Normalisation of the map

where (map>-1d15) map=map/ecartt
meanmap=mean/ecartt
minval=minval/ecartt
maxval=maxval/ecartt


!==========================================================================
!            Minkowski functionals of the map
!========================================================================== 
   
nameout='mf.dat'
  
call mink(map,nameout,p,ecartt,maxnu,minval,meanmap,umpix)
  
deallocate(map)



end program main
 






subroutine mink(dt,nameout,p,ecartt,maxnu,minval,meanmap,umpix)

  !Minkowski functionals of the map. Loop on threshold values, functionals are
  ! calculated at each threshold and written on an ASCII file.
  

  USE healpix_types
  USE fitstools
  USE CND_REG2D_mod
  USE pix_tools
  
  implicit none
  integer(I4B),intent(in)          :: p,umpix 
  integer(I4B),parameter           :: npix=49152,nside=64
  real(DP), dimension(0:npix-1)    :: dt
  real(DP)                         :: level,ecartt,maxnu,meanmap,minval
  integer(I4B)                     :: j,k
  type(EXCRS_STAT)                 :: stats
  integer(I4B)                     :: nvoids,ordering
  integer(I4B)                     :: result_volume,result_genus  
  integer(I4B)                     :: result_sides, result_maskedsides  
  real(DP)                         :: result_length
  real(DP)                         :: v0,v1,v2,fsky
  character*80,intent(in)          :: nameout
  TYPE(CND_CNTRL_TYPE)             :: CND_CNTRL
  

!===========================================================================
!            Initialisation
!=========================================================================== 

  ordering =  CND_CNTRL%RING !  ordering =  CND_CNTRL%NESTED
  CND_CNTRL%CONNECT = CND_CNTRL%STAT !  CND_CNTRL%CONNECT = CND_CNTRL%EXACT 
  

  fsky=umpix
  fsky=fsky/npix 
  write(*,*)'fsky=',fsky
  
  
  open(1,file=nameout)
  write(1,*)'initial rms of the map:',ecartt


!===========================================================================
!            Loop on thresholds
!=========================================================================== 

  stats%do_sides =.false.

  do j=0,p-1

     level=-maxnu+j*2.d0*maxnu/(p-1)  !map normalized by ecartt

     if ( CND_CNTRL%CONNECT == CND_CNTRL%STAT ) then
        if ( level < meanmap ) then 
           CND_CNTRL%EXCURSION = CND_CNTRL%BELOW 
        else   
           CND_CNTRL%EXCURSION = CND_CNTRL%ABOVE 
        endif
     endif

     
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
     
     
     write(1,'(f5.2,f13.5,f11.6,f11.7,i7)')level,v2,v1,v0,nvoids
     
     
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





end subroutine mink  


    
