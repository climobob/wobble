 program readnc_for_sbaam
!
! Old Sigma and NEMSIS files save vertical levels from low (higher pressure) to high (lower pressure)
! The new NetCDF files save vertical levels from high (lower pressure) to low (higher pressure)
! need to reverse the heigh order before calling the SBAAM subroutine
!

 use netcdf
 implicit none

 integer, parameter :: n2d=2, n3d=14
 integer, parameter :: MAX_ATT_LEN = 80
 character (len = *), parameter :: UNITS = "units"
 character(len=300) :: input_file
 character*200      :: argument    
 character(len=7)   :: var2d(n2d)
 character(len=7)   :: var3d(n3d)
 real(4), allocatable  :: pressfc(:,:)
 real(4), allocatable  :: dummy3d(:,:,:)
 real(4), allocatable  :: u(:,:,:), v(:,:,:)
 real(4)               :: dtor, rad, sum
 real(4)               :: time
 real(4), allocatable  ::  ak(:), bk(:)
 real(4), allocatable  ::  sl(:), lat(:), clat(:), slat(:), wlat(:)
 character*(MAX_ATT_LEN) :: time_units

 data var2d/'hgtsfc ','pressfc'/
 data var3d/ 'tmp    ', 'ugrd   ', 'vgrd   ', 'dpres  ', 'dzdt   ', &
 &           'delz   ', 'o3mr   ', 'spfh   ', 'cld_amt', 'clwmr  ', &
 &           'rwmr   ', 'snmr   ', 'icmr   ', 'grle   '/

 integer  :: error, ncid, id_var,dimid
 integer  :: xtype, len, attnum
 integer  :: im,jm,km, n, i, j, k,nargs, itime
 integer, dimension(5) :: mdate, ndate
!real(4)  :: sf, ofs
 character(len=10) :: dim_nam   

 nargs = iargc()   
 if (nargs<1) then
  print*, "use: readnc_for_sbaam file.nc"
  stop
 endif 
 call getarg(1,argument)   
 input_file=trim(argument)
 print*,'input_file=',input_file

 error=nf90_open(input_file,nf90_nowrite,ncid)
 error=nf90_inq_dimid(ncid,"grid_xt",dimid)
 error=nf90_inquire_dimension(ncid,dimid,dim_nam,im)
 error=nf90_inq_dimid(ncid,"grid_yt",dimid)
 error=nf90_inquire_dimension(ncid,dimid,dim_nam,jm)
 error=nf90_inq_dimid(ncid,"pfull",dimid)
 error=nf90_inquire_dimension(ncid,dimid,dim_nam,km)
 print*, "im,jm,km:",im,jm,km  

 error=nf90_inq_varid(ncid, 'time', id_var)
!print *, ' nf90_inq_varid error = ', error
 error=nf90_get_var(ncid, id_var, time)
!print *, ' nf90_get_var error = ', error
 print*,'time: ',time
 error=nf90_get_att(ncid, id_var, UNITS, time_units)
 time_units=trim(time_units)
 print *, ' time_units = ', time_units
 read(time_units,100) mdate
 100 format(12x,i4,1x,i2,1x,i2,1x,i2,1x,i2)
 print *, ' mdate = ', mdate
 itime = time
 call date_offset(mdate,itime*100,ndate)
 print *, ' ndate = ', ndate
 

 error=nf90_inquire_attribute(ncid, NF90_GLOBAL, "ak", xtype, len, attnum)
 print *, ' ak xtype = ', xtype,'   len = ', len, '   attnum = ', attnum
 allocate (ak(len))
 print *, 'Is AK allocated ? ', allocated(ak)
 error=nf90_get_att(ncid, NF90_GLOBAL, 'ak', ak)
 error=nf90_inquire_attribute(ncid, NF90_GLOBAL, "bk", xtype, len, attnum)
 print *, ' bk xtype = ', xtype,'   len = ', len, '   attnum = ', attnum
 allocate (bk(len))
 print *, 'Is BK allocated ? ', allocated(bk)
 error=nf90_get_att(ncid, NF90_GLOBAL, 'bk', bk)
!print *, ' ak = ', ak
!print *, ' bk = ', bk
 allocate (sl(len))
 print *, 'Is SL allocated ? ', allocated(sl)
 do k = 1, len
    sl(len-k+1) = bk(k) + ak(k) / 100000.0
 enddo
 print *, ' sl = ', sl

 allocate (lat(jm))
 print *, 'Is LAT allocated ? ', allocated(lat)
 allocate (clat(jm))
 print *, 'Is CLAT allocated ? ', allocated(clat)
 allocate (slat(jm))
 print *, 'Is SLAT allocated ? ', allocated(slat)
 allocate (wlat(jm))
 print *, 'Is WLAT allocated ? ', allocated(wlat)
 error=nf90_inq_varid(ncid, 'grid_yt', id_var)
 error=nf90_get_var(ncid, id_var, lat)

!---get latitude weights
!dtor = 3.1415927/180.0
!do j = 1, jm
!   rad = lat(j)*dtor
!   slat = sin(rad)
!   clat = cos(rad)
!enddo

 call splat(4,jm,slat,wlat)
 clat = sqrt(1.0 - slat*slat)

!sum = 0.0
!do i = 1, jm
!   print *, 'lat_index ',i,' =', lat(i), slat(i), clat(i), wlat(i)
!   sum = sum + wlat(i)
!enddo
!print *,' Sum of wlat = ', sum

! -- 2d variables
 allocate (pressfc(im,jm))
 print *, 'Is PRESSFC allocated ? ', allocated(pressfc)
 error=nf90_inq_varid(ncid, 'pressfc', id_var)
 error=nf90_get_var(ncid, id_var, pressfc)
 print*,'pressfc',' max min:',maxval(pressfc),minval(pressfc)

! -- 3d variables
 allocate (dummy3d(im,jm,km))
 print *, 'Is PRESSFC allocated ? ', allocated(pressfc)
 allocate (u(im,jm,km))
 print *, 'Is UGRD allocated ? ', allocated(u)
 allocate (v(im,jm,km))
 print *, 'Is VGRD allocated ? ', allocated(v)
 error=nf90_inq_varid(ncid, 'ugrd', id_var)
 error=nf90_get_var(ncid, id_var, dummy3d)
 do k = 1, km
   u(:,:,km-k+1) = dummy3d(:,:,k)
 enddo
!do k = 1, km
!  print *,'dummy',' max min:',maxval(dummy3d(:,:,k)),minval(dummy3d(:,:,k))
!enddo
!do k = 1, km
!  print *,'u',' max min:',maxval(u(:,:,k)),minval(u(:,:,k))
!enddo

 error=nf90_inq_varid(ncid, 'vgrd', id_var)
 error=nf90_get_var(ncid, id_var, dummy3d)
 do k = 1, km
   v(:,:,km-k+1) = dummy3d(:,:,k)
 enddo
!do k = 1, km
!  print *,'dummy',' max min:',maxval(dummy3d(:,:,k)),minval(dummy3d(:,:,k))
!enddo
!do k = 1, km
!  print *,'v',' max min:',maxval(v(:,:,k)),minval(v(:,:,k))
!enddo

 deallocate(dummy3d, ak, bk)

 open(52,form='formatted',status='unknown',position='append')
 write(52,200) mod(mdate(1),100), mdate(2), mdate(3), mdate(4), itime, mod(ndate(1),100), ndate(2), ndate(3), ndate(4)
 200 format(1X,4I4,1X,I3,1X,4I2.2)
 call energy(pressfc,u,v,im,jm,sl,wlat,slat,clat,km)

 deallocate(sl)
 deallocate(lat, clat, slat, wlat)
 deallocate(u, v, pressfc)

end
