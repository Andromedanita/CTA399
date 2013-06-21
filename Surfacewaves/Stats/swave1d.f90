! requires -r8 because of the small search step size
integer, parameter :: nx=1024*4,iscale=8,ny=iscale*nx
real, dimension(ny) :: dx,r1,r2,rhot
complex, dimension(ny) :: cdx
!real, dimension(nx,ny) :: rho
real, dimension(nx) :: rho1,dz

pi=4*atan(1.)
call random_seed()
call random_number(r1)
call random_number(r2)
dx=sqrt(-2*log(r1))*cos(2*pi*r2)

write(*,*)'f1',ny,sum(dx**2)/ny
cdx=dx/ny
call four1(cdx,ny,1)
do j=1,ny/2+1
   ak=j-1
   if (j .eq. 1) then
	cdx(j)=0
	cycle
   endif
   cdx(j)=abs(ak)**(-0.5)*exp(-(ak)**2/(ny/320./2)**2/2)*cdx(j)
enddo
cdx(ny/2+2:)=conjg(cdx(ny/2:2:-1))
call four1(cdx,ny,-1)
dx=real(cdx)
sdx=sqrt(sum(dx**2)/ny)
write(*,*)'f2',sdx
dx=dx*20*2/sdx

rho1=0
sigma=2
open(10,file='sheet.dat')
do j=1,ny
   x0=j/iscale+dx(j)
   write(10,*) x0
   i0=max(1.,x0-4*sigma)
   i1=min(nx*1.,x0+4*sigma)
   do i=i0,i1
      rho1(i)=rho1(i)+exp(-(i-x0)**2/2/sigma**2)/sigma
   end do
enddo
rho1=rho1*nx/sum(rho1)
open(10,file='rho.dat')
do i=1,nx
   write(10,*) rho1(i)
enddo
open(10,file='dt.dat')
rho1=rho1*1000!*10
dz=cshift(rho1,1)-rho1
do i=1,nx-1
   write(10,*) i+dz(i)
enddo
rhot=0
open(10,file='test.dat')
do xi=1,ny,0.001
   i=xi
   ii=xi/iscale+1
   if (ii+1>nx) cycle
   y=xi/iscale+1
   jm=y
   jp=jm+1
   w=jp-y
   dz1=w*dz(jm)+(1-w)*dz(jp)
   jj=nint(dz1+y)
   if (jj .ne. nx/2) cycle
   ak2=dz(jp)-dz(jm)
   rhot(i)=1/(1+ak2)
!   write(10,*) xi,y,rhot(i)
end do
open(10,file='flux.dat')
do i=1,ny,iscale
   write(10,*) sum(rhot(i:i+iscale-1))/(count(rhot(i:i+iscale-1) .ne. 0)+1.e-10)
end do
contains

  subroutine pmap(fn,rmap1,nx,ny,iscale)
  real rmap(nx,ny),rmap1(nx,ny)
  integer*2, dimension(nx,ny) :: imap
  integer*1, dimension(nx,ny) :: imap1
  character(len=*):: fn
  integer npix,mypos

  npix=min(ny/2-1,nx/2-1,300)
  
  
  rmap=rmap1
!  write(*,*) 'rms=',sqrt(sum(rmap(nx/2-npix:nx/2+npix,ny/2-npix:ny/2+npix)**2)/npix**2/4)
  do while (iscale > 1)      
     rmap=sign((sqrt(abs(rmap))),rmap)
     iscale=iscale-1
  end do
  rmax=maxval(rmap)
  rmin=minval(rmap)
  write(*,*) trim(fn),rmax,rmin
  imap=255*(rmap-rmin)/(rmax-rmin)
  imap1=127*(rmap-rmin)/(rmax-rmin)
  open(10,file=fn)
  write(10,'(2hP5)')
  write(10,*)nx,ny
  write(10,*) 255
!  write(10,*) 127
  INQUIRE(UNIT=10, POS=mypos)
  close(10)
  open(10,file=fn, access='stream',position='append')
!  write(10,pos=mypos) int(imap,1)
  write(10) int(imap,1)
  close(10)
end subroutine pmap


end program
