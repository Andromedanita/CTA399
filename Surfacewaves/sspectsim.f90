integer, parameter :: n=1024,nf=4096*4
complex, dimension(n,n) :: cspect
real, dimension(n,n) :: rspect
real, dimension(nf) :: rphase
complex ctmp
integer nn(2)

nn=n

pi=4*atan(1.)
call random_number(rphase)
rphase=rphase*2*pi

cspect=0
open(10,file='flux.dat')
do i=1,nf
   read(10,*) amu
   if (amu .eq. 0) cycle
   x=(i-nf/2.)/(nf/2.)
   ix=n/4*x+1
   iy=n/2*x**2+1
   ctmp=exp((0.,1.)*rphase(i))*sqrt(abs(amu))
   cspect(ix,iy)=ctmp
enddo
rspect=cshift(abs(cspect),n/2)
rspect(:,:n/2)=rspect(:,n/2:1:-1)
call pmap('sspect0.pgm',rspect,n,n/2,1)
call fourn(cspect,nn,2,-1)
rspect=abs(cspect)**2
rspect=rspect*n**2/sum(rspect)
cspect=rspect
rspect=transpose(rspect)
call pmap('rspect.pgm',rspect,n,n,2)
call fourn(cspect,nn,2,1)
cspect(1,1)=0
rspect=cshift(abs(cspect),n/2)
rspect(:,:n/2)=rspect(:,n/2:1:-1)
call pmap('sspect.pgm',rspect,n,n/2,1)

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

end
