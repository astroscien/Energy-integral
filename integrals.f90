# Energy-integral

module integrals
use typedefs
use background
use forcing
use solution
implicit none

contains

subroutine response_integrals(x,energy, work_Udrho, work_zahn)
implicit none
real (dp), intent (in) :: x(:)
real (dp), intent (out) :: energy, work_Udrho, work_zahn
type (response) :: r
double precision :: a,cfluid,cgrav_in,cgrav_out
double precision :: rim1,rip1,dh
integer :: i,nr

call unpack_response(x,r)

nr=bg%nr

a=0.d0
cfluid=0.d0
cgrav_in=0.d0
cgrav_out=0.d0
work_Udrho=0.d0

do i=1,nr

  if (i>1) then
    rim1=0.5*(bg%r(i)+bg%r(i-1))
  else
    rim1=0.d0
  endif
  if (i<nr) then
    rip1=0.5*(bg%r(i)+bg%r(i+1))
  else
    rip1=bg%r(nr)
  endif

  a=a + bg%rho(i) * (rip1**3-rim1**3)/3.0 &
  * ( r%xir(i)**2 + f%l * (f%l+1)*r%xih(i)**2 )

  dh=r%dp(i)/bg%rho(i)

  cfluid=cfluid + (rip1**3-rim1**3)/3.0 * bg%rho(i) &
  * ( dh**2/bg%csq(i) + bg%nsq(i)*r%xir(i)**2 )
  ! * dh**2/bg%csq(i)
  ! * bg%nsq(i)*r%xir(i)**2

  cgrav_in=cgrav_in - (rip1**3-rim1**3)/3.0 &
  * ( r%dg(i)**2 + f%l*(f%l+1)*(r%dphi(i)/bg%r(i))**2 )/(4.0*pi*cgrav)

  work_Udrho = work_Udrho - bg%rho(i)*(rip1**3-rim1**3)/3.0 &
  * f%U(i) * ( dh/bg%csq(i) + bg%nsq(i)/bg%g(i)*r%xir(i) )

  !print *,bg%r(i)/bg%r(bg%nr),f%omega**2 * a/bg%Eeq,cfluid/bg%Eeq,cgrav_in/bg%Eeq,cgrav_out/bg%Eeq

end do

cgrav_out= - (f%l+1)*r%dphi(nr)**2 * bg%r(nr)/(4.0*pi*cgrav)

!print *,1.d0,f%omega**2 * a/bg%Eeq,cfluid/bg%Eeq,cgrav_in/bg%Eeq,cgrav_out/bg%Eeq

energy = f%omega**2 * a + cfluid + cgrav_in + cgrav_out

! zahn (1970) nice trick for work integral. get by integration by parts.
! should agree with work_Udrho, i.e. this is a check.
work_zahn=(2*f%l+1)*f%U(nr)*r%dphi(nr) * bg%r(nr)/(4.0*pi*cgrav)

if (dbg) then
  print *
  print *,"subroutine: response_integrals"
  print *,"l=",f%l
  print *,"omega/omega0=",f%omega/bg%omega0
  print *,"energy/E0=",energy/bg%E0
  print *,"work_Udrho/E0=",work_Udrho/bg%E0
  print *,"work_zahn/E0=",work_zahn/bg%E0
  print *,"end subroutine response_integrals"
  print *
endif

end subroutine response_integrals



end module integrals
