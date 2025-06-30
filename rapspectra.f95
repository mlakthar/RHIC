!this program tell us spectra of rapidity of  particles in limited phase space.
program  rapidity_spectra
implicit none
!declaration of variable

real,dimension(20)::y,a,B,M
real ,parameter::pi=4.0*atan(1.0)
integer::i,g
real::T,alpha,beta,m0,cpot,h
!character(14)::myfilename="rapspectra.dat "
!M=dn/dy
!units are taken in GeV,GeV/c^2
!cpot :chenical potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Print*,"value of T,cpot,g,m0"
read*,T,cpot,g,m0
!g=4;m0=938;
h=0.005
beta=1/T
alpha=EXP(beta*cpot)!fugacity
!open(1,file='rapspectra.dat', action='write')
do i=1,20
y(i)=(i*h)
a(i)=beta*cosh(y(i))
B(i)=(1/((2*pi)**2))*g*alpha*cosh(y(i))
M(i)=B(i)*EXP(-a(i)*m0)*((m0**2/a(i))+(2*m0/a(i)**2)-(2/a(i)**3))
Print* ,"y=",y(i) ," M=",M(i)
!write(1,*)   y(i) , M(i)
!output data into a file
open(1,file='rapspectra.dat', action='write')
!Print*,    y(i)  ,         M(i)
!Print*, "numb=",i,  "y=",y(i), " M=",M(i)
write(1,*)   M(i)
end do
end program  rapidity_spectra
