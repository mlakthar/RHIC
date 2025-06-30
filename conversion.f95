program rapdistrb_pseudorapdistrb
implicit none
integer::i
real,dimension(20)::M,N
real::beta,m0,E

character(14)::myfilename='rapspectra.dat'
Print*,"value of m0 and E"
read*, m0,E
!M stands for pseudo rap distrb ,E centre of mass energy
!N standas for rap distrb,m0 rest mass
open(2,file='rapspectra.dat',action='read')
read(2,*) (M(i),i=1,20)
close(2)
print*,(M(i),i=1,20)
beta=sqrt(1-(m0/(E*10**3))**2)
do i=1,20
N(i)=beta*M(i)
end do
print*,N
!!!!!!!!!!!!!!!!!!!!!
!data file for output 
open(1,file='pseudospectra.dat', action='write')
do i=1,20
  write(1,*)  M(i),N(i)
  
  end do
end program rapdistrb_pseudorapdistrb