!=====without external file,,data is inserted manually=====
!______________________________________________
!!        implicit none
!!        integer :: n,nn
!!        real :: Ecm(20), sig (20), elint, eN
!!        data (Ecm(nn),nn=1,4)/0.4, 0.6, 0.8, 1.0/ 
!!        data (sig(nn),nn=1,4)/1, 2, 3, 4/
!!        do n = 1, 4 
!!        elint = ((Ecm (n)/10)**2) 
!!        eN = elint*sig(n)
!!        print*, eN, sig(n), elint*10**(-5)
!!        end do
!_____________________________________________   
!----------------------------------------------

!======== with external data file x,y===========
        program scale_change
        implicit none
        integer::n
        real :: Ecm(20), sig(20), elint, eN
        open(unit=10, file='ttbar.dat')
        write(*,*) "eN,sig,elint"
        do n = 1, 20
        read (10,*) Ecm (n), sig(n)
        end do
        close (10)
        do n = 1, 20
        elint = ((Ecm (n)/10)**2)*10**(-5.0)
        eN=elint*sig(n)
        print*,Ecm(n),sig(n),eN,elint 
        end do 
        end! program scale_change   
