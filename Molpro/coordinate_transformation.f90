        module para
        implicit none
        real*8,parameter::pi=dacos(-1d0)
        real*8::r1,r2,r3,th1,th2,th3
        real*8::xn,yn,zn,xh1,yh1,zh1,xh2,yh2,zh2,xh3,yh3,zh3
        real*8::h1h2,h2h3,h3h1,R
        real*8::alpha1,alpha2,alpha3,h1m1,h2m2,h3m3,nm3
        real*8::h1h3m3,h2h3m3,nh2h1,nh3m3
        integer*8::i
        contains
         subroutine angle_and_length(r1,r2,r3,th1,th2,th3,h1h2,h2h3,h3h1,alpha1,alpha2,alpha3)
         implicit none
         real*8::h1h2,h2h3,h3h1,r1,r2,r3,th1,th2,th3
         real*8::alpha1,alpha2,alpha3
         h1h2=dsqrt(r1**2+r2**2-2*r1*r2*dcos(pi*th1/180d0))
         h2h3=dsqrt(r2**2+r3**2-2*r2*r3*dcos(pi*th2/180d0))
         h3h1=dsqrt(r3**2+r1**2-2*r3*r1*dcos(pi*th3/180d0))
         alpha1=(dacos((h1h2**2+h3h1**2-h2h3**2)/(2*h1h2*h3h1)))*180/pi
         alpha2=(dacos((h1h2**2+h2h3**2-h3h1**2)/(2*h1h2*h3h1)))*180/pi
         alpha3=(dacos((h2h3**2+h3h1**2-h1h2**2)/(2*h2h3*h3h1)))*180/pi
         !write(*,*)alpha1,alpha2,alpha3,h1h2,h2h3,h3h1
         end subroutine

         !median_line1:to get the coordinates of H1 H2 and H3
         subroutine median_line1(h1h2,h2h3,h3h1,alpha1,alpha2,alpha3,h1m1,h2m2,h3m3,h2h3m3,h1h3m3)
         implicit none
         real*8::h1h2,h2h3,h3h1,alpha1,alpha2,alpha3
         real*8::h1m1,h2m2,h3m3,h2h3m3,h1h3m3
         h3m3=dsqrt((h1h2/2d0)**2+h2h3**2-2*(h1h2/2d0)*h2h3*dcos(pi*alpha2/180))
         h2h3m3=(dacos((h3m3**2+h2h3**2-(h1h2/2d0)**2)/(2*h3m3*h2h3)))*180/pi
         h1h3m3=alpha3-h2h3m3
         !write(*,*)h3m3,h2h3m3,h1h3m3
         end subroutine
                
         !median_line2:to get nm3 and nh3m3 and R
         subroutine median_line2(r1,r2,h1h2,nh2h1,nm3,h3m3,nh3m3,R)
         implicit none
         real*8::nh2h1,h1h2,nm3,nh3m3,h3m3
         real*8::r1,r2,R
         nh2h1=(dacos((h1h2**2+r2**2-r1**2)/(2*h1h2*r2)))*180/pi
         nm3=dsqrt((h1h2/2d0)**2+r2**2-2*(h1h2/2d0)*r2*dcos(pi*nh2h1/180))
         nh3m3=(dacos((r3**2+h3m3**2-nm3**2)/(2*r3*h3m3)))*180/pi
         R=dsqrt(r3**2+((2d0/3d0)*h3m3)**2-2*r3*(2d0/3d0)*h3m3*dcos(pi*nh3m3/180))
         !write(*,*)nh2h1,nm3,nh3m3,h3m3,R
         end subroutine
        end module
        
        program main
        use para
        implicit none
        !real*8::r1,r2,r3,th1,th2,th3
        !real*8::xn,yn,zn,xh1,yh1,zh1,xh2,yh2,zh2,xh3,yh3,zh3
        !real*8::h1h2,h2h3,h3h1,R
        !real*8::alpha1,alpha2,alpha3,h1m1,h2m2,h3m3,nm3
        !real*8::h1h3m3,h2h3m3,nh2h1,nh3m3
        !integer::i
        open(1,file='NH3.dat')
        open(2,file='xyz_NH3.dat')
        do i=1,35
        read(1,*)r1,r2,r3,th1,th2,th3
        call angle_and_length(r1,r2,r3,th1,th2,th3,h1h2,h2h3,h3h1,alpha1,alpha2,alpha3)
        call median_line1(h1h2,h2h3,h3h1,alpha1,alpha2,alpha3,h1m1,h2m2,h3m3,h2h3m3,h1h3m3)
        !write(*,*)h3m3,h1h2,h2h3,h3h1,h2h3m3,h1h3m3
        call median_line2(r1,r2,h1h2,nh2h1,nm3,h3m3,nh3m3,R)
        !write(*,*)h3m3,h2h3m3
        xh3=(2d0/3d0)*h3m3
        yh3=0d0
        zh3=0d0
        xh2=(2d0/3d0)*h3m3-h2h3*dcos(pi*h2h3m3/180d0)
        yh2=h2h3*dsin(pi*h2h3m3/180d0)
        zh2=0d0
        xh1=(2d0/3d0)*h3m3-h3h1*dcos(pi*h1h3m3/180d0)
        yh1=-h3h1*dsin(pi*h1h3m3/180d0)
        zh1=0d0
        xn=(R**2-r3**2+xh3**2)/(2*xh3)
        yn=(R**2-2*xh2*xn+xh2**2+yh2**2-r2**2)/(2*yh2)
        zn=dsqrt(abs(R**2-xn**2-yn**2))
        write(2,"(3f12.4)")r3
        write(2,"(3f14.8)")xh1,yh1,zh1
        write(2,"(3f14.8)")xh2,yh2,zh2
        write(2,"(3f14.8)")xh3,yh3,zh3
        write(2,"(3f14.8)")xn,yn,zn
        !write(*,*)R
        end do
        end program
