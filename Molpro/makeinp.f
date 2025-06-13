        implicit real*8(a-h,o-z)
        integer ik
        parameter (ndd=10000)
        parameter (ncn=70,nnn=70)
        character*80 line,ff1*15,ff2*16
        open(1,file='h2s.dat')
        open(222,file='cmfy.in')
        open(255,file='job.sh')
        do 100 ik=1,330
        read(1,*)theta,rhs1,rhs2
        call fname(ik,ff1,ff2)
        open(20,file=ff1)
        write(255,*)' molpro23 -n 8  -W ./ -d tmp/ ',ff1
      !  write(30,*)' rm -f tmp/',ff2
        rewind(222)
        read(222,1001)line
        write(20,1001)line
 1001   format(a)
        read(222,1001)line
        line(8:19)=ff2
        write(20,1001)line
        
        read(222,1001)line
        write(20,1001)line
        read(222,1001)line
        write(20,1001)line
        read(222,1001)line
        read(222,1001)line
        read(222,1001)line
        write(20,1001)line

        read(222,1001)line
        read(222,1001)line
        read(222,1001)line
        write(20,1002)'th=',theta
        write(20,1002)'rhs1=',rhs1,'  ang'
        write(20,1002)'rhs2=',rhs2,'  ang'
 1002   format(a,f10.5,a)
        do ii=1,120
        read(222,1001)line
        write(20,1001)line
        enddo
        close(20)
 100    continue
         
        end

        subroutine fname(ik,ff1,ff2)
        implicit real*8(a-h,o-z)
        character*5 ci,ff1*15,ff2*16
           rewind(777)
            write(777,21)ik
 21       format(i5)
             rewind(777)
            read(777,22)ci
 22      format(a)
            do 2 i=1,5
           if (ci(i:i).ne.' ') then
           kk=i
           goto 3
           endif
  2      continue
  3       continue
          ff1='h2s-'//ci(KK:5)//'.in'
          ff2='h2s-'//ci(KK:5)//'.wfu'
          end

