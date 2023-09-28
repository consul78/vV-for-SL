        program vVmodes
c    this program computes the van Vleck modes starting from TM-O lengths
c    it first computes the average length of the TM-O bond
c    then it computes the J-T modes using the average as a reference
c
        implicit none
c
        integer ind(18), j, k
        double precision x1p(18),y1p(18),z1p(18)
        double precision x2p(18),y2p(18),z2p(18)
        double precision lxp(18),lyp(18),lzp(18)
        double precision dxp(18),dyp(18),dzp(18)
        double precision q1p(18),q2p(18),q3p(18)
        double precision x1m(18),y1m(18),z1m(18)
        double precision x2m(18),y2m(18),z2m(18)
        double precision lxm(18),lym(18),lzm(18)
        double precision dxm(18),dym(18),dzm(18)
        double precision q1m(18),q2m(18),q3m(18)
        double precision sumtmp,avgtmp,sumtmm,avgtmm
        character(len = 10) hash, comment1, comment2
        character(len = 10) com1, com2, com3, com4, com5, com6
        character(len = 10) w1, w2, w3, w4, w5, w6
c
        open(unit=11,file='TM_xx.dat',status='old')
        open(unit=21,file='vV123_modes.dat')
        open(unit=22,file='vV231_modes.dat')
        open(unit=23,file='vV312_modes.dat')
c
c    initialisation of the variables
c
        do j=1,18
           ind(j) = 0 ! the layer index
c
           lxp(j) = 0 !
           lyp(j) = 0 ! the O-Mn-O length - chain 'p'
           lzp(j) = 0 !
c
           dxp(j) = 0 !
           dyp(j) = 0 ! the difference of said length wrt average
           dzp(j) = 0 !
c
c          same as above - chain 'm'
c
           lxm(j) = 0
           lym(j) = 0
           lzm(j) = 0
           dxm(j) = 0
           dym(j) = 0
           dzm(j) = 0
c
c          the coordinates as read from TM_xx.dat
c
           x1p(j) = 0
           y1p(j) = 0
           z1p(j) = 0
           x2p(j) = 0
           y2p(j) = 0
           z2p(j) = 0
c
c          an attempt to read the coordinates recursively
c
c           do k=1,2
c              xp(j,k) = 0 ! coordinates - chain 'p' & 'm'
c              yp(j,k) = 0
c              zp(j,k) = 0
c              xm(j,k) = 0
c              ym(j,k) = 0
c              zm(j,k) = 0
c           end do
c
c          the output variables
c
           q1p(j) = 0 ! breathing on chain 'p'
           q2p(j) = 0 ! JT Q2 on chain 'p'
           q3p(j) = 0 ! JT Q3 on chain 'p'
           q1m(j) = 0 ! breathing on chain 'm'
           q2m(j) = 0 ! JT Q2 on chain 'm'
           q3m(j) = 0 ! JT Q2 on chain 'm'
        end do
c
c       useful sums and averages
c
        sumtmp = 0
        avgtmp = 0
        sumtmm = 0
        avgtmm = 0
c
c    computation of the TM-O bond length in all crystal directions
c
        read(11,*) hash
c        write(6,*) hash
        read(11,*) hash 
c        write(6,*) hash
        read(11,*) hash
c        write(6,*) hash
        read(11,*) hash
c        write(6,*) hash
        write(6,*) 'computing van Vleck modes ...'
        do j=1,18
           read(11,*) ind(j),x1p(j),y1p(j),z1p(j),x2p(j),y2p(j),z2p(j),
     &x1m(j),y1m(j),z1m(j),x2m(j),y2m(j),z2m(j)
c           write(6,*) "j=",j,"j+1 = ",j+1,"ind = ",ind(j)
           lxp(j) = (x1p(j) + x2p(j))
           lyp(j) = (y1p(j) + y2p(j))
           lzp(j) = (z1p(j) + z2p(j))
           lxm(j) = (x1m(j) + x2m(j))
           lym(j) = (y1m(j) + y2m(j))
           lzm(j) = (z1m(j) + z2m(j))
c           sumtm = sumtm + (lx(j)*ly(j)*lz(j))
        end do
c
c    computation of the average TM-O bond length
c    line 0 and n are the same
c
        do j=1,18
           sumtmp = sumtmp + ((lxp(j)*lyp(j)*lzp(j))**(1.d0/3.d0))
           sumtmm = sumtmm + ((lxm(j)*lym(j)*lzm(j))**(1.d0/3.d0))
        end do
c
        avgtmp = sumtmp/18
        avgtmm = sumtmm/18
c
c    computation of the van Vleck modes
c
        write(6,*) 'writing on disk ... 123'
        do j=1,18
           dxp(j) = lxp(j) - avgtmp
           dyp(j) = lyp(j) - avgtmp
           dzp(j) = lzp(j) - avgtmp
           dxm(j) = lxm(j) - avgtmm
           dym(j) = lym(j) - avgtmm
           dzm(j) = lzm(j) - avgtmm
           q1p(j) = ( dxp(j) + dyp(j) + dzp(j))/sqrt(3.d0)
           q2p(j) = ( dxp(j) - dyp(j))/sqrt(2.d0)
           q3p(j) = (-dxp(j) - dyp(j) + 2*dzp(j))/sqrt(6.d0)
           q1m(j) = ( dxm(j) + dym(j) + dzm(j))/sqrt(3.d0)
           q2m(j) = ( dxm(j) - dym(j))/sqrt(2.d0)
           q3m(j) = (-dxm(j) - dym(j) + 2*dzm(j))/sqrt(6.d0)
           write(21,101) ind(j),q1p(j),q2p(j),q3p(j),
     &q1m(j),q2m(j),q3m(j)
        end do
        write(21,101) 18,q1p(1),q2p(1),q3p(1),
     &q1m(1),q2m(1),q3m(1)
cc
cc    re-initialisation of the van Vleck modes
cc
c        do j=1,18
c           dxp(j) = 0
c           dyp(j) = 0
c           dzp(j) = 0
c           dxm(j) = 0
c           dym(j) = 0
c           dzm(j) = 0
c           q1p(j) = 0
c           q2p(j) = 0
c           q3p(j) = 0
c           q1m(j) = 0
c           q2m(j) = 0
c           q3m(j) = 0
c        end do
cc
cc    re-computation of the van Vleck modes
cc
c        write(6,*) 'writing on disk ... 231'
c        do j=1,18
c           dxp(j) = lxp(j) - avgtmp
c           dyp(j) = lyp(j) - avgtmp
c           dzp(j) = lzp(j) - avgtmp
c           dxm(j) = lxm(j) - avgtmm
c           dym(j) = lym(j) - avgtmm
c           dzm(j) = lzm(j) - avgtmm
c           q1p(j) = ( dyp(j) + dzp(j) + dxp(j))/sqrt(3.d0)
c           q2p(j) = ( dyp(j) - dzp(j))/sqrt(2.d0)
c           q3p(j) = (-dyp(j) - dzp(j) + 2*dxp(j))/sqrt(6.d0)
c           q1m(j) = ( dym(j) + dzm(j) + dxm(j))/sqrt(3.d0)
c           q2m(j) = ( dym(j) - dzm(j))/sqrt(2.d0)
c           q3m(j) = (-dym(j) - dzm(j) + 2*dxm(j))/sqrt(6.d0)
c           write(22,101) ind(j),q1p(j),q2p(j),q3p(j),
c     &q1m(j),q2m(j),q3m(j)
c        end do
c        write(22,101) 18,q1p(1),q2p(1),q3p(1),
c     &q1m(1),q2m(1),q3m(1)
cc
cc    re-re-initialisation of the van Vleck modes
cc
c        do j=1,18
c           dxp(j) = 0
c           dyp(j) = 0
c           dzp(j) = 0
c           dxm(j) = 0
c           dym(j) = 0
c           dzm(j) = 0
c           q1p(j) = 0
c           q2p(j) = 0
c           q3p(j) = 0
c           q1m(j) = 0
c           q2m(j) = 0
c           q3m(j) = 0
c        end do
cc
cc    re-re-computation of the van Vleck modes
cc
c        write(6,*) 'writing on disk ... 312'
c        do j=1,18
c           dxp(j) = lxp(j) - avgtmp
c           dyp(j) = lyp(j) - avgtmp
c           dzp(j) = lzp(j) - avgtmp
c           dxm(j) = lxm(j) - avgtmm
c           dym(j) = lym(j) - avgtmm
c           dzm(j) = lzm(j) - avgtmm
c           q1p(j) = ( dzp(j) + dxp(j) + dyp(j))/sqrt(3.d0)
c           q2p(j) = ( dzp(j) - dxp(j))/sqrt(2.d0)
c           q3p(j) = (-dzp(j) - dxp(j) + 2*dyp(j))/sqrt(6.d0)
c           q1m(j) = ( dzm(j) + dxm(j) + dym(j))/sqrt(3.d0)
c           q2m(j) = ( dzm(j) - dxm(j))/sqrt(2.d0)
c           q3m(j) = (-dzm(j) - dxm(j) + 2*dym(j))/sqrt(6.d0)
c           write(23,101) ind(j),q1p(j),q2p(j),q3p(j),
c     &q1m(j),q2m(j),q3m(j)
c        end do
c        write(23,101) 18,q1p(1),q2p(1),q3p(1),
c     &q1m(1),q2m(1),q3m(1)

        write(6,*) 'done!'
c
c
101     format(4x,i3,5x,f11.5,7x,f11.5,7x,f11.5,8x,
     &f11.5,7x,f11.5,7x,f11.5)
        close(11)
        close(21)
        close(22)
        close(23)
c
        stop
        end
