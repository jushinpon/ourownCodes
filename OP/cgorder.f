c     This code was mainly developed by Prof. Shin-Pon Ju at NSYSU. It can only be used 
c     in Ju's group or with Ju's permission.

      Program main
!      include 'cxml_include.f90'
      implicit real*8(a-h,o-z)
      
      dimension molecule(500000),x(1000000),y(1000000),z(1000000)
     &	,TOP(200)
	dimension amass(1000000),rotiner(3,3),Q(3,3),op(1000),gy(3,3),
     &fmolmasscenterx(100000),fmolmasscentery(100000)
     &,fmolmasscenterz(100000)
	dimension fenx(100000),feny(100000),fenz(100000)!inertia eigenvector
      dimension fnenx(1000),fneny(1000),fnenz(1000) !nematic eigenvector
      dimension gyrx(1000),gyry(1000),gyrz(1000),gyeig(3) !gyration eigenvector
      character*2 el(1000000)
      CHARACTER JOBZ, UPLO
      INTEGER INFO, LDA, LWORK, N
	parameter (lwork=10*3)
      DOUBLE PRECISION A(3,3), W(3), WORK(lwork), evec(3)
	logical PBC

      PBC=.True.
c	PBC=.False.

      JOBZ='V'
	UPLO='U'
	N=3               !!!¶¥¼Æ
c	NEVEC=1
	LDA=3  



      open(119,file='1.dat',status='unknown')
	read(119,*)ifile
	write(*,*)"*************************************"
	write(*,*)"total file numbers in arc file:",ifile
	read(119,*)imolecule
      write(*,*)"total molecules in a file:",imolecule
	natom = 0
	do 1 i=1,imolecule
	   read(119,*)molecule(i)
c         write(*,*)"total atoms in molecule",i, ": ",molecule(i)
	   natom = natom+molecule(i)
1     continue
      write(*,*)"*************************************"

      nall = natom*ifile

c********read PBC lengths if needed 

      if(PBC) then
        read(119,*)boxlx
        read(119,*)boxly
        read(119,*)boxlz
        halflx = boxlx/2.      
        halfly = boxly/2.
        halflz = boxlz/2.
	endif
C**********   

      do 2 i=1,nall 
         read(119,*)x(i),y(i),z(i),el(i)
c         write(*,*)x(i),y(i),z(i),el(i)
c	pause
2     continue
      close(119)

C*********** handle PBC
      If(PBC) then

       do 41 i=1,ifile
	   nfatom=natom*(i-1) ! natom: total atoms in a file
	   nfmolecule= imolecule*(i-1)
	   nmolecule = 0
	     do 51 j=1,imolecule ! molecule counter in a file
	       
              if(j.gt.1)nmolecule=nmolecule+molecule(j-1) ! because not molecule(0) if j=1
                     index = nfatom+nmolecule+1

                 do 61 k=2,molecule(j)    ! atom counter in a chain     
      	           index = nfatom+nmolecule+k
                     xtemp = x(index-1)
                     ytemp = y(index-1)
                     ztemp = z(index-1)


c c                    xmc=xmc+amass(index)*x(index)
c                     ymc=ymc+amass(index)*y(index) 
c                     zmc=zmc+amass(index)*z(index) 
c                     totmass = totmass+amass(index)

                    delx=x(index)-xtemp
                    dely=y(index)-ytemp
                    delz=z(index)-ztemp

				  if(delx .lt. -halflx) x(index)=x(index)+boxlx
	              if(delx .gt. halflx) x(index)=x(index)-boxlx
	              if(dely.lt. -halfly)  y(index)=y(index)+boxly
	              if(dely .gt. halfly) y(index)=y(index)-boxly
                    if(delz.lt. -halflz)  z(index)=z(index)+boxlz
	              if(delz .gt. halflz) z(index)=z(index)-boxlz	

61               continue
51        continue             
41     continue             



	endif

C********************** 
      


      do 3 i=1,nall      
         if(el(i) .eq. "C2")then
	     amass(i)=29.0620
	   elseif(el(i) .eq. "C3")then
           amass(i)=42.0810
	   elseif(el(i) .eq. "C61")then
	     amass(i)=76.0980
         elseif(el(i) .eq. "C62")then
	     amass(i)=76.0980
	   elseif(el(i) .eq. "CN")then
	     amass(i)=26.0180
	   endif
3     continue

      do 4 i=1,ifile
	   nfatom=natom*(i-1)
	   nfmolecule= imolecule*(i-1)
	   nmolecule = 0
	fmolmasscenterx=0.d0 ! get mass center of each molecule
      fmolmasscentery=0.d0
	fmolmasscenterz=0.d0

	     do 5 j=1,imolecule
	        xmc = 0.d0
              ymc = 0.d0
	        zmc = 0.d0
	        totmass = 0.d0
              if(j.gt.1)nmolecule=nmolecule+molecule(j-1)
                 do 6 k=1,molecule(j)         
      	           index = nfatom+nmolecule+k
                     xmc=xmc+amass(index)*x(index)
                     ymc=ymc+amass(index)*y(index) 
                     zmc=zmc+amass(index)*z(index) 
                     totmass = totmass+amass(index)
6                continue
              xmc=xmc/totmass
              ymc=ymc/totmass
              zmc=zmc/totmass 
	        fmolmasscenterx(j)= xmc
			fmolmasscentery(j)= ymc
			fmolmasscenterz(j)= zmc
			 
              rotiner= 0.d0
              gy =0.d0
                 do 7 k=1,molecule(j)
                    index = nfatom+nmolecule+k
	              disx=(x(index)-xmc)
                    disy=(y(index)-ymc)
                    disz=(z(index)-zmc)
	              dissq = disx*disx+disy*disy+disz*disz      	
                    rotiner(1,1)= rotiner(1,1)+amass(index)*
     &     				        (dissq-disx*disx)
                    rotiner(2,2)= rotiner(2,2)+amass(index)*
     &                            (dissq-disy*disy)
                    rotiner(3,3)= rotiner(3,3)+amass(index)*
     &                            (dissq-disz*disz)
                    rotiner(1,2)= rotiner(1,2)+amass(index)*
     &                            (-disx*disy)
                    rotiner(1,3)= rotiner(1,3)+amass(index)*
     &                            (-disx*disz)
	              rotiner(2,3)= rotiner(2,3)+amass(index)*
     &                            (-disy*disz)

                    if(imolecule.eq.1)then
	                 gy(1,1)= gy(1,1)+disx*disx
	                 gy(2,2)= gy(2,2)+disy*disy
	                 gy(3,3)= gy(3,3)+disz*disz
	                 gy(1,2)= gy(1,2)+disx*disy
	                 gy(1,3)= gy(1,3)+disx*disz
	                 gy(2,3)= gy(2,3)+disy*disz
	              endif 
7                continue
              rotiner(2,1) =rotiner(1,2)
              rotiner(3,1) =rotiner(1,3)
              rotiner(3,2) =rotiner(2,3)
              A=rotiner
              call DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO ) 
              molindex = nfmolecule+j
              fenx(molindex)=A(1,1)
              feny(molindex)=A(2,1)
              fenz(molindex)=A(3,1)
              if(imolecule.eq.1)then
	           gy(2,1) =gy(1,2)
                 gy(3,1) =gy(1,3)
                 gy(3,2) =gy(2,3)
	           write(*,*)"molindex = ",molindex
                 A=gy/dble(molecule(molindex))
              call DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
                 molindex = nfmolecule+j
                 gyrx(molindex)=A(1,3)
                 gyry(molindex)=A(2,3)
                 gyrz(molindex)=A(3,3)
                 trgy1x=A(1,2)
	           trgy1y =A(2,2)
	           trgy1z =A(3,2)
                 trgy2x=A(1,1)
	           trgy2y =A(2,1)
	           trgy2z =A(3,1)
                 gyeig(1) = w(3)
                 gyeig(2) = w(2)
                 gyeig(3) = w(1)
	        endif
5          continue      
4     continue

c     getting nematic direction for each file!!!

      do 10 i = 1,ifile
        nfmolecule= imolecule*(i-1)
	  Q=0.d0
	   do 11 j= 1,imolecule
           molindex = nfmolecule+j
           Q(1,1)= Q(1,1)+1.5d0*fenx(molindex)*fenx(molindex)-0.5d0
           Q(2,2)= Q(2,2)+1.5d0*feny(molindex)*feny(molindex)-0.5d0
           Q(3,3)= Q(3,3)+1.5d0*fenz(molindex)*fenz(molindex)-0.5d0
	     Q(1,2)= Q(1,2)+1.5d0*(fenx(molindex)*feny(molindex))
           Q(1,3)= Q(1,3)+1.5d0*(fenx(molindex)*fenz(molindex))
           Q(2,3)= Q(2,3)+1.5d0*(feny(molindex)*fenz(molindex))
11       continue
        Q(2,1) =Q(1,2)
        Q(3,1) =Q(1,3)
        Q(3,2) =Q(2,3)
        A=Q/dble(imolecule)
        call DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
	
c	write(*,*)
      
        fnenx(i)=A(1,3)
        fneny(i)=A(2,3)
        fnenz(i)=A(3,3)

	  tran1x = A(1,2)
	  tran1y = A(2,2)
	  tran1z = A(3,2) 
        tran2x = A(1,1)
	  tran2y = A(2,1)
	  tran2z = A(3,1)
	  open(227,file='majoraxis.dat',status='unknown')
	  write(227,*)"major axis from moment of inertia"
	  write(227,*)fnenx(i),fneny(i),fnenz(i)
	  write(227,*)"transver1"
	  write(227,*)tran1x,tran1y,tran1z
        write(227,*)"transver2"
        write(227,*)tran2x,tran2y,tran2z 
10    continue

c     getting order parameter for each file!!
      op=0.d0
      open(120,file='orpara.dat',status='unknown')

	if(imolecule.eq.1)open(127,file='majoraxis.dat',status='unknown')
      
      do 13 i=1,ifile
         nfmolecule= imolecule*(i-1)
         cossq=0.d0	 
           do 14 j = 1,imolecule
              molindex = nfmolecule+j
              cosval=fenx(molindex)*fnenx(i)+feny(molindex)*fneny(i)+
     &               fenz(molindex)*fnenz(i)
              cossq=cossq+cosval*cosval
14         continue
         cossq=cossq/dble(imolecule)
         op(i)=1.5d0*cossq-0.5d0
	   write(120,*)i,op(i), "Orientational OP"	
         if(imolecule.eq.1)then
	   write(127,*)"major axis from moment of inertia"
	   write(127,*)fnenx(i),fneny(i),fnenz(i)
         write(127,*)"transver1"
         write(127,*)tran1x,tran1y,tran1z
         write(127,*)"transver2"
         write(127,*)tran2x,tran2y,tran2z
	   endif
C     translational OP here (THE JOURNAL OF CHEMICAL PHYSICS 130, 234501 .2009.)
        
	   TOP = 0
      do 555 itop =1,200
         nfatom=natom*(i-1)
	   nfmolecule= imolecule*(i-1)
	   nmolecule = 0
	   index = 0
	   d= 0.25*dble(itop)
         TOPcos=0.
	   TOPsin = 0. !! sumation of translational OP
       do 511 j1=1,imolecule ! molecule ID
	        
              
                         
C      	          fmolmasscenterz(j)

	              TOPcrossprod =fnenx(i)*fmolmasscenterx(j1)
     &				  +fneny(i)*fmolmasscentery(j1)
     &				  +fnenz(i)*fmolmasscenterz(j1)

	           TOPcos= TOPcos+cos(6.2831852*TOPcrossprod/d)
                 TOPsin= TOPsin+sin(6.2831852*TOPcrossprod/d)                     

     
511      continue  
         TOP(itop)=dsqrt( (TOPcos)*(TOPcos)+
     &    (TOPsin)*(TOPsin) )/dble(imolecule)                                       )
555     continue ! end of different d trial test

      write(120,*)i,maxval(top), "Translational OP"
c         write(*,*)"TOP test"
c	   do io = 1,200
c         write(*,*)io,' ',TOP(io)
c	   enddo
c	   write(*,*)maxval(top)
c	   write(*,*)fnenx(i),fneny(i),fnenz(i)
c	   pause 
c	   write(127,*)"************  "
c         write(127,*)"major axis from gyration tensor"
c	   write(127,*)gyrx(i),gyry(i),gyrz(i)
c         write(127,*)"transver1"
c         write(127,*)trgy1x,trgy1y,trgy1z
c         write(127,*)"transver2"
c         write(127,*)trgy2x,trgy2y,trgy2z 
c	   write(127,*)"************  "
c	   write(127,*) "eigenvalues from gyration tensor"
c	   write(127,*)gyeig(1),gyeig(2),gyeig(3)
c         write(127,*)"************  "
c	   write(127,*) "axial lengthes"
c	   write(127,*)dsqrt(abs(gyeig(1)))*2.d0,dsqrt(abs(gyeig(2)))
c     &         	   *2.d0,dsqrt(abs(gyeig(3)))*2.d0

	   	
13    continue
      close(120)

      if(imolecule.eq.1)close(127)
      open(121,file='opsummary.dat',status='unknown')
      opmax =maxval(op)
      opmin =opmax
	opave = 0.d0
	do 15 i=1,ifile
         opave=opave+op(i)
         if(op(i).lt.opmin)opmin=op(i)
15    continue
      opave=opave/dble(ifile)
c	write(*,*)opmax,opmin
      write(121,*)"Averaged OP :",opave
      write(121,*)"Largest OP :",opmax
      write(121,*)"Smallest OP :",opmin

      topmax = maxval(TOP)
	id = 0
      do i=1,200
      if(topmax .eq. TOP(i)) id=i 
	enddo


      write(121,*)"*****"
      write(121,*)"Translational OP :",topmax
      write(121,*)"The layer spacing of TOP :",0.25*dble(id)
	close(121)
      stop
	end
