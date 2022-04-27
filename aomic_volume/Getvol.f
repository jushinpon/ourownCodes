C     This code do the postprocess of lammps dumped cfg for getting the atomic vol and then the system vol by sum the atomic vol.
C     Developed by Prof. Shin-Pon Ju 2016/10/19

c      V_i=4�k/3*a_i^3,a_i=1/2*(sum(rij^(-1) )/sum(rij^(-2)): atomic vol equation from the paper below 
c     D. Srolovitz, K. Maeda, V. Vitek, and T. Egami, "Structural defects in amorphous solids statistical analysis of a computer model
C    ," Philosophical Magazine A, vol. 44, pp. 847-866, 1981.
	
C     data files you should prepare
c.    1. input.cfg
c.    2. finputdata.dat (parameters for this code)

C.    output file:cfg_vol.dat

C..... You need use the same title for each file and use a Perl script to go through every cfg file (rename your cfg to be input cfg)

	include 'information.f'
	PROGRAM Getvol
	USE information
	implicit real*8(a-h,o-z)      
	
      call general_para ! make parameters workable		   
	              
      open(112,file="input.cfg",status='old')
c........... read cfg file		    
      read(112,*)
	read(112,*)ntime
	read(112,*)
      read(112,*)natom
 
      allocate(x(natom))
	allocate(y(natom))
      allocate(z(natom))
      allocate(distance(natom,max_neigh))! assume the maxmial coordinate number is 16
	allocate(Itype(natom))
	allocate(bin(natom))
      allocate(ncount(natom))			  			  			 
     
      read(112,*)
      read(112,*)xlo, xhi
	xl=(xhi-xlo)
	half_xl = xl/2.
      read(112,*)ylo, yhi
	yl=(yhi-ylo)
	half_yl = yl/2.
      read(112,*)zlo, zhi
	zl=(zhi-zlo)
	half_zl = zl/2.	
      read(112,*)

      do 117 j2=1,natom
        read(112,*)ID,Itype(ID),x(ID),y(ID),z(ID) !! ******** You should check the format for your casec       
C .............. move the origin to (0,0,0)
        x(ID) = x(ID) -xlo
	  y(ID) = y(ID) -ylo
	  z(ID) = z(ID) -zlo     
117   continue      
	close(112)
c      get the sum of each atomic vol
           
      call volume
c #     write(*,*)"Vol", volsum
      open(112,file="cfg_vol.dat",status='unknown')
 	write(112,*)volsum
	close(112)

	END

	
	include 'volume.f'