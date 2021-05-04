!     This code do the postprocess of lammps dumped cfg for getting the atomic vol and then the system vol by sum the atomic vol.
!     Developed by Prof. Shin-Pon Ju 2016/10/19
!
!      V_i=4Pi/3*a_i^3,a_i=1/2*(sum(rij^(-1) )/sum(rij^(-2)): atomic vol equation from the paper below 
!     D. Srolovitz, K. Maeda, V. Vitek, and T. Egami, "Structural defects in amorphous solids statistical analysis of a computer model
!    ," Philosophical Magazine A, vol. 44, pp. 847-866, 1981.
!	
!     data files you should prepare
!.    1. input.cfg
!.    2. finputdata.dat (parameters for this code)
!
!.    output file:cfg_vol.dat
!
!..... You need use the same title for each file and use a Perl script to go through every cfg file (rename your cfg to be input cfg)

include 'information.f90'
PROGRAM Getvol
USE information
implicit real*8(a-h,o-z)      


call general_para ! make parameters workable		   

allocate(CNAVol(CNAtypeNo))
allocate(CNAVolmax(CNAtypeNo))
allocate(CNAVolmin(CNAtypeNo))
allocate(CNAtype_counter(CNAtypeNo))
allocate(CN_CNAtype_counter(CNAtypeNo))
allocate(atomtype_of_CNA(atomtypeNo)) ! atom type fraction of a CNA type
              
open(112,file=trim(cfg_file),status='old')
!. read cfg file		    
read(112,*)
read(112,*)ntime
read(112,*)
read(112,*)natom
 
allocate(x(natom))
allocate(y(natom))
allocate(z(natom))
allocate(distance(natom,maxNB))! assume the maxmial coordinate number is maxNB
!allocate(CN_CNAtype(natom,CNAtypeNo))! neighbour atom CNA type number
allocate(Itype(natom)) !atom type in lammps
allocate(bin(natom)) ! for list bin
allocate(ncount(natom))			  			  			 
allocate(CNAtype(natom)) !CNA type of each atom			  			  			 
allocate(CN_ID(natom,maxNB)) !the first neighbour atom IDs of each ref atom			  			  			 
!allocate(CN_ID(natom,atomtypeNo,CNAtypeNo))			  			  			 
allocate(atomVol(natom)) !for atomic volume term			  			  			 
   
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
	read(112,*)ID,Itype(ID),x(ID),y(ID),z(ID),CNAtype(ID) !! ******** You should check the format for your casec       
	CNAtype(ID) = CNAtype(ID) + 1 ! make the lowest value 1
!..... move the origin to (0,0,0)
    x(ID) = x(ID) -xlo
	y(ID) = y(ID) -ylo
	z(ID) = z(ID) -zlo     
117   continue      
close(112)

!get the volume of each atom           
call volume
CNAVolmax = -100.d0
CNAVolmin = 100.d0

CNAVol = 0.d0
CNAtype_counter = 0
CN_CNAtype_counter = 0 !neighbour atom CNA counter
CN_CNAtype = 0 !CNA type Number of neighbour atoms 
ICNA_type = 0 ! counter for a specific CNA type
atomtype_of_CNA = 0 ! atom types of a specific CNA type
do i = 1,natom
	nty = CNAtype(i)
	CNAVol(nty) = CNAVol(nty) + atomVol(i)
	CNAtype_counter(nty) = CNAtype_counter(nty) + 1
	if(atomVol(i) .gt. CNAVolmax(nty)) CNAVolmax(nty) = atomVol(i)  
	if(atomVol(i) .lt. CNAVolmin(nty)) CNAVolmin(nty) = atomVol(i)  
if(nty .eq. 2) then	
	ICNA_type =  ICNA_type + 1 
	atomtype_of_CNA(itype(i))  = atomtype_of_CNA(itype(i)) + 1
	do j = 1, ncount(i)
		jid = CN_ID(i,j)
		jcnatype = CNAtype(jid)
		CN_CNAtype_counter(jcnatype) = CN_CNAtype_counter(jcnatype) + 1
	enddo
endif
	
enddo
!     write(*,*)"Vol", volsum
open(112,file=trim(output_file),status='unknown')
write(112,*)"CNA analysis summary by rlist ",rlist," ->0:other,1:FCC, HCP:2,BCC:3,ICO:4" 
write(112,*)
write(112,*)"CNAtype frac aveVol maxVol minVol"

do ity = 1, CNAtypeNo 
	aveVol = dble(CNAVol(ity))/dble(CNAtype_counter(ity))
	frac = dble(CNAtype_counter(ity))/dble(natom)
	write(112,*)(ity - 1), frac, aveVol, CNAVolmax(ity),CNAVolmin(ity)
enddo

write(112,*)
write(112,*) "Average CNAtype number of the first neighbour atoms of a ref FCC atom"

do ity1 = 1, CNAtypeNo	
	aveCNA = dble(CN_CNAtype_counter(ity1))/dble(CNAtype_counter(2))
	write(112,*)(ity1 -1),aveCNA
enddo

write(112,*)
write(112,*) "Atom type fraction of an FCC atom"

do iaty = 1, atomtypeNo	
	atomfrac = dble(atomtype_of_CNA(iaty))/dble(ICNA_type)
	write(112,*)iaty,atomfrac
enddo

close(112)

END
	
include 'volume.f90'
