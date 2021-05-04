MODULE INFORMATION
implicit real*8(a-h,o-z)

character*30 string
character*72 cfg_file,output_file
real,allocatable::x(:),y(:),z(:),distance(:,:),atomVol(:)      
integer,allocatable::binpnt(:),bin(:),Itype(:),ncount(:),neighb(:,:)
integer natom,ID,CNAtypeNo,atomtypeNo,maxNB
real*8 xl,yl,zl, half_xl,half_yl,half_zl
real*8 cutoff,prefactor,volsum
logical pbcx, pbcy, pbcz
real*8 rlist,rlistsq
integer,allocatable::CN_ID(:,:),CNAtype(:) !,CN_CNAtype(:,:)
real*8,allocatable:: CNAVol(:),CNAVolmax(:),CNAVolmin(:)
integer,allocatable::CNAtype_counter(:),CN_CNAtype_counter(:),atomtype_of_CNA(:)

contains
!*********************************************
subroutine general_para

!******** PBC conditions for your system *********

open(1,file="./Finputdata.dat",status='old')
read(1,*)string,cfg_file 
read(1,*)string,output_file 
read(1,*)string,atomtypeNo ! the CNA type number you consider
read(1,*)string,CNAtypeNo ! the CNA type number you consider
read(1,*)string,rlist ! from the rdf profile (you may get it by ovito)
read(1,*)string,pbcx
read(1,*)string,pbcy
read(1,*)string,pbcz
close(1)	
       
rlistsq = rlist*rlist
prefactor = 4.d0*3.1415926/(3.d0*8.d0)
maxNB = 20
return
endsubroutine general_para 
!***************************************************
	 
ENDMODULE
