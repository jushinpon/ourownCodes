	MODULE INFORMATION
	implicit real*8(a-h,o-z)

	character*30 name
	
      real*8,allocatable::x(:),y(:),z(:),distance(:,:)      
	integer,allocatable::binpnt(:),bin(:),Itype(:),ncount(:)
     &,neighb(:,:)
	integer natom,ID,max_neigh
	real*8 xl,yl,zl, half_xl,half_yl,half_zl
	real*8 cutoff,prefactor,volsum
	logical pbcx, pbcy, pbcz
      real*8 rlistx,rlisty,rlistz,rlistsq

    

	contains
c*********************************************
      subroutine general_para

C******** PBC conditions for your system *********

      open(1,file="inputdata.dat",status='old')
      read(1,*)name,cutoff ! from the rdf profile (you may get it by ovito)
      read(1,*)name,pbcx
	read(1,*)name,pbcy
	read(1,*)name,pbcz
	read(1,*)name,max_neigh
	close(1)	
       
      rlistsq = cutoff*cutoff
	prefactor = 4.d0*3.1415926/(3.d0*8.d0)
	return
	endsubroutine general_para 
c***************************************************
	 
	ENDMODULE
