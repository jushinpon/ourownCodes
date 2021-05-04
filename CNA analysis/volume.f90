!aint(): truncates a real number by removing its fractional part
!The function NINT rounds to the nearest whole number    	
!The function INT will:Leave an integer unchanged.Round a real number towards 0.

!first, second, third, fourth neighbour IDs 
!CN_No(i,j)? i: atom id, j: number of the first neighbor atom type 
!CN_ID(i,j,k)? i: atom id, j: neighbor type, k: counter --> ID --> CNA_type
! allocate in information.f95
subroutine volume
use information
implicit real*8(a-h,o-z)     

nxcell = dint(xl/rlist) !get the integer, and remove the decimal.
temp = xl - dble(nxcell) * rlist
if (temp .gt. 0.0)  nxcell = nxcell + 1
nycell = dint(yl/rlist)
temp = yl - dble(nycell) * rlist
if (temp .gt. 0.0)  nycell = nycell + 1
nzcell=dint(zl/rlist)
temp = zl - dble(nzcell) * rlist
if (temp .gt. 0.0)  nzcell = nzcell + 1    

if(nxcell .lt. 3 .or. nycell .lt. 3 .or.nzcell .lt. 3) then
    write(*,*)"Cells are fewer than 3. STOP here!!"
	write(*,*)"nxcell",nxcell,"nycell",nycell,"nzcell",nzcell
	stop
endif

rlistsq = rlist*rlist
rlistinv = 1.d0/rlist
allocate(binpnt(nxcell*nycell*nzcell))
       
binpnt=0 ! array for each cell ID
bin =0

do 2 i=1,natom

	ix = aint(x(i)*rlistinv) + 1
    if(ix .gt. nxcell)ix = nxcell

    iy = aint(y(i)*rlistinv) + 1
    if(iy .gt. nycell)iy = nycell 

	iz=aint(z(i)*rlistinv) + 1
    if(iz .gt. nzcell)iz = nzcell

    ib=(iz-1)*(nxcell*nycell) + (iy-1)*nxcell + ix
!ib is the label for a specific cell  atom i belongs to......
   	 bin(i) = binpnt(ib)
	 binpnt(ib) = i ! finally, binpnt(ib) keeps the lagest ID in this cell
!example, atoms 2,5,7 in cell 1 (ib = 1)
!ib=1,i = 2 --> bin(2) = binpnt(1) = 0 (first time), binpnt(1) = 2
!ib=1,i = 5 --> bin(5) = binpnt(1) = 2 (from the above value), binpnt(1) = 5
!ib=1,i = 7 --> bin(7) = binpnt(1) = 5 (from the above value), binpnt(1) = 7
! The final value of binpnt(1) after the natom loop will keep the larger ID in cell 1
 2	continue

!CN_No(i,j)? i: atom id, j: atom type of the first neighbor atom 
!CN_ID(i,j,k)? i: atom id, j: r type, k: counter --> ID --> CNA_type
!distance(i,j) i: atom id, j: neighbor counter ID
CN_No = 0 ! number of each atom type
CN_ID = 0 ! atom id of the first neighbour atoms of each type 
distance = 0.d0 !distance between the ref atom i and its first neighbour atom j
do 3 i=1,natom

    index = 0 ! the first neighbour atom counter
	 
!c	 npoint(i)=nabors+1
	xtmp = x(i)													 
	ytmp = y(i)
	ztmp = z(i)
	
	ixx = aint(xtmp*rlistinv) + 1
	if(ixx .gt. nxcell) ixx = nxcell
	iyy = aint(ytmp*rlistinv) + 1
	if(iyy .gt. nycell) iyy = nycell
	izz=aint(ztmp*rlistinv) + 1
	if(izz .gt. nzcell) izz = nzcell

!c       write(*,*)"ixx,iyy,izz= ",ixx,iyy,izz
! go through nearest cells
	  do 4 k1=-1,1
	    do 5 k2=-1,1
	      do 6 k3=-1,1

	       ix = ixx + k1 ! get j's on i's cell id
	       if(ix.lt.1) ix = nxcell
	       if(ix.gt.nxcell) ix=1
           iy = iyy + k2
	       if(iy.lt.1) iy = nycell
	       if(iy.gt.nycell) iy = 1
	       iz = izz + k3
           if(iz.lt.1) iz = nzcell
	       if(iz.gt.nzcell) iz = 1
           !bin ID of atom j's cell
     	   ib=(iz-1)*(nxcell)*(nycell)+(iy-1)*nxcell+ix
	       j = binpnt(ib)	! the higest ID in cell ib          
!cell without atoms or atoms in this cell have been gone through
10            if(j.ne.0)then      
                 if(j.ne.i)then !in the same cell
                    dx=xtmp-x(j)
                    dy=ytmp-y(j)
                    dz=ztmp-z(j) 
					if (abs(dx) > half_xl .and. pbcx) dx = dx - sign(xl,dx)
					if (abs(dy) > half_yl .and. pbcy) dy = dy - sign(yl,dy)
					if (abs(dz) > half_zl .and. pbcz) dz = dz - sign(zl,dz)
                    rsq = dx*dx + dy*dy + dz*dz
                    if(rsq .le. rlistsq)then ! the first neighbour atom          
                        index = index + 1
                        if(index .gt. maxNB)then
							write(*,*)"the maximal neighbour atoms are more than current setting"
							stop
						endif
						r = dsqrt(rsq)
   	                     distance(i,index) = r
   	                   !  natom_type = Itype(j)
   	                     CN_ID(i,index) = j
                         !CN_CNAtype(i,index) = CNAtype(j)   	                     
                    endif
                  endif
                   j = bin(j)
				 goto 10
	            endif !j ne 0
 6	      continue !kx
 5	    continue !ky
 4	  continue !kz

      ncount(i)=index
 3	continue

 deallocate(binpnt)
!..................... finish cell list	
atomVol = 0.d0 
do 771 i=1,natom
   sumrsqinv = 0.d0
   sumrinv = 0.d0
        do 772 j=1,ncount(i)

			disinv = 1.d0/distance(i,j)         
			dissqinv = disinv*disinv
	
			sumrsqinv=sumrsqinv+dissqinv
			sumrinv = sumrinv+disinv

772     continue
     
	if( sumrsqinv .NE. 0.d0)then ! no neighbour atom
      ai2 = sumrinv/sumrsqinv
	  atomVol(i) = prefactor*ai2*ai2*ai2
	endif
	
771   continue

return
end  
