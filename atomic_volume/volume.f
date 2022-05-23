  	subroutine volume
	use information
 	implicit real*8(a-h,o-z)     
	real*8 dx,dy,dz
	
	rlisttest=cutoff
      nlistx=aint(xl/rlisttest) !get the integer, and remove the decimal.
	nlisty=aint(yl/rlisttest)
	nlistz=aint(zl/rlisttest)    


      rlistx=rlisttest+(xl-dble(nlistx)*rlisttest)/dble(nlistx) ! get the cell with integer times
      rlisty=rlisttest+(yl-dble(nlisty)*rlisttest)/dble(nlisty)
	rlistz=rlisttest+(zl-dble(nlistz)*rlisttest)/dble(nlistz)

	rlistsqx=rlistx*rlistx
	rlistsqy=rlisty*rlisty
	rlistsqz=rlistz*rlistz
	
	nxcell=nint(xl/rlistx) ! get the closest integer
	nycell=nint(yl/rlisty)
	nzcell=nint(zl/rlistz)

      allocate(binpnt(nxcell*nycell*nzcell))
       
      binpnt=0
      bin =0

	rlistxinv=1.d0/rlistx
	rlistyinv=1.d0/rlisty
	rlistzinv=1.d0/rlistz	
		

	do 2 i=1,natom
	 ix=aint(x(i)*rlistxinv)+1
       if(ix .gt. nxcell)ix=ix-1
	 iy=aint(y(i)*rlistyinv)+1
       if(iy .gt. nycell)iy=iy-1 
	 iz=aint(z(i)*rlistzinv)+1
       if(iz .gt. nzcell)iz=iz-1

	 ib=(iz-1)*(nxcell*nycell)+(iy-1)*nxcell+ix
!.... ib is the label for a specific cell  atom i belongs to......
   	 bin(i)=binpnt(ib)
	 binpnt(ib)=i ! finally, binpnt(ib) keeps the lagest ID in this cell
 2	continue
	
!	pause
!      if(npre.eq.165)write(*,*)"force -3"

	do 3 i=1,natom

       index=0
	 
!	 npoint(i)=nabors+1
	 xtmp=x(i)													 
	 ytmp=y(i)
	 ztmp=z(i)
	
	 ixx=aint(xtmp*rlistxinv)+1
	 if(ixx .gt. nxcell)ixx=ixx-1
	 iyy=aint(ytmp*rlistyinv)+1
	 if(iyy .gt. nycell)iyy=iyy-1
	 izz=aint(ztmp*rlistzinv)+1
	 if(izz .gt. nzcell)izz=izz-1


!       write(*,*)"ixx,iyy,izz= ",ixx,iyy,izz
	  do 4 k1=-1,1
	    do 5 k2=-1,1
	      do 6 k3=-1,1

	       ix=ixx+k1
	       if(ix.lt.1) ix=nxcell
	       if(ix.gt.nxcell) ix=1

             iy=iyy+k2
	       if(iy.lt.1) iy=nycell
	       if(iy.gt.nycell) iy=1

	       iz=izz+k3
             if(iz.lt.1) iz=nzcell
	       if(iz.gt.nzcell) iz=1

     	       ib=(iz-1)*(nxcell)*(nycell)+(iy-1)*nxcell+ix
	       j=binpnt(ib)	! the higest ID in cell ib          

10            if(j.ne.0)then        

                 if(j.ne.i)then
                    dx=xtmp-x(j)
                    dy=ytmp-y(j)
                    dz=ztmp-z(j) 
      			if (abs(dx) > half_xl .and. pbcx) then
				  dx = dx - sign(xl,dx)
				endif
				if (abs(dy) > half_yl .and. pbcy) then
					dy = dy - sign(yl,dy)
				endif
	            if (abs(dz) > half_zl .and. pbcz) then 
					dz = dz - sign(zl,dz)
				endif
!                  if (abs(dy) > half_yl .and. pbcy) dy = dy - sign(yl,dy)
!                  if (abs(dz) > half_zl .and. pbcz) dz = dz - sign(zl,dz)    
                    rsq=dx*dx + dy*dy + dz*dz

                    if(rsq .le. rlistsq)then           

                         index=index+1

                         r=dsqrt(rsq)
!                        neighb(i,index)=j	
					if(index .gt. max_neigh)then
		write(*,*)''
		write(*,*)'***********************'
		write(*,*)'The index parameter over maximum you assign'
		write(*,*)'maximum you assign: ',max_neigh
		write(*,*)'current value: ',index
		write(*,*)''
					stop
					endif	
 	                  distance(i,index)=r
	                                         
                    endif
                  

                  endif

                   j=bin(j)
				 goto 10

	            endif
 6	      continue
 5	    continue
 4	  continue

      ncount(i)=index
 3	continue

      deallocate(binpnt)
!..................... finish cell list	
 
      volsum =0.d0
      
	do 771 i=1,natom
        sumrsqinv = 0.d0
	  sumrinv = 0.d0
	 
        do 772 j=1,ncount(i)

        disinv = 1.d0/distance(i,j)         
        dissqinv = disinv*disinv

        sumrsqinv=sumrsqinv+dissqinv
        sumrinv = sumrinv+disinv

772     continue
     
	if( sumrsqinv .NE. 0.d0)then
      ai2 = sumrinv/sumrsqinv
	atomvol = prefactor*ai2*ai2*ai2
!	write(*,*)"atom ",i," vol = ",atomvol
!	write(*,*)"ai = ",ai2/2.
      volsum = volsum+atomvol
	endif
!      write(*,*)i,volsum
	
771   continue
	
!      write(*,*)"volsum",volsum


	return
	end  