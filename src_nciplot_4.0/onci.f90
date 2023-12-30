module onci
   use reader
   use props
   use tools_io
   use param

   implicit none

   public 

contains

   subroutine write_cube_header_onci(lu, l1, l2, ntotal, xinit, xinc, nstep, m, nfile)

      integer, intent(in) :: lu, nfile
      character*(*), intent(in) :: l1, l2
      integer,intent(in)   :: ntotal, nstep(3)
      real*8, intent(in)   :: xinit(3), xinc(3)
      type(molecule),intent(in) :: m
      integer :: i, j

      write (lu, *) trim(l1)
      write (lu, *) trim(l2)
      write (lu, '(I5,3(F12.6))') ntotal, xinit

      write (lu, '(I5,3(F12.6))') nstep(1), xinc(1), 0d0, 0d0
      write (lu, '(I5,3(F12.6))') nstep(2), 0d0, xinc(2), 0d0
      write (lu, '(I5,3(F12.6))') nstep(3), 0d0, 0d0, xinc(3)
      do i = 1, nfile
         do j = 1, m%n
            write (lu, '(I4,F5.1,F11.6,F11.6,F11.6)') m%z(j), 0d0, m%x(:, j)
         end do
      enddo

   end subroutine write_cube_header_onci

   subroutine write_cube_body_onci(lu, n, c)

      integer, intent(in) :: lu
      integer, intent(in) :: n(3)
      real*8, intent(in) :: c(0:n(1) - 1, 0:n(2) - 1, 0:n(3) - 1)
      
      integer :: i, j, k

      do i = 0, n(1) - 1
         do j = 0, n(2) - 1
            write (lu, '(6(1x,e12.5))') (c(i, j, k), k=0, n(3) - 1)
         enddo
      enddo
      close (lu)

   end subroutine write_cube_body_onci

subroutine calcsab(x,m,rho,l,p,csab)
     
       
     real*8, intent(in)       :: rho, x(3)
     type(molecule),intent(in):: m
     integer, intent(in)      :: l, p
     real*8, intent(out)      :: csab
     
     real*8                   ::store(2,10), sab_3D

     !Value of the orbital l and p of molecule m at point x
     call phi_sab(x,m,l,p,store)

     sab_3D=(32.d0*dot_product(store(1,1)*store(1,2:4),store(2,1)*store(2,2:4)))/((abs(rho)/100.d0)**(8.d0/3.d0))
     csab=sab_3D

end subroutine calcsab

subroutine orb_density(xinit, xinc, n, mol, nmol, rho_orb)
        use reader
        use tools_io
        use tools_math
        use param
        use props 

        real*8, intent(in) :: xinit(3), xinc(3)
        integer, intent(in) :: n(3), nmol
        type(molecule) :: mol(nmol)
        real*8, allocatable, dimension(:,:,:,:), intent(out) :: rho_orb
        real*8, allocatable, dimension(:, :) :: dx, dy, dz, d2, rhoaux
      real*8, allocatable, dimension(:) :: tp, maxc
      integer :: nmcent, istat
      real*8, allocatable :: chi(:, :), phi(:, :), hess(:, :, :)
      logical, allocatable :: ldopri(:, :)

      integer :: i, j, k, m, iat, ip, jp, kp, nn, l(3)
      integer :: ityp, ipri, ipria, ix, imo
      real*8 :: ex, xl2, xl(3, 0:2), x0(3), al
      real*8 :: wk1(3), wk2(3), hvecs(3, 3), grad2, heigs(3)
      real*8 :: rho53, ebose, df

      real*8, parameter :: cutoff_pri = 1d-10
      real*8, parameter :: fothirds = 4d0/3d0


      nmcent = 0
      do i = 1, nmol
         nmcent = max(nmcent, mol(i)%n)
      end do

      allocate (rho_orb(0:n(1)-1,0:n(2)-1,0:n(3)-1,mol(1)%nmo), stat=istat)
      if (istat /= 0) call error('calcprops', 'error allocating distance matrix', faterr)
      allocate (dx(n(1), nmcent), dy(n(1), nmcent), dz(n(1), nmcent), d2(n(1), nmcent), stat=istat)
      if (istat /= 0) call error('calcprops', 'error allocating distance matrix', faterr)
      allocate (hess(n(1), 3, 3), stat=istat)
      if (istat /= 0) call error('calcprops', 'error allocating hessian', faterr)
      allocate (tp(n(1)), rhoaux(n(1),mol(1)%nmo),  stat=istat)
      if (istat /= 0) call error('calcprops', 'error allocating rhoaux', faterr)

      rho_orb = 0d0
      hess = 0d0
      rhoaux = 0d0
      tp = 0d0
      dx(:,:) = 0d0
      dy(:,:) = 0d0
      dz(:,:) = 0d0
      d2(:,:) = 0d0

      do j = 0, n(2) - 1
         jp = j + 1
         do k = 0, n(3) - 1
            kp = k + 1
            ! zero in-line hessian and tp
            hess = 0d0
            rhoaux = 0d0
            tp = 0d0

            ! run over molecules
            ! Original
            do m = 1, nmol

               ! allocate primitive evaluation array
               allocate (chi(mol(m)%npri, 10), phi(mol(m)%nmo, 10))

               ! identify the max coefficient
               allocate (maxc(mol(m)%npri), ldopri(mol(m)%npri, 10))
               maxc = 0d0
               do imo = 1, mol(m)%nmo
                  do ipri = 1, mol(m)%npri
                     maxc(ipri) = max(maxc(ipri), abs(mol(m)%c(imo, ipri)))
                  enddo
               enddo

               ! calculate distances

               do iat = 1, mol(m)%n
                  do i = 0, n(1) - 1
                     ip = i + 1
                     dx(ip, iat) = xinit(1) + i*xinc(1) - mol(m)%x(1, iat)
                     dy(ip, iat) = xinit(2) + j*xinc(2) - mol(m)%x(2, iat)
                     dz(ip, iat) = xinit(3) + k*xinc(3) - mol(m)%x(3, iat)
                     d2(ip, iat) = dx(ip, iat)*dx(ip, iat) + dy(ip, iat)*dy(ip, iat) + dz(ip, iat)*dz(ip, iat)
                  enddo
               enddo

               ! calculate primitives at the points
               do i = 0, n(1) - 1
                  ip = i + 1
                  nn = 0
                  do ityp = 1, mol(m)%maxntyp
                     do ipria = nn + 1, nn + mol(m)%ntyp(ityp)
                        ipri = mol(m)%intyp(ipria)
                        iat = mol(m)%icenter(ipri)
                        al = mol(m)%e(ipri)
                        ex = exp(-al*d2(ip, iat))

                        x0 = (/dx(ip, iat), dy(ip, iat), dz(ip, iat)/)

                        call index0(ityp, l)
                        do ix = 1, 3
                           if (l(ix) == 0) then
                              xl(ix, 0) = 1d0
                              xl(ix, 1) = 0d0
                              xl(ix, 2) = 0d0
                           else if (l(ix) == 1) then
                              xl(ix, 0) = x0(ix)
                              xl(ix, 1) = 1d0
                              xl(ix, 2) = 0d0
                           else if (l(ix) == 2) then
                              xl(ix, 0) = x0(ix)*x0(ix)
                              xl(ix, 1) = 2d0*x0(ix)
                              xl(ix, 2) = 2d0
                           else if (l(ix) == 3) then
                              xl(ix, 0) = x0(ix)*x0(ix)*x0(ix)
                              xl(ix, 1) = 3d0*x0(ix)*x0(ix)
                              xl(ix, 2) = 6d0*x0(ix)
                           else if (l(ix) == 4) then
                              xl2 = x0(ix)*x0(ix)
                              xl(ix, 0) = xl2*xl2
                              xl(ix, 1) = 4d0*xl2*x0(ix)
                              xl(ix, 2) = 12d0*xl2
                           else
                              call error('pri012', 'power of L not supported', faterr)
                           end if
                        end do

                        chi(ipri, 1) = xl(1, 0)*xl(2, 0)*xl(3, 0)*ex
                        chi(ipri, 2) = (xl(1, 1) - 2*al*x0(1)**(l(1) + 1))*xl(2, 0)*xl(3, 0)*ex
                        chi(ipri, 3) = (xl(2, 1) - 2*al*x0(2)**(l(2) + 1))*xl(1, 0)*xl(3, 0)*ex
                        chi(ipri, 4) = (xl(3, 1) - 2*al*x0(3)**(l(3) + 1))*xl(1, 0)*xl(2, 0)*ex
                        chi(ipri, 5) = (xl(1, 2) - 2*al*(2*l(1) + 1)*xl(1, 0) &
                                        + 4*al*al*x0(1)**(l(1) + 2))*xl(2, 0)*xl(3, 0)*ex
                        chi(ipri, 6) = (xl(2, 2) - 2*al*(2*l(2) + 1)*xl(2, 0) &
                                        + 4*al*al*x0(2)**(l(2) + 2))*xl(3, 0)*xl(1, 0)*ex
                        chi(ipri, 7) = (xl(3, 2) - 2*al*(2*l(3) + 1)*xl(3, 0) &
                                        + 4*al*al*x0(3)**(l(3) + 2))*xl(1, 0)*xl(2, 0)*ex
                        chi(ipri, 8) = (xl(1, 1) - 2*al*x0(1)**(l(1) + 1))* &
                                       (xl(2, 1) - 2*al*x0(2)**(l(2) + 1))*xl(3, 0)*ex
                        chi(ipri, 9) = (xl(1, 1) - 2*al*x0(1)**(l(1) + 1))* &
                                       (xl(3, 1) - 2*al*x0(3)**(l(3) + 1))*xl(2, 0)*ex
                        chi(ipri, 10) = (xl(3, 1) - 2*al*x0(3)**(l(3) + 1))* &
                                        (xl(2, 1) - 2*al*x0(2)**(l(2) + 1))*xl(1, 0)*ex

                        do ix = 1, 10
                           ldopri(ipri, ix) = (abs(chi(ipri, ix))*maxc(ipri) > cutoff_pri)
                        enddo
                     enddo ! ipria = nn+1, nn+ntyp

                     nn = nn + mol(m)%ntyp(ityp)
                  enddo ! ityp = 1, maxntyp

                  ! build the MO values at the point
                  phi = 0d0
                  do ix = 1, 10
                     do ipri = 1, mol(m)%npri
                        if (.not. ldopri(ipri, ix)) cycle
                        do imo = 1, mol(m)%nmo
                           phi(imo, ix) = phi(imo, ix) + mol(m)%c(imo, ipri)*chi(ipri, ix)
                        enddo
                     enddo
                  enddo

                 ! contribution to the density, etc.
                  do imo = 1, mol(m)%nmo
                     rhoaux(ip,imo) = rhoaux(ip,imo) + mol(m)%occ(imo)*phi(imo, 1)*phi(imo, 1)
                  enddo
               enddo ! i = 0, n-1

               deallocate (chi, phi, maxc, ldopri)

            enddo ! m = 1, nmol

            ! accumulate intermediate variables
            do i = 0, n(1)-1
               ip=i+1
                do imo = 1, mol(1)%nmo  
                ! rho and grad
                   rhoaux(ip,imo) = max(rhoaux(ip,imo), 1d-30)
               !if (istat /= 0) call error('calcprops', 'Error diagonalizing hessian.', faterr)
               rho_orb(i, j, k, imo) = rhoaux(ip,imo)
            enddo
           enddo
         enddo ! k = 0, n(3)-1
      enddo ! j = 0, n(2)-1

      deallocate (dx, dy, dz, d2, hess, tp)

end subroutine orb_density

subroutine interorb(m,xinit,xinc,nstep,crho,cgrad,rhocut,dimplot,csab,interaction_index,checkpoint)
        use reader
        use props 

        type(molecule),intent(in) ::m
        real*8,intent(in)         ::xinit(3), xinc(3)
        integer,intent(in)        ::nstep(3)
        real*8,intent(in)         ::crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
        real*8,intent(in)         ::cgrad(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
        real*8,intent(in)         ::rhocut, dimplot
        real*8,allocatable,intent(out):: interaction_index(:,:), csab(:,:,:,:)
        integer,intent(out)      ::checkpoint 

        integer                   :: l, p, i, j, k, check
        real*8                    ::csab_temp(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
        real*8                    :: maximum, minimum
        real*8                    :: x(3)
        real*8,allocatable        :: csab1(:,:,:,:), csab2(:,:,:,:) 
        real*8,allocatable        :: interaction_index1(:,:), interaction_index2(:,:)



        csab_temp=0.0
        check=-1
        do l=1, m%nmo-1
           do p=l+1, m%nmo
              csab_temp=0.0
              maximum=-10.d0
              minimum=0.d0
              do k=0, nstep(3)-1
                 do j=0, nstep(2)-1
                    do i=0, nstep(1)-1
                       x=xinit + (/i,j,k/)*xinc 
                       if (abs(crho(i,j,k))/100.d0<rhocut .and. cgrad(i,j,k)<dimplot) then
                               call calcsab(x,m,crho(i,j,k),l,p,csab_temp(i,j,k))
                               minimum=min(csab_temp(i,j,k),minimum)
                       else
                               csab_temp(i,j,k)=0.d0
                       endif
                    enddo
                 enddo
              enddo
              if (minimum < 0.0 ) then 
                      check=check+1
                      if (.not. allocated(csab1)) then
                              allocate(csab1(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1,0:check))
                              allocate(interaction_index1(0:check,0:3))
                              if (check/=0) then
                                      csab1(:,:,:,0:check-1)=csab2
                                      csab1(:,:,:,check)=csab_temp
                                      interaction_index1(0:check-1,:)=interaction_index2
                                      interaction_index1(check,0)=check
                                      interaction_index1(check,1)=l
                                      interaction_index1(check,2)=p
                                      interaction_index1(check,3)=minimum
                                      deallocate(csab2)
                                      deallocate(interaction_index2)
                              else
                                      csab1(:,:,:,check)=csab_temp
                                      interaction_index1(check,0)=check
                                      interaction_index1(check,1)=l
                                      interaction_index1(check,2)=p
                                      interaction_index1(check,3)=minimum
                              endif
                      else if (.not. allocated(csab2)) then
                              allocate(csab2(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1,0:check))
                              allocate(interaction_index2(0:check,0:3))
                              if (check/=0) then
                                      csab2(:,:,:,0:check-1)=csab1
                                      csab2(:,:,:,check)=csab_temp
                                      interaction_index2(0:check-1,:)=interaction_index1
                                      interaction_index2(check,0)=check
                                      interaction_index2(check,1)=l
                                      interaction_index2(check,2)=p
                                      interaction_index2(check,3)=minimum
                                      deallocate(csab1)
                                      deallocate(interaction_index1)
                              else
                                      csab2(:,:,:,check)=csab_temp
                                      interaction_index2(check,0)=check
                                      interaction_index2(check,1)=l
                                      interaction_index2(check,2)=p
                                      interaction_index2(check,3)=minimum
                              endif
                      else
                              call error('nciplot','problem during dynamics allocation',faterr)
                      endif 
              endif
          enddo
        enddo
        print*, 'CHECK: ', check
        allocate(interaction_index(0:check,0:3))
        allocate(csab(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1,0:check))
        if (allocated(interaction_index1)) then
                csab=csab1
                interaction_index=interaction_index1
        else if (allocated(interaction_index2)) then
                csab=csab2
                interaction_index=interaction_index2
        else
                call error('nciplot','impossible to allocate interaction indexes',faterr)
        endif
        checkpoint=check

end subroutine interorb

subroutine write_interorb(rhocut,dimplot,oname, nstep, crho, cgrad,csab, checkpoint) 

        real*8,intent(in)           :: rhocut, dimplot
        integer,intent(in)          :: nstep(3), checkpoint
        real*8,intent(in)           ::crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
        real*8,intent(in)           ::cgrad(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
        real*8,intent(in)           ::csab(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1,0:checkpoint)
        character*(100),intent(in)  :: oname
        integer                     ::i,j,k

        open(unit=23, file=trim(oname)//"s_ab.dat")
        do k=0,nstep(3)-1
           do j=0, nstep(2)-1
              do i=0, nstep(1)-1
                 if (abs(crho(i,j,k))/100.d0<rhocut .and. cgrad(i,j,k)<dimplot) then 
                 write(23,*) crho(i,j,k)/100.d0, csab(i,j,k,:)
                 endif
              enddo
           enddo
        enddo
        close(23)

end subroutine write_interorb

subroutine sort_index(interaction_index,checkpoint,maxvalue)

        real*8,intent(inout)        ::interaction_index(0:checkpoint,0:3)
        integer,intent(in)          :: checkpoint
        real*8,intent(out)          ::maxvalue

        real*8                      ::inter_sort(0:checkpoint,0:3)
        real*8                      ::maximumvalues(0:checkpoint)
        integer                     ::minpos, i

        maximumvalues(:)=interaction_index(0:checkpoint,3)
        do i=0, checkpoint
           minpos=minloc(maximumvalues,1)
           inter_sort(i,:)=interaction_index(minpos-1,:)
           maximumvalues(minpos-1)=100000.d0
        enddo
        interaction_index=inter_sort 
        print*, '-----------------------------------------------------------------------------------------'
        print*, 'Orbital N°     with   Orital N°         Minimum'
        do i=0, checkpoint
           print*, int(interaction_index(i,1)), int(interaction_index(i,2)), interaction_index(i,3)
        enddo
        print*, '-----------------------------------------------------------------------------------------'
        print*, 'What minimum value cutoff do you want to use for the rest of the calculation?'
        !read(*,*) maxvalue
        maxvalue=0.0
        print* , maxvalue

end subroutine sort_index

subroutine write_index(nstep,csab,interaction_index,checkpoint,maxvalue, ntotal, xinit,xinc,m,nfile)

        integer, intent(in)         :: nstep(3),nfile
        real*8, intent(in)          :: interaction_index(0:checkpoint,0:3)
        real*8,intent(in)           :: csab(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1,0:checkpoint)
        real*8,intent(in)           :: maxvalue
        integer, intent(in)         :: checkpoint, ntotal
        type(molecule),intent(in)   :: m 
        real*8, intent(in)          ::xinc(3), xinit(3) 

        character(len=12)           :: interindex
        integer                     :: i,j,k,l
        real*8,allocatable          ::cs2local(:,:,:)

        open(unit=6000, file='interaction_index.dat')
        do l=0, checkpoint
           if(interaction_index(l,3)<maxvalue) then
                   allocate(cs2local(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1))
                   do k=0, nstep(3)-1
                      do j=0, nstep(2)-1
                         do i=0, nstep(1)-1
                            cs2local(i,j,k)=csab(i,j,k,int(interaction_index(l,0)))
                         enddo
                      enddo
                    enddo
           endif
           if (allocated(cs2local)) then
                   
                   write(interindex,'(F10.0)') interaction_index(l,0)
                   open(unit=102, file=trim(adjustl(interindex))//"_interaction_S2.cube")
                   call write_cube_header_onci(102,'S2_cube','3D plot pairwise orbital RDG',ntotal, xinit, xinc, nstep, m, nfile)
                   call write_cube_body_onci(102,nstep,cs2local)
                   close(102)
                   write(6000,*) 'Interaction N°',int(interaction_index(l,0))+1,': OM ',int(interaction_index(l,1)), 'with OM', &
                           int(interaction_index(l,2)),  'minimum value', interaction_index(l,3)
                   deallocate(cs2local)
           endif
        enddo

end subroutine write_index
   
subroutine dataGeom_sab(sum_sab_vol, xinc, nstep, checkpoint, csab,rmbox_coarse, nfiles)
      real*8, intent(inout) :: sum_sab_vol(0:2,0:checkpoint)
      real*8, intent(in) :: xinc(3)
      integer, intent(in) :: nstep(3), checkpoint
      real*8, intent(in) :: csab(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 0:checkpoint)
      logical, intent(in) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      integer :: i, j, k, n, l
      integer, intent(in) :: nfiles
      integer :: i1(2), j1(2), k1(2), negative, positive
      real*8 :: sum_sab
      
      sum_sab_vol=0
      negative = 0
      positive = 0
      ! integral of rho^n over the volume of cubes
      do l=0, checkpoint
      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               if (.not. rmbox_coarse(i, j, k)) then
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
                  sum_sab_vol(0,l)=sum_sab_vol(0,l)+sum(abs(csab(i1,j1,k1,l)))*xinc(1)*xinc(2)*xinc(3)/8
                  sum_sab_vol(1,l) = sum_sab_vol(1,l) + sum(abs(csab(i1, j1, k1,l))**2) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_sab_vol(2,l) = sum_sab_vol(2,l) + sum(abs(csab(i1, j1, k1,l))**3) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
               end if
            end do
         end do
      end do
      enddo

end subroutine dataGeom_sab

subroutine dataGeom_orb_inte(sum_orb_vol, xinc, nstep, mol, crho, rmbox_coarse, nfiles)

      real*8, intent(in) :: xinc(3)
      integer, intent(in) :: nstep(3)
      type(molecule), intent(in) :: mol(nfiles)
      real*8, intent(in) :: crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, mol(1)%nmo)
      logical, intent(in) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      integer :: i, j, k, n, l
      integer, intent(in) :: nfiles
      integer :: i1(2), j1(2), k1(2), negative, positive
      real*8 :: sum_vol
      real*8 :: rho_temp(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1)
      real*8, intent(inout) :: sum_orb_vol(0:2,mol(1)%nmo)

      
      rho_temp=0
      sum_orb_vol=0
      negative = 0
      positive = 0
      
      ! integral of rho^n over the volume of cubes
      do l=1, mol(1)%nmo
         rho_temp=0
         do i=0, nstep(1)-1 
            do j=0, nstep(2)-1
               do k=0, nstep(3)-1
                     rho_temp(i,j,k)=crho(i,j,k,l)
               enddo
            enddo
         enddo         
       do i = 0, nstep(1) - 2
        do j = 0, nstep(2) - 2
          do k = 0, nstep(3) - 2
            if (.not. rmbox_coarse(i, j, k)) then
               i1 = (/i, i + 1/)
               j1 = (/j, j + 1/)
               k1 = (/k, k + 1/)
              ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
               sum_orb_vol(0,l) = sum_orb_vol(0,l) + sum(abs((rho_temp(i1, j1, k1)))) &
                                 *xinc(1)*xinc(2)*xinc(3)/8.d0
               sum_orb_vol(1,l) = sum_orb_vol(1,l) + sum(abs(((rho_temp(i1, j1, k1))))**2) &
                                 *xinc(1)*xinc(2)*xinc(3)/8.d0
               sum_orb_vol(2,l) = sum_orb_vol(2,l) + sum(abs(((rho_temp(i1, j1, k1))))**3) &
                                 *xinc(1)*xinc(2)*xinc(3)/8.d0
            end if
         end do
       end do
     end do
    enddo

end subroutine dataGeom_orb_inte

subroutine dataGeom_gonci(sum_gonci_vol, orb_dens,interaction_index, m, xinc, nstep, checkpoint, csab,rmbox_coarse, nfiles)
      integer, intent(in) :: nstep(3), checkpoint
      real*8, intent(inout) :: sum_gonci_vol(0:checkpoint)
      type(molecule), intent(in) :: m(nfiles)
      real*8, intent(in) :: interaction_index(0:checkpoint,0:3)
      real*8, intent(in) :: orb_dens(0:nstep(1)-1, 0:nstep(2)-1, 0:nstep(3)-1, m(1)%nmo)
      real*8, intent(in) :: xinc(3) 
      real*8 :: rho_temp(0:nstep(1)-1, 0:nstep(2)-1, 0:nstep(3)-1)
      real*8, intent(in) :: csab(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 0:checkpoint)
      logical, intent(in) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      integer :: i, j, k, n, l
      integer, intent(in) :: nfiles
      integer :: i1(2), j1(2), k1(2), negative, positive

      sum_gonci_vol=0
      negative = 0
      positive = 0
      ! integral of rho^n over the volume of cubes
      do l=0, checkpoint
      rho_temp=0
         do i=0, nstep(1)-1
            do j=0, nstep(2)-1
               do k=0, nstep(3)-1
                     rho_temp(i,j,k)=csab(i,j,k,int(interaction_index(l,0)))*orb_dens(i,j,k,int(interaction_index(l,1)))*&
                             orb_dens(i,j,k,int(interaction_index(l,2)))
               enddo
            enddo
         enddo

      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               if (.not. rmbox_coarse(i, j, k)) then
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
                  sum_gonci_vol(l)=sum_gonci_vol(l)+sum(abs(rho_temp(i1,j1,k1)))*xinc(1)*xinc(2)*xinc(3)/8
               end if
            end do
         end do
      end do
      enddo

end subroutine dataGeom_gonci


end module onci
       


 













