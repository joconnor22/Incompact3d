!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module tbl

  use decomp_2d
  use variables
  use param

  implicit none

  integer :: fs
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_tbl, boundary_conditions_tbl, postprocess_tbl, visu_tbl, visu_tbl_init, cuspline

contains

  subroutine init_tbl (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d_io
    use param , only : zptwofive
    use MPI


    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ierror,ii,is,it,code

    integer, dimension (:), allocatable :: seed

    if (iscalar==1) then

       phi1(:,:,:,:) = zptwofive !change as much as you want
          if ((nclyS1==2).and.(xstart(2)==1)) then
             !! Generate a hot patch on bottom boundary
             phi1(:,1,:,:) = one
          endif
          if ((nclySn==2).and.(xend(2)==ny)) THEN
             phi1(:,xsize(2),:,:) = zptwofive
          endif

    endif
    ux1=zero;uy1=zero;uz1=zero

    !a blasius profile is created in ecoule and then duplicated for the all domain
    call blasius()

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_tbl
  !********************************************************************
  subroutine boundary_conditions_tbl (ux,uy,uz,phi)

    use navier, only : tbl_flrt
    use param , only : zero, zptwofive

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype) :: x, y, z, um
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx

    integer :: i, j, k, is

    !INFLOW with an update of bxx1, byy1 and bzz1 at the inlet

    call blasius()
    !INLET FOR SCALAR, TO BE CONSISTENT WITH INITIAL CONDITION
    if (iscalar==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(1,:,:,:)=zptwofive
             if ((xstart(2)==1)) then
                phi(:,1,:,:) = one
             endif
             if ((xend(2)==ny)) THEN
                phi(:,xsize(2),:,:) = zptwofive
             endif
          enddo
       enddo
    endif

    !OUTFLOW based on a 1D convection equation

    udx=one/dx
    udy=one/dy
    udz=one/dz
    uddx=half/dx
    uddy=half/dy
    uddz=half/dz

    do k=1,xsize(3)
       do j=1,xsize(2)

          cx=ux(nx,j,k)*gdt(itr)*udx

          if (cx<zero) cx=zero
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
          if (iscalar==1) phi(nx,:,:,:) =  phi(nx,:,:,:) - cx*(phi(nx,:,:,:)-phi(nx-1,:,:,:))
          enddo
    enddo

    !! Bottom Boundary
    if (ncly1 == 2) then
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i, k) = zero
          byy1(i, k) = zero
          byz1(i, k) = zero
        enddo
      enddo
      if (iblow == 1) then 
         call wall_blowing()
      end if
    endif
    !! Top Boundary
    if (nclyn == 2) then
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = uy(i, xsize(2) - 1, k)
             byzn(i, k) = uz(i, xsize(2) - 1, k)
          enddo
       enddo
    endif

    !SCALAR   
    if (itimescheme/=7) then
    if (iscalar/=0) then
          if ((nclyS1==2).and.(xstart(2)==1)) then
             !! Generate a hot patch on bottom boundary
             phi(1,1,:,:) = one
          endif
          if ((nclySn==2).and.(xend(2)==ny)) THEN
             phi(1,xsize(2),:,:) = phi(1,xsize(2)-1,:,:)
          endif
    endif
    endif

    !update of the flow rate (what is coming in the domain is getting out)
    call tbl_flrt(ux,uy,uz)

    return
  end subroutine boundary_conditions_tbl

  !********************************************************************
  !********************************************************************
  subroutine blasius()

    use decomp_2d_io
    use MPI
    use param, only : zero, zptwo, zpeight, one, nine
    use dbg_schemes, only: exp_prec, sqrt_prec

    implicit none

    real(mytype) :: eta_bl, f_bl, g_bl, x_bl,h_bl
    real(mytype) :: delta_int, delta_eta, eps_eta


    real(mytype) :: x, y, z
    integer :: i, j, k, is

    do k=1,xsize(3)
       do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)

          eta_bl=y*real(4.91,mytype)/nine

          !OLD POLYNOMIAL FITTING

          delta_eta=zero
          eps_eta=zero
          delta_int=zptwo

          if (eta_bl>=(real(7.5,mytype)/nine)) then
             delta_eta=eta_bl-real(7.5,mytype)/nine
             eta_bl=real(7.5,mytype)/nine
             eps_eta=real(0.00015,mytype)
          end if

          f_bl=1678.64209592595000_mytype*eta_bl**14-11089.69250174290_mytype*eta_bl**13 &
               +31996.435014067000_mytype*eta_bl**12-52671.52497797990_mytype*eta_bl**11 &
               +54176.169116766700_mytype*eta_bl**10-35842.82047060970_mytype*eta_bl**9  &
               +15201.308887124000_mytype*eta_bl**8 -4080.171379356480_mytype*eta_bl**7  &
               +702.12963452810300_mytype*eta_bl**6 -56.20639258053180_mytype*eta_bl**5  &
               -17.018112827391400_mytype*eta_bl**4 +0.819582894357566_mytype*eta_bl**3  &
               -0.0601348202321954_mytype*eta_bl**2 +2.989739912704050_mytype*eta_bl**1

          f_bl=f_bl+(1-exp_prec(-delta_eta/delta_int))*eps_eta



          if (eta_bl >= (7.15_mytype/nine)) then
             delta_int=zpeight
             delta_eta=eta_bl-7.15_mytype/nine
             eta_bl   =       7.15_mytype/nine
             eps_eta  =     0.0005_mytype
          end if

          g_bl=4924.052847797540_mytype*eta_bl**14-34686.2970972733000_mytype*eta_bl**13 &
               +108130.253843618_mytype*eta_bl**12-195823.099139525000_mytype*eta_bl**11 &
               +227305.908339065_mytype*eta_bl**10-176106.001047617000_mytype*eta_bl**9  &
               +92234.5885895112_mytype*eta_bl**8 -32700.3687158807000_mytype*eta_bl**7  &
               +7923.51008739107_mytype*eta_bl**6 -1331.09245288739000_mytype*eta_bl**5  &
               +130.109496961069_mytype*eta_bl**4 -7.64507811014497000_mytype*eta_bl**3  &
               +6.94303207046209_mytype*eta_bl**2 -0.00209716712558639_mytype*eta_bl**1 ! &

          g_bl=g_bl+(1-exp_prec(-delta_eta/delta_int))*eps_eta



          x_bl=one/(4.91_mytype**2*xnu)

          bxx1(j,k)=f_bl/1.0002014996204402_mytype/1.0000000359138641_mytype !To assure 1.0 in infinity
          bxy1(j,k)=g_bl*sqrt_prec(xnu/x_bl)/1.000546554_mytype
          bxz1(j,k)=zero
       enddo
    enddo

    !STORE VALUE F_BL_INF G_BL_INF (ONLY ONE MORE TIME)------------------

    y=yly
    eta_bl=y*4.91_mytype/nine  !The 9 is due to interpolation

    delta_eta=zero
    eps_eta=zero
    delta_int=zptwo

    if (eta_bl>=(7.5_mytype/nine)) then
       delta_eta=eta_bl-7.5_mytype/nine
       eta_bl   =       7.5_mytype/nine
       eps_eta  =   0.00015_mytype
    end if

    !To assure 1.0 in infinity
    f_bl_inf=1678.6420959259500_mytype*eta_bl**14-11089.69250174290_mytype*eta_bl**13 &
            +31996.435014067000_mytype*eta_bl**12-52671.52497797990_mytype*eta_bl**11 &
            +54176.169116766700_mytype*eta_bl**10-35842.82047060970_mytype*eta_bl**9  &
            +15201.308887124000_mytype*eta_bl**8 -4080.171379356480_mytype*eta_bl**7  &
            +702.12963452810300_mytype*eta_bl**6 -56.20639258053180_mytype*eta_bl**5  &
            -17.018112827391400_mytype*eta_bl**4 +0.819582894357566_mytype*eta_bl**3  &
            -0.0601348202321954_mytype*eta_bl**2 +2.989739912704050_mytype*eta_bl**1


    f_bl_inf=f_bl_inf+(1-exp_prec(-delta_eta/delta_int))*eps_eta
    f_bl_inf=f_bl_inf/1.0002014996204402_mytype/1.0000000359138641_mytype !To assure 1.0 in infinity

#ifdef DEBG
    if (nrank == 0) write(*,*)'f_bl_inf ', f_bl_inf
#endif

    if (eta_bl>= (7.15_mytype/nine)) then
       delta_int=zpeight
       delta_eta=eta_bl-7.15_mytype/nine
       eta_bl   =       7.15_mytype/nine
       eps_eta  =     0.0005_mytype
    end if

    g_bl_inf=4924.05284779754_mytype*eta_bl**14-34686.2970972733000_mytype*eta_bl**13 &
            +108130.253843618_mytype*eta_bl**12-195823.099139525000_mytype*eta_bl**11 &
            +227305.908339065_mytype*eta_bl**10-176106.001047617000_mytype*eta_bl**9  &
            +92234.5885895112_mytype*eta_bl**8 -32700.3687158807000_mytype*eta_bl**7  &
            +7923.51008739107_mytype*eta_bl**6 -1331.09245288739000_mytype*eta_bl**5  &
            +130.109496961069_mytype*eta_bl**4 -7.64507811014497000_mytype*eta_bl**3  &
            +6.94303207046209_mytype*eta_bl**2 -0.00209716712558639_mytype*eta_bl**1


    g_bl_inf=g_bl_inf+(1-exp_prec(-delta_eta/delta_int))*eps_eta
    g_bl_inf=g_bl_inf/1.000546554_mytype
#ifdef DEBG
    if (nrank == 0) write(*,*)'g_bl_inf ', g_bl_inf
#endif

    return
  end subroutine blasius

   !############################################################################
   subroutine wall_blowing()

      implicit none

      real(mytype) :: amp_time, amp_space, x
      integer :: i, j, k

      amp_time = exp_blend(t, blow_ramp / two, blow_ramp)

      if (xstart(2) == 1) then
         do i = 1, nx
            x = dx * (i - 1)
            if (x >= blow_x(1) .and. x <= blow_x(size(blow_x))) then
               do j = 1, size(blow_x) - 1
                  if (x >= blow_x(j) .and. x <= blow_x(j+1)) then
                     k = j
                     exit
                  end if
               end do
               amp_space = cusplint(blow_x(k:k+1), blow_amp(k:k+1), blow_ampd2(k:k+1), x)
               byy1(i,:) = amp_time * amp_space
            end if
         enddo
      endif
   end subroutine wall_blowing

   ! Cubic spline (2nd derivative)
   function cuspline(loc, val) result(d2val)

      implicit none

      real(mytype), dimension(:), intent(in) :: loc, val
      real(mytype), dimension(size(loc)) :: d2val
      real(mytype), dimension(size(loc)) :: rhs
      real(mytype) :: dy1, dyn, a, b, c, d
      integer i

      d2val = 0.0_mytype
      rhs = 0.0_mytype

      ! Set 1st derivative BCs
      dy1 = 0.0_mytype
      dyn = 0.0_mytype

      ! Thomas algorithm for tridiagonal systems (with BCs for spline)
      d2val(1) = 0.5_mytype
      rhs(1) = (3.0_mytype / (loc(2) - loc(1))) * ((val(2) - val(1)) / (loc(2) - loc(1)) - dy1);
      do i = 2, size(loc) - 1
          a = loc(i) - loc(i-1)
          b = 2.0_mytype * (loc(i+1) - loc(i-1))
          c = loc(i+1) - loc(i)
          d = 6.0_mytype * ((val(i+1) - val(i)) / c - (val(i) - val(i-1)) / a)
          d2val(i) = c / (b - a * d2val(i-1))
          rhs(i) = (d - a * rhs(i-1)) / (b - a * d2val(i-1))
      end do
      a = loc(size(loc)) - loc(size(loc)-1);
      b = 2.0_mytype * (loc(size(loc)) - loc(size(loc)-1));
      d = 6.0_mytype * (dyn - (val(size(loc)) - val(size(loc)-1)) / (loc(size(loc)) - loc(size(loc)-1)));
      rhs(size(loc)) = (d - a * rhs(size(loc)-1)) / (b - a * d2val(size(loc)-1));
      d2val(size(loc)) = rhs(size(loc));
      do i = size(loc) - 1, 1, -1
          d2val(i) = rhs(i) - d2val(i) * d2val(i+1);
      end do
   end function cuspline

   ! Cubic spline (interpolate)
   function cusplint(loc, val, d2val, intloc) result(intval)

      implicit none

      real(mytype), dimension(2), intent(in) :: loc, val, d2val
      real(mytype), intent(in) :: intloc
      real(mytype) :: intval
      real(mytype) :: A, B, C, D

      A = (loc(2) - intloc) / (loc(2) - loc(1))
      B = (intloc - loc(1)) / (loc(2) - loc(1))
      C = (A**3 - A) * (loc(2) - loc(1))**2 / 6.0_mytype
      D = (B**3 - B) * (loc(2) - loc(1))**2 / 6.0_mytype
      intval = A * val(1) + B * val(2) + C * d2val(1) + D * d2val(2)
   end function cusplint

   !############################################################################
   function quintic_blend(x, centre, width) result(y)

      implicit none

      real(mytype), intent(in) :: x, centre, width
      real(mytype) :: xdash, y

      xdash = (x - centre) / width
      if (xdash < -half) then
         y = zero
      else if (xdash > half) then
         y = one
      else
         y = 0.5_mytype + 1.875_mytype * xdash - 5.0_mytype * xdash**3 + 6.0_mytype * xdash**5
      end if
   end function

   !############################################################################
   function septic_blend(x, centre, width) result(y)

      implicit none

      real(mytype), intent(in) :: x, centre, width
      real(mytype) :: xdash, y

      xdash = (x - centre) / width
      if (xdash < -half) then
         y = zero
      else if (xdash > half) then
         y = one
      else
         y = 0.5_mytype + 2.1875_mytype * xdash - 8.75_mytype * xdash**3 + 21.0_mytype * xdash**5 - 20.0_mytype * xdash**7
      end if
   end function

   !############################################################################
   function exp_blend(x, centre, width) result(y)

      implicit none

      real(mytype), intent(in) :: x, centre, width
      real(mytype) :: xdash, y
      real(mytype), parameter :: TOL = 1.0e-6_mytype

      xdash = (x - centre) / width
      if (xdash <= -half + TOL) then
         y = zero
      else if (xdash >= half - TOL) then
         y = one
      else
         y = one / (one + exp(one / (xdash - 0.5_mytype) + one / (xdash + 0.5_mytype)))
      end if
   end function

  !############################################################################
  subroutine postprocess_tbl(ux1,uy1,uz1,ep1)

    USE MPI
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE ibm_param
    
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    character(len=30) :: filename

  end subroutine postprocess_tbl

  subroutine visu_tbl_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)

    visu_initialised = .true.

  end subroutine visu_tbl_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_tbl
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs TBL-specific visualization
  !!
  !############################################################################
  subroutine visu_tbl(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    character(len=32), intent(in) :: num

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)
    call write_field(di1, ".", "vort", trim(num), flush=.true.) ! Reusing temporary array, force flush

  end subroutine visu_tbl

end module tbl
