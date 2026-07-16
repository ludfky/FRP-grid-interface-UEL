! ============================================================================
! Four-node three-dimensional zero-thickness interface UEL for Abaqus/Standard
!
! Node ordering:
!   1: matrix/concrete side, start of element
!   2: FRP-grid side,       start of element
!   3: FRP-grid side,       end of element
!   4: matrix/concrete side,end of element
!
! Each node has three translational degrees of freedom (1,2,3).
! Two Gauss points are used along the element length. The branch width is
! integrated analytically through the tributary area b_f * l_e.
! ============================================================================

subroutine UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS, &
               PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE, &
               TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP, &
               ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD, &
               PNEWDT,JPROPS,NJPROP,PERIOD)

  include 'ABA_PARAM.INC'

  integer :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE
  integer :: KSTEP, KINC, JELEM, NDLOAD, NPREDF, MLVARX, MDLOAD
  integer :: NJPROP
  integer :: JDLTYP(MDLOAD,*), LFLAGS(*), JPROPS(*)
  real*8  :: RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), SVARS(NSVARS)
  real*8  :: ENERGY(8), PROPS(NPROPS), COORDS(MCRD,NNODE)
  real*8  :: U(NDOFEL), DU(MLVARX,*), V(NDOFEL), A(NDOFEL)
  real*8  :: TIME(2), DTIME, PARAMS(*), ADLMAG(MDLOAD,*)
  real*8  :: PREDEF(2,NPREDF,NNODE), DDLMAG(MDLOAD,*), PNEWDT, PERIOD

  integer, parameter :: ngp = 2
  integer, parameter :: nsv_gp = 9
  integer :: i, j, k, gp, io
  real*8 :: xi(ngp), wt(ngp), n1, n2, bf, le, jacw
  real*8 :: bmat(3,12), dmat(3,3), traction(3)
  real*8 :: delta(3), delta_old(3), uold(12), fint(12)
  real*8 :: eu(3), ev(3), ew(3), xs(3), xe(3)
  real*8 :: state(nsv_gp)

  call zero_vector(RHS(1,1), MLVARX)
  call zero_matrix(AMATRX, NDOFEL, NDOFEL)
  call zero_vector(fint, 12)
  call zero_vector(ENERGY, 8)

  if (NNODE /= 4 .or. NDOFEL /= 12) then
     write(7,*) 'FRP-grid UEL error: expected NNODE=4 and NDOFEL=12.'
     write(7,*) 'JELEM=', JELEM, ' NNODE=', NNODE, ' NDOFEL=', NDOFEL
     PNEWDT = 0.25d0
     return
  end if
  if (NPROPS < 33) then
     write(7,*) 'FRP-grid UEL error: at least 33 real properties are required.'
     write(7,*) 'JELEM=', JELEM, ' NPROPS=', NPROPS
     PNEWDT = 0.25d0
     return
  end if
  if (NSVARS < ngp*nsv_gp) then
     write(7,*) 'FRP-grid UEL error: at least 18 state variables are required.'
     write(7,*) 'JELEM=', JELEM, ' NSVARS=', NSVARS
     PNEWDT = 0.25d0
     return
  end if

  bf = PROPS(1)
  eu = PROPS(25:27)
  ev = PROPS(28:30)
  ew = PROPS(31:33)
  call orthonormalize_frame(eu, ev, ew)

  ! Element length measured between the paired-node midpoints.
  do i = 1, 3
     xs(i) = 0.5d0*(COORDS(i,1) + COORDS(i,2))
     xe(i) = 0.5d0*(COORDS(i,4) + COORDS(i,3))
  end do
  le = norm3(xe-xs)
  if (le <= 1.0d-12) then
     write(7,*) 'FRP-grid UEL error: zero or invalid element length. JELEM=', JELEM
     PNEWDT = 0.25d0
     return
  end if

  do i = 1, 12
     uold(i) = U(i) - DU(i,1)
  end do

  xi(1) = -1.0d0/sqrt(3.0d0)
  xi(2) =  1.0d0/sqrt(3.0d0)
  wt(1) = 1.0d0
  wt(2) = 1.0d0

  do gp = 1, ngp
     n1 = 0.5d0*(1.0d0-xi(gp))
     n2 = 0.5d0*(1.0d0+xi(gp))
     call build_b_matrix(n1, n2, eu, ev, ew, bmat)

     delta     = matmul(bmat, U)
     delta_old = matmul(bmat, uold)

     io = (gp-1)*nsv_gp
     state = SVARS(io+1:io+nsv_gp)
     call interface_update(delta, delta_old, PROPS, state, traction, dmat)
     SVARS(io+1:io+nsv_gp) = state

     jacw = 0.5d0*bf*le*wt(gp)
     fint = fint + matmul(transpose(bmat), traction)*jacw
     AMATRX = AMATRX + matmul(transpose(bmat), matmul(dmat,bmat))*jacw
  end do

  do i = 1, 12
     RHS(i,1) = -fint(i)
  end do

contains

  subroutine interface_update(del, del_old, p, svars_gp, t, d)
    real*8, intent(in)    :: del(3), del_old(3), p(*)
    real*8, intent(inout) :: svars_gp(nsv_gp)
    real*8, intent(out)   :: t(3), d(3,3)

    real*8 :: kuu0, tumax, dup, duf, tures
    real*8 :: kvc, kvt, dv0, kuv0, ac, at
    real*8 :: betaK, betaT, lambdaR, kww, omega
    real*8 :: eta_ku, eta_kv, eta_kuv, eta_tu, eta_dv
    real*8 :: eta_bk, eta_bt
    real*8 :: du, dv, dw, x, xold, sgn, oldmax, tol
    real*8 :: mcycle, chi, residual, kunl, tenv, kenv
    real*8 :: tbase, kbase, tvbase, kvv, kuv, tmax_m
    real*8 :: line_t
    logical :: was_unloading, was_reloading

    kuu0   = p(2)
    tumax  = p(3)
    dup    = p(4)
    duf    = p(5)
    tures  = p(6)
    kvc    = p(7)
    kvt    = p(8)
    dv0    = p(9)
    kuv0   = p(10)
    ac     = p(11)
    at     = p(12)
    betaK  = p(13)
    betaT  = p(14)
    lambdaR= p(15)
    kww    = p(16)
    eta_ku = p(17)
    eta_kv = p(18)
    eta_kuv= p(19)
    eta_tu = p(20)
    eta_dv = p(21)
    eta_bk = p(22)
    eta_bt = p(23)
    omega  = max(0.0d0,min(1.0d0,p(24)))

    ! Equivalent node-region interpolation p^e = p^S[1+w(eta-1)].
    kuu0  = kuu0 *(1.0d0 + omega*(eta_ku -1.0d0))
    kvc   = kvc  *(1.0d0 + omega*(eta_kv -1.0d0))
    kvt   = kvt  *(1.0d0 + omega*(eta_kv -1.0d0))
    kuv0  = kuv0 *(1.0d0 + omega*(eta_kuv-1.0d0))
    tumax = tumax*(1.0d0 + omega*(eta_tu -1.0d0))
    dv0   = dv0  *(1.0d0 + omega*(eta_dv -1.0d0))
    betaK = betaK*(1.0d0 + omega*(eta_bk -1.0d0))
    betaT = betaT*(1.0d0 + omega*(eta_bt -1.0d0))

    du = del(1); dv = del(2); dw = del(3)
    x = abs(du); xold = abs(del_old(1))
    if (du >= 0.0d0) then
       sgn = 1.0d0
    else
       sgn = -1.0d0
    end if
    tol = 1.0d-10

    if (svars_gp(5) < 1.0d0) then
       svars_gp = 0.0d0
       svars_gp(1) = x
       svars_gp(2) = abs(dv)
       svars_gp(5) = 1.0d0
       svars_gp(6) = 0.0d0
       svars_gp(7) = betaK
       svars_gp(8) = betaT
       svars_gp(9) = omega
    end if

    oldmax = max(0.0d0,svars_gp(1))
    mcycle = max(1.0d0,svars_gp(5))
    chi = svars_gp(6)
    was_unloading = (chi < -0.5d0)
    was_reloading = (chi > 0.5d0 .and. chi < 1.5d0)

    if (x < xold-tol) then
       ! Unloading begins or continues.
       if (.not. was_unloading) svars_gp(3) = xold
       chi = -1.0d0
       residual = lambdaR*max(oldmax,xold)
       svars_gp(4) = residual
    else if (x > xold+tol) then
       if (was_unloading) chi = 1.0d0
       if ((chi > 0.5d0 .and. chi < 1.5d0) .and. x > oldmax+1.0d-6) then
          mcycle = mcycle + 1.0d0
          chi = 2.0d0
       end if
       if (x > oldmax) oldmax = x
    end if

    tmax_m = betaT**(mcycle-1.0d0)*tumax
    call tangential_envelope(x, kuu0, tmax_m, dup, duf, tures, tenv, kenv)
    kunl = max(1.0d-12,betaK*kuu0)
    residual = lambdaR*max(oldmax,x)

    if (chi < -0.5d0 .or. (chi > 0.5d0 .and. chi < 1.5d0)) then
       line_t = kunl*max(0.0d0,x-residual)
       tbase = min(line_t,tenv)
       if (line_t <= tenv .and. x > residual) then
          kbase = kunl
       else if (line_t <= tenv) then
          kbase = 0.0d0
       else
          kbase = kenv
       end if
    else
       tbase = tenv
       kbase = kenv
    end if

    if (dv <= 0.0d0) then
       tvbase = kvc*dv
       kvv = kvc
    else if (dv <= dv0 .and. dv0 > 1.0d-14) then
       tvbase = kvt*dv*(1.0d0-dv/dv0)
       kvv = kvt*(1.0d0-2.0d0*dv/dv0)
    else
       tvbase = 0.0d0
       kvv = 0.0d0
    end if

    if (dv0 > 1.0d-14) then
       kuv = kuv0*(1.0d0 + ac*max(0.0d0,-dv/dv0) &
                         - at*max(0.0d0, dv/dv0))
    else
       kuv = kuv0
    end if
    kuv = max(0.0d0,kuv)

    ! Mixed-mode traction update consistent with the symmetric algorithmic matrix.
    t(1) = sgn*tbase + kuv*dv
    t(2) = tvbase + kuv*du
    t(3) = kww*dw

    d = 0.0d0
    d(1,1) = kbase
    d(1,2) = kuv
    d(2,1) = kuv
    d(2,2) = kvv
    d(3,3) = kww

    svars_gp(1) = max(oldmax,x)
    svars_gp(2) = max(svars_gp(2),abs(dv))
    svars_gp(4) = residual
    svars_gp(5) = mcycle
    svars_gp(6) = chi
    svars_gp(7) = betaK
    svars_gp(8) = betaT
    svars_gp(9) = omega
  end subroutine interface_update

  subroutine tangential_envelope(x,k0,tmax,dp,df,tres,tval,kt)
    real*8, intent(in)  :: x,k0,tmax,dp,df,tres
    real*8, intent(out) :: tval,kt
    if (x <= dp) then
       tval = k0*x
       kt = k0
    else if (x <= df .and. df > dp+1.0d-14) then
       kt = -(tmax-tres)/(df-dp)
       tval = tmax + kt*(x-dp)
    else
       tval = tres
       kt = 0.0d0
    end if
  end subroutine tangential_envelope

  subroutine build_b_matrix(n1,n2,eu,ev,ew,b)
    real*8, intent(in)  :: n1,n2,eu(3),ev(3),ew(3)
    real*8, intent(out) :: b(3,12)
    integer :: a, c
    real*8 :: fac(4)
    fac(1)= n1; fac(2)=-n1; fac(3)=-n2; fac(4)= n2
    b = 0.0d0
    do a=1,4
       do c=1,3
          b(1,3*(a-1)+c)=fac(a)*eu(c)
          b(2,3*(a-1)+c)=fac(a)*ev(c)
          b(3,3*(a-1)+c)=fac(a)*ew(c)
       end do
    end do
  end subroutine build_b_matrix

  subroutine orthonormalize_frame(eu,ev,ew)
    real*8, intent(inout) :: eu(3),ev(3),ew(3)
    real*8 :: ne, nv
    ne = norm3(eu)
    if (ne <= 1.0d-14) then
       eu=(/1.0d0,0.0d0,0.0d0/)
    else
       eu=eu/ne
    end if
    ev = ev-dot_product(ev,eu)*eu
    nv = norm3(ev)
    if (nv <= 1.0d-14) then
       if (abs(eu(1)) < 0.9d0) then
          ev=(/1.0d0,0.0d0,0.0d0/)
       else
          ev=(/0.0d0,1.0d0,0.0d0/)
       end if
       ev=ev-dot_product(ev,eu)*eu
       ev=ev/norm3(ev)
    else
       ev=ev/nv
    end if
    call cross3(ev,eu,ew)
    ew=-ew
    ew=ew/norm3(ew)
    call cross3(eu,ew,ev)
    ev=ev/norm3(ev)
  end subroutine orthonormalize_frame

  subroutine cross3(a,b,c)
    real*8, intent(in) :: a(3),b(3)
    real*8, intent(out):: c(3)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
  end subroutine cross3

  real*8 function norm3(a)
    real*8, intent(in) :: a(3)
    norm3=sqrt(max(0.0d0,dot_product(a,a)))
  end function norm3

  subroutine zero_vector(a,n)
    integer, intent(in) :: n
    real*8, intent(out) :: a(n)
    a=0.0d0
  end subroutine zero_vector

  subroutine zero_matrix(a,n,m)
    integer, intent(in) :: n,m
    real*8, intent(out) :: a(n,m)
    a=0.0d0
  end subroutine zero_matrix

end subroutine UEL
