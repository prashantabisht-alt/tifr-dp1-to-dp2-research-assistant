!=====================================================================
! Jerky chiral Active Brownian Particle (jcABP) — Euler-Maruyama
!
! Full EOM (Eq. 3 of Jose & Löwen 2025):
!   λ r⃗''' + m r⃗'' + γ r⃗' = γ v0 n̂(t) + sqrt(2Dγ²) η⃗(t)
!   θ̇ = ω₀ + sqrt(2 D_r) η_r(t)
!
! First-order form (7-variable state):
!   dx/dt  = vx
!   dvx/dt = ax
!   dax/dt = -ax/τ_J - vx/τ_F² + v0 cos(θ)/τ_F² + sqrt(2D)/τ_F² · ξ_x
!   (same for y with sin θ)
!   dθ/dt  = ω₀ + sqrt(2 D_r) · ξ_r
!
! where τ_J = λ/m, τ_F = sqrt(λ/γ), D_r = 1/τ_P, ω₀ = 1/τ_C
!
! Outputs: t/τ_unit  <x>  <y>  <x²+y²>
!
! Compile: gfortran -O3 -o jcABP_sim jcABP_sim.f90
!=====================================================================
program jcABP_simulation
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)

  ! ---------- physical parameters (common) ----------
  real(dp), parameter :: v0 = 1.0_dp
  real(dp), parameter :: m  = 1.0_dp
  real(dp), parameter :: D  = 0.5_dp

  ! ---------- simulation controls ----------
  integer  :: seed
  seed = 24681357
  call sgrnd(seed)

  ! ====================================================================
  ! SECTION 3: No chirality (ω₀ = 0, D_r = 1/τ_P)
  ! Fig 1+2 panel (a): τ_J=0.2, τ_F=0.28, τ_P=1.0
  ! ====================================================================
  print *, "=== Section 3 panel (a): tau_J=0.2, tau_F=0.28 ==="
  call run_one_case("jcABP_sec3a.dat",          &
                    tau_J=0.2_dp, tau_F=0.28_dp, &
                    tau_P=1.0_dp, tau_C=1.0d15,  &
                    Tmax=20.0_dp, dt=0.002_dp,   &
                    Ntraj=20000, save_every=25)

  ! Fig 1+2 panel (b): τ_J=20, τ_F=1.0, τ_P=1.0
  print *, "=== Section 3 panel (b): tau_J=20, tau_F=1.0 ==="
  call run_one_case("jcABP_sec3b.dat",           &
                    tau_J=20.0_dp, tau_F=1.0_dp,  &
                    tau_P=1.0_dp,  tau_C=1.0d15,  &
                    Tmax=200.0_dp, dt=0.02_dp,    &
                    Ntraj=20000,   save_every=50)

  ! ====================================================================
  ! SECTION 4.1: Chirality, D_r = 0
  ! Fig 3(e)/4(a): τ_J=0.2, τ_F=0.2, τ_C=1.0, τ_P → ∞
  ! ====================================================================
  print *, "=== Section 4.1: tau_J=0.2, tau_F=0.2, tau_C=1.0 ==="
  call run_one_case("jcABP_sec41.dat",           &
                    tau_J=0.2_dp,  tau_F=0.2_dp,  &
                    tau_P=1.0d15,  tau_C=1.0_dp,  &
                    Tmax=25.0_dp,  dt=0.002_dp,   &
                    Ntraj=10000,   save_every=25)

  ! ====================================================================
  ! SECTION 4.2: Full model (D_r ≠ 0, ω₀ ≠ 0)
  ! Fig 8(a): τ_J=0.1, τ_F=0.1, τ_P=10, τ_C=0.5, v0=0.05
  ! Fig 8(b): τ_J=5.0, τ_F=1.0, τ_P=10, τ_C=0.5, v0=0.05
  ! (Note: v0=0.05 for Sec 4.2 — hardcoded in subroutine call below)
  ! ====================================================================
  print *, "=== Section 4.2 panel (a): tau_J=0.1, tau_F=0.1 ==="
  call run_one_case_v0("jcABP_sec42a.dat",       &
                    tau_J=0.1_dp,  tau_F=0.1_dp,  &
                    tau_P=10.0_dp, tau_C=0.5_dp,  &
                    v0_in=0.05_dp,                &
                    Tmax=20.0_dp,  dt=0.001_dp,   &
                    Ntraj=20000,   save_every=50)

  print *, "=== Section 4.2 panel (b): tau_J=5.0, tau_F=1.0 ==="
  call run_one_case_v0("jcABP_sec42b.dat",       &
                    tau_J=5.0_dp,  tau_F=1.0_dp,  &
                    tau_P=10.0_dp, tau_C=0.5_dp,  &
                    v0_in=0.05_dp,                &
                    Tmax=200.0_dp, dt=0.02_dp,    &
                    Ntraj=20000,   save_every=25)

  print *, "All cases done."

contains

  ! ================================================================
  ! Run one simulation case (uses module-level v0)
  ! ================================================================
  subroutine run_one_case(outfile, tau_J, tau_F, tau_P, tau_C, &
                          Tmax, dt, Ntraj, save_every)
    implicit none
    character(len=*), intent(in) :: outfile
    real(dp), intent(in) :: tau_J, tau_F, tau_P, tau_C, Tmax, dt
    integer,  intent(in) :: Ntraj, save_every

    call run_one_case_v0(outfile, tau_J, tau_F, tau_P, tau_C, &
                         v0, Tmax, dt, Ntraj, save_every)
  end subroutine run_one_case

  ! ================================================================
  ! Run one simulation case (explicit v0)
  ! ================================================================
  subroutine run_one_case_v0(outfile, tau_J, tau_F, tau_P, tau_C, &
                             v0_in, Tmax, dt, Ntraj, save_every)
    implicit none
    character(len=*), intent(in) :: outfile
    real(dp), intent(in) :: tau_J, tau_F, tau_P, tau_C
    real(dp), intent(in) :: v0_in, Tmax, dt
    integer,  intent(in) :: Ntraj, save_every

    integer  :: Nsteps, Nsave, n, it, uout
    real(dp), allocatable :: sum_x(:), sum_y(:), sum_r2(:)

    Nsteps = nint(Tmax / dt)
    Nsave  = Nsteps / save_every

    allocate(sum_x(0:Nsave));   sum_x  = 0.0_dp
    allocate(sum_y(0:Nsave));   sum_y  = 0.0_dp
    allocate(sum_r2(0:Nsave));  sum_r2 = 0.0_dp

    ! ---------- ensemble loop ----------
    do n = 1, Ntraj
      call one_trajectory(sum_x, sum_y, sum_r2, Nsave, Nsteps, &
                          save_every, dt, tau_J, tau_F, tau_P, tau_C, v0_in)
      if (mod(n, 500) == 0) print *, "  Completed", n, "of", Ntraj
    end do

    ! ---------- normalize ----------
    sum_x  = sum_x  / real(Ntraj, dp)
    sum_y  = sum_y  / real(Ntraj, dp)
    sum_r2 = sum_r2 / real(Ntraj, dp)

    ! ---------- write output ----------
    open(newunit=uout, file=outfile, status="replace", action="write")
    write(uout,*) "# t   <x>   <y>   <x^2+y^2>"
    do it = 0, Nsave
      write(uout,'(4(1X,ES20.12))') it * save_every * dt, &
                                     sum_x(it), sum_y(it), sum_r2(it)
    end do
    close(uout)
    print *, "  Saved: ", trim(outfile)

    deallocate(sum_x, sum_y, sum_r2)
  end subroutine run_one_case_v0

  ! ================================================================
  ! Single trajectory — Euler-Maruyama for the 7-variable jcABP
  ! ================================================================
  subroutine one_trajectory(sum_x, sum_y, sum_r2, Nsave, Nsteps, &
                            save_every, dt, tau_J, tau_F, tau_P, tau_C, v0_in)
    implicit none
    real(dp), intent(inout) :: sum_x(0:Nsave), sum_y(0:Nsave), sum_r2(0:Nsave)
    integer,  intent(in)    :: Nsave, Nsteps, save_every
    real(dp), intent(in)    :: dt, tau_J, tau_F, tau_P, tau_C, v0_in

    ! --------- state ---------
    real(dp) :: x, vx, ax
    real(dp) :: y, vy, ay
    real(dp) :: theta

    ! --------- derived constants ---------
    real(dp) :: inv_tJ, inv_tF2, omega0, D_r
    real(dp) :: noise_a, noise_r, sqrt_dt
    real(dp) :: cos_th, sin_th
    real(dp) :: dax, day
    real(dp) :: gx, gy, gr

    integer  :: it, isave

    ! --------- precompute ---------
    inv_tJ  = 1.0_dp / tau_J
    inv_tF2 = 1.0_dp / (tau_F * tau_F)
    sqrt_dt = sqrt(dt)

    ! omega_0 = 1/tau_C  (0 if tau_C very large)
    if (tau_C > 1.0d10) then
      omega0 = 0.0_dp
    else
      omega0 = 1.0_dp / tau_C
    end if

    ! D_r = 1/tau_P  (0 if tau_P very large)
    if (tau_P > 1.0d10) then
      D_r = 0.0_dp
    else
      D_r = 1.0_dp / tau_P
    end if

    noise_a = sqrt(2.0_dp * D) * inv_tF2   ! prefactor for translational noise
    noise_r = sqrt(2.0_dp * D_r)            ! prefactor for rotational noise

    ! --------- initial conditions (all zero) ---------
    x = 0.0_dp;  vx = 0.0_dp;  ax = 0.0_dp
    y = 0.0_dp;  vy = 0.0_dp;  ay = 0.0_dp
    theta = 0.0_dp

    ! save t=0
    sum_x(0)  = sum_x(0)  + x
    sum_y(0)  = sum_y(0)  + y
    sum_r2(0) = sum_r2(0) + x*x + y*y

    isave = 1

    ! --------- time stepping ---------
    do it = 1, Nsteps

      ! Gaussian random kicks
      call gaussian(gx)
      call gaussian(gy)
      call gaussian(gr)

      cos_th = cos(theta)
      sin_th = sin(theta)

      ! jerk (acceleration) update: dax = deterministic*dt + stochastic*sqrt(dt)
      dax = (-inv_tJ*ax - inv_tF2*vx + v0_in*inv_tF2*cos_th) * dt  &
            + noise_a * sqrt_dt * gx
      day = (-inv_tJ*ay - inv_tF2*vy + v0_in*inv_tF2*sin_th) * dt  &
            + noise_a * sqrt_dt * gy

      ! Euler-Maruyama update (order: position, velocity, acceleration, angle)
      x  = x  + vx * dt
      y  = y  + vy * dt
      vx = vx + ax * dt
      vy = vy + ay * dt
      ax = ax + dax
      ay = ay + day
      theta = theta + omega0 * dt + noise_r * sqrt_dt * gr

      ! accumulate ensemble averages
      if (mod(it, save_every) == 0 .and. isave <= Nsave) then
        sum_x(isave)  = sum_x(isave)  + x
        sum_y(isave)  = sum_y(isave)  + y
        sum_r2(isave) = sum_r2(isave) + x*x + y*y
        isave = isave + 1
      end if

    end do

  end subroutine one_trajectory

  ! ================================================================
  ! Box-Muller Gaussian (polar method) — uses MT19937 grnd()
  ! ================================================================
  subroutine gaussian(z)
    implicit none
    real(dp), intent(out) :: z
    real(dp) :: x1, x2, w
    real(dp), external :: grnd
    w = 2.0_dp
    do while (w > 1.0_dp .or. w == 0.0_dp)
      x1 = 1.0_dp - 2.0_dp * grnd()
      x2 = 1.0_dp - 2.0_dp * grnd()
      w  = x1*x1 + x2*x2
    end do
    z = sqrt(-2.0_dp * log(w) / w) * x1
  end subroutine gaussian

end program jcABP_simulation
include 'mt.f90'
