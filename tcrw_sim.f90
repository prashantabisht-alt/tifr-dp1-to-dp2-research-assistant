!=====================================================================
! TCRW — Topological Chiral Random Walker on 2D square lattice
!
! Model rules (Osat, Meyberg, Metson & Speck, arXiv:2602.12020):
!   State: (x, y, d) where d in {0,1,2,3} = {up, right, down, left}
!
!   At each discrete time step:
!     With prob D_r:  NOISE STEP
!       - Walker stays at (x,y)
!       - Director rotates: CCW (d-1 mod 4) with prob omega
!                           CW  (d+1 mod 4) with prob (1-omega)
!
!     With prob (1-D_r):  CHIRAL STEP
!       - Walker translates one step in direction d
!       - Director rotates: CW  (d+1 mod 4) with prob omega
!                           CCW (d-1 mod 4) with prob (1-omega)
!       (rotation chirality OPPOSITE to noise step)
!
!       OBC: if translation would leave grid, entire step is blocked.
!
!   Direction encoding: d=0:up(y+1), 1:right(x+1), 2:down(y-1), 3:left(x-1)
!   Convention (matches tcrw_core.py):
!     With encoding above, (d-1) mod 4 = ↑→←→↓→→→↑ (CCW, angle increases)
!                          (d+1) mod 4 = ↑→→→↓→←→↑ (CW,  angle decreases)
!
! Outputs PBC: t  MSD(t)
! Outputs OBC: visit histogram P(x,y) and current fields J(x,y)
!
! Compile: gfortran -O3 -o tcrw_sim tcrw_sim.f90 mt.f90
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_simulation
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)

  ! Direction vectors: DX(0:3), DY(0:3) for d = 0(up),1(right),2(down),3(left)
  integer, parameter :: DX(0:3) = (/ 0, 1, 0, -1 /)
  integer, parameter :: DY(0:3) = (/ 1, 0, -1, 0 /)

  integer :: seed
  seed = 24681357
  call sgrnd(seed)

  ! ====================================================================
  ! Case 1: PBC — MSD for achiral walker (omega=0.5)
  ! ====================================================================
  print *, "=== PBC: omega=0.5, D_r=0.001 ==="
  call run_pbc("tcrw_pbc_achiral.dat", &
               omega=0.5_dp, D_r=0.001_dp, L=200, &
               T_steps=200000, N_traj=200, save_every=100)

  ! ====================================================================
  ! Case 2: PBC — MSD for chiral walker (omega=0.9)
  ! ====================================================================
  print *, "=== PBC: omega=0.9, D_r=0.1 ==="
  call run_pbc("tcrw_pbc_chiral.dat", &
               omega=0.9_dp, D_r=0.1_dp, L=200, &
               T_steps=200000, N_traj=200, save_every=100)

  ! ====================================================================
  ! Case 3: PBC — D scan at fixed omega for D(omega) plot
  ! ====================================================================
  print *, "=== PBC: D vs omega scan ==="
  call run_D_vs_omega("tcrw_D_vs_omega.dat", &
                      D_r=0.1_dp, L=200, &
                      T_steps=500000, N_traj=500, save_every=100, &
                      omega_min=0.0_dp, omega_max=1.0_dp, N_omega=21)

  ! ====================================================================
  ! Case 4: OBC — steady-state visit histogram & currents
  ! ====================================================================
  print *, "=== OBC: omega=0.9, D_r=0.1, L=20 ==="
  call run_obc("tcrw_obc_visits.dat", "tcrw_obc_currents.dat", &
               omega=0.9_dp, D_r=0.1_dp, L=20, &
               T_steps=1000000, N_traj=100)

  ! ====================================================================
  ! Case 5: PBC — single trajectory for plotting
  ! ====================================================================
  print *, "=== PBC: single trajectory, omega=0.9 ==="
  call run_single_traj("tcrw_trajectory.dat", &
                       omega=0.9_dp, D_r=0.05_dp, L=1000, &
                       T_steps=10000)

  print *, "All cases done."

contains

  ! ================================================================
  ! PBC simulation — compute MSD averaged over N_traj walkers
  ! ================================================================
  subroutine run_pbc(outfile, omega, D_r, L, T_steps, N_traj, save_every)
    implicit none
    character(len=*), intent(in) :: outfile
    real(dp), intent(in) :: omega, D_r
    integer,  intent(in) :: L, T_steps, N_traj, save_every

    integer :: Nsave, n, t, isave, uout, dir
    integer :: x, y, ux, uy
    real(dp), allocatable :: sum_r2(:)
    real(dp) :: r_step, r_rot, r2
    real(dp), external :: grnd

    Nsave = T_steps / save_every
    allocate(sum_r2(0:Nsave))
    sum_r2 = 0.0_dp

    ! ---------- ensemble loop ----------
    do n = 1, N_traj

      ! Initial state: start at (L/2, L/2), random director
      x = L / 2
      y = L / 2
      dir = mod(int(grnd() * 4.0_dp), 4)
      ux = 0   ! unwrapped x displacement
      uy = 0   ! unwrapped y displacement

      ! save t=0
      sum_r2(0) = sum_r2(0) + 0.0_dp

      isave = 1

      ! ---------- time loop ----------
      do t = 1, T_steps
        r_step = grnd()
        r_rot  = grnd()

        if (r_step < D_r) then
          ! --- Noise step: stay put, rotate ---
          ! Convention matches Python: prob omega -> CCW (d-1)
          if (r_rot < omega) then
            dir = mod(dir + 3, 4)       ! CCW (d-1 mod 4 = d+3 mod 4)
          else
            dir = mod(dir + 1, 4)       ! CW  (d+1 mod 4)
          end if
        else
          ! --- Chiral step: translate then rotate ---
          ! Translate
          x  = mod(x + DX(dir) + L, L)
          y  = mod(y + DY(dir) + L, L)
          ux = ux + DX(dir)
          uy = uy + DY(dir)

          ! Rotate OPPOSITE to noise: prob omega -> CW (d+1)
          if (r_rot < omega) then
            dir = mod(dir + 1, 4)       ! CW
          else
            dir = mod(dir + 3, 4)       ! CCW
          end if
        end if

        ! Record MSD
        if (mod(t, save_every) == 0 .and. isave <= Nsave) then
          r2 = real(ux, dp)**2 + real(uy, dp)**2
          sum_r2(isave) = sum_r2(isave) + r2
          isave = isave + 1
        end if
      end do

      if (mod(n, 50) == 0) print *, "  Completed", n, "of", N_traj
    end do

    ! ---------- normalize & write ----------
    sum_r2 = sum_r2 / real(N_traj, dp)

    open(newunit=uout, file=outfile, status="replace", action="write")
    write(uout, '(A)') "# t   MSD"
    do isave = 0, Nsave
      write(uout, '(I12, 1X, ES20.12)') isave * save_every, sum_r2(isave)
    end do
    close(uout)
    print *, "  Saved: ", trim(outfile)

    deallocate(sum_r2)
  end subroutine run_pbc

  ! ================================================================
  ! OBC simulation — visit histogram & current fields
  ! ================================================================
  subroutine run_obc(visit_file, current_file, omega, D_r, L, T_steps, N_traj)
    implicit none
    character(len=*), intent(in) :: visit_file, current_file
    real(dp), intent(in) :: omega, D_r
    integer,  intent(in) :: L, T_steps, N_traj

    integer :: n, t, dir, x, y, nx, ny, uout
    integer(8), allocatable :: visits(:,:)
    real(dp),   allocatable :: Jx(:,:), Jy(:,:)
    real(dp) :: r_step, r_rot, total_steps
    real(dp), external :: grnd
    integer :: ix, iy

    allocate(visits(0:L-1, 0:L-1))
    allocate(Jx(0:L-1, 0:L-1))
    allocate(Jy(0:L-1, 0:L-1))
    visits = 0
    Jx = 0.0_dp
    Jy = 0.0_dp

    ! ---------- ensemble loop ----------
    do n = 1, N_traj

      ! Initial state: random position, random director
      x   = mod(int(grnd() * real(L, dp)), L)
      y   = mod(int(grnd() * real(L, dp)), L)
      dir = mod(int(grnd() * 4.0_dp), 4)

      do t = 1, T_steps
        ! Count visit before step
        visits(x, y) = visits(x, y) + 1

        r_step = grnd()
        r_rot  = grnd()

        if (r_step < D_r) then
          ! --- Noise step ---  prob omega -> CCW (d-1 = d+3)
          if (r_rot < omega) then
            dir = mod(dir + 3, 4)       ! CCW
          else
            dir = mod(dir + 1, 4)       ! CW
          end if
        else
          ! --- Chiral step with OBC ---
          nx = x + DX(dir)
          ny = y + DY(dir)

          if (nx >= 0 .and. nx < L .and. ny >= 0 .and. ny < L) then
            ! Can move: record current, translate, rotate (opposite chirality)
            Jx(x, y) = Jx(x, y) + real(DX(dir), dp)
            Jy(x, y) = Jy(x, y) + real(DY(dir), dp)
            x = nx
            y = ny
            if (r_rot < omega) then
              dir = mod(dir + 1, 4)     ! CW
            else
              dir = mod(dir + 3, 4)     ! CCW
            end if
          end if
          ! Blocked: do nothing (no move, no rotation)
        end if
      end do

      ! Count final position
      visits(x, y) = visits(x, y) + 1

      if (mod(n, 20) == 0) print *, "  OBC traj", n, "of", N_traj
    end do

    total_steps = real(T_steps, dp) * real(N_traj, dp)

    ! ---------- write visit histogram ----------
    open(newunit=uout, file=visit_file, status="replace", action="write")
    write(uout, '(A)') "# x  y  P(x,y)"
    do ix = 0, L-1
      do iy = 0, L-1
        write(uout, '(2I6, 1X, ES20.12)') ix, iy, &
              real(visits(ix, iy), dp) / (total_steps + real(N_traj, dp))
      end do
      write(uout, *)   ! blank line for gnuplot
    end do
    close(uout)
    print *, "  Saved: ", trim(visit_file)

    ! ---------- write current field ----------
    open(newunit=uout, file=current_file, status="replace", action="write")
    write(uout, '(A)') "# x  y  Jx  Jy"
    do ix = 0, L-1
      do iy = 0, L-1
        write(uout, '(2I6, 2(1X,ES20.12))') ix, iy, &
              Jx(ix, iy) / total_steps, Jy(ix, iy) / total_steps
      end do
      write(uout, *)
    end do
    close(uout)
    print *, "  Saved: ", trim(current_file)

    deallocate(visits, Jx, Jy)
  end subroutine run_obc

  ! ================================================================
  ! Measure D(omega) — scan omega at fixed D_r
  ! ================================================================
  subroutine run_D_vs_omega(outfile, D_r, L, T_steps, N_traj, save_every, &
                            omega_min, omega_max, N_omega)
    implicit none
    character(len=*), intent(in) :: outfile
    real(dp), intent(in) :: D_r, omega_min, omega_max
    integer,  intent(in) :: L, T_steps, N_traj, save_every, N_omega

    integer  :: iw, n, t, isave, Nsave, dir, x, y, ux, uy, uout
    integer  :: i0, npts
    real(dp) :: omega, r_step, r_rot, r2, slope, D_eff
    real(dp), allocatable :: sum_r2(:)
    real(dp), allocatable :: t_fit(:), msd_fit(:)
    real(dp) :: sum_t, sum_msd, sum_t2, sum_t_msd
    real(dp), external :: grnd

    Nsave = T_steps / save_every
    allocate(sum_r2(0:Nsave))

    open(newunit=uout, file=outfile, status="replace", action="write")
    write(uout, '(A)') "# omega   D_eff"

    do iw = 0, N_omega - 1
      if (N_omega > 1) then
        omega = omega_min + (omega_max - omega_min) * real(iw, dp) / real(N_omega - 1, dp)
      else
        omega = omega_min
      end if

      sum_r2 = 0.0_dp

      ! --- ensemble ---
      do n = 1, N_traj
        x = L / 2;  y = L / 2
        dir = mod(int(grnd() * 4.0_dp), 4)
        ux = 0;  uy = 0
        isave = 1

        do t = 1, T_steps
          r_step = grnd()
          r_rot  = grnd()

          if (r_step < D_r) then
            ! Noise: prob omega -> CCW (d-1)
            if (r_rot < omega) then
              dir = mod(dir + 3, 4)     ! CCW
            else
              dir = mod(dir + 1, 4)     ! CW
            end if
          else
            x  = mod(x + DX(dir) + L, L)
            y  = mod(y + DY(dir) + L, L)
            ux = ux + DX(dir)
            uy = uy + DY(dir)
            ! Chiral: prob omega -> CW (d+1)
            if (r_rot < omega) then
              dir = mod(dir + 1, 4)     ! CW
            else
              dir = mod(dir + 3, 4)     ! CCW
            end if
          end if

          if (mod(t, save_every) == 0 .and. isave <= Nsave) then
            r2 = real(ux, dp)**2 + real(uy, dp)**2
            sum_r2(isave) = sum_r2(isave) + r2
            isave = isave + 1
          end if
        end do
      end do

      sum_r2 = sum_r2 / real(N_traj, dp)

      ! --- linear fit to last half: MSD = 4D*t ---
      i0 = Nsave / 2
      npts = Nsave - i0 + 1

      ! Least-squares slope: slope = (N*sum(t*msd) - sum(t)*sum(msd)) / (N*sum(t^2) - (sum(t))^2)
      sum_t = 0.0_dp;  sum_msd = 0.0_dp
      sum_t2 = 0.0_dp; sum_t_msd = 0.0_dp
      do isave = i0, Nsave
        r2 = real(isave * save_every, dp)
        sum_t     = sum_t     + r2
        sum_msd   = sum_msd   + sum_r2(isave)
        sum_t2    = sum_t2    + r2 * r2
        sum_t_msd = sum_t_msd + r2 * sum_r2(isave)
      end do
      slope = (real(npts,dp)*sum_t_msd - sum_t*sum_msd) / &
              (real(npts,dp)*sum_t2 - sum_t*sum_t)
      D_eff = slope / 4.0_dp

      write(uout, '(2(1X,ES20.12))') omega, D_eff
      print '(A,F6.3,A,ES12.4)', "  omega=", omega, "  D=", D_eff
    end do

    close(uout)
    print *, "  Saved: ", trim(outfile)

    deallocate(sum_r2)
  end subroutine run_D_vs_omega

  ! ================================================================
  ! Single trajectory output (PBC, unwrapped coords)
  ! ================================================================
  subroutine run_single_traj(outfile, omega, D_r, L, T_steps)
    implicit none
    character(len=*), intent(in) :: outfile
    real(dp), intent(in) :: omega, D_r
    integer,  intent(in) :: L, T_steps

    integer :: t, dir, x, y, ux, uy, uout
    real(dp) :: r_step, r_rot
    real(dp), external :: grnd

    x = L / 2;  y = L / 2
    dir = mod(int(grnd() * 4.0_dp), 4)
    ux = 0;  uy = 0

    open(newunit=uout, file=outfile, status="replace", action="write")
    write(uout, '(A)') "# t  ux  uy"
    write(uout, '(3I12)') 0, ux, uy

    do t = 1, T_steps
      r_step = grnd()
      r_rot  = grnd()

      if (r_step < D_r) then
        ! Noise: prob omega -> CCW (d-1)
        if (r_rot < omega) then
          dir = mod(dir + 3, 4)         ! CCW
        else
          dir = mod(dir + 1, 4)         ! CW
        end if
      else
        x  = mod(x + DX(dir) + L, L)
        y  = mod(y + DY(dir) + L, L)
        ux = ux + DX(dir)
        uy = uy + DY(dir)
        ! Chiral: prob omega -> CW (d+1)
        if (r_rot < omega) then
          dir = mod(dir + 1, 4)         ! CW
        else
          dir = mod(dir + 3, 4)         ! CCW
        end if
      end if

      write(uout, '(3I12)') t, ux, uy
    end do

    close(uout)
    print *, "  Saved: ", trim(outfile)
  end subroutine run_single_traj

end program tcrw_simulation
include 'mt.f90'
