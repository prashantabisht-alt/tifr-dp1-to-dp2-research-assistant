!=====================================================================
! tcrw_fig1c.f90 — Fig 1(c): MSD(t) for 4 chirality values
!
! Paper: Osat, Meyberg, Metson, Speck,
!        "Topological chiral random walker"
!        arXiv:2602.12020v1 [cond-mat.stat-mech], 12 Feb 2026.
! Panel: Fig. 1(c) — ⟨|r(t)|²⟩ vs t on a log–log plot, showing normal
!        diffusion (MSD ∝ t) at long times for several ω.
!
! Parameters (cited directly from the paper):
!   D_r      = 1.0e-3                    ! Fig 1 caption: "D_r = 10^{-3}"
!   T        = 1_000_000                 ! Fig 1(c) x-axis: t ∈ [10^0, 10^6]
!   ω values ∈ { 0.5, 0.7, 0.9, 1.0 }    ! Fig 1(c) legend (paper, p. 2)
!   N_traj   = 1000                      ! our choice (paper silent);
!                                        ! gives visibly smooth log–log MSD
!   n_tpts   = 100                       ! log-spaced sample points t ∈ [1, T]
!
! Physics we expect:
!   - ballistic at very short t:  ⟨r²⟩ ~ t²   (before any rotation)
!   - crossover: at ω ≈ 1 the walker makes tight CW loops → caged until
!     a noise event (≈ every 1/D_r = 10³ steps) frees it
!   - long-t normal diffusion: ⟨r²⟩ = 4 D(ω) t
!   - at ω = 0.5 the walker is effectively achiral; diffusion kicks in
!     almost immediately and the curve is closest to the dashed MSD ∝ t.
!
! Boundary condition:
!   Same equivalence argument as Fig 1(b): PBC with unwrapped positions
!   equals unbounded for a single walker's MSD. We use tcrw_step_unbounded.
!
! Step rule (same as Fig 1(b), see tcrw_step.f90).
!
! Output :  tcrw_fig1c_msd_wX.X.txt          (one file per ω)
!   Columns:  t     <|r|^2>
!   Header:   metadata + column names
!
! Cost estimate:
!   N_traj × T × n_omega = 1000 × 10^6 × 4 = 4 × 10^9 steps total.
!   ~2 RNG calls per step; Mersenne Twister ~50 ns/call → ~100 ns/step.
!   Expect ~7–15 min total on a modern Mac (depending on -O level).
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig1c.f90 -o tcrw_fig1c
! Run   :  ./tcrw_fig1c
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig1c
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)

   ! ---- Fig 1(c) parameters (single source of truth) ----
   real(dp), parameter :: D_r      = 1.0e-3_dp
   integer,  parameter :: T_steps  = 1000000
   integer,  parameter :: N_traj   = 1000
   integer,  parameter :: n_tpts   = 100
   real(dp), parameter :: omega_list(4) = (/ 0.5_dp, 0.7_dp, 0.9_dp, 1.0_dp /)
   integer,  parameter :: seed     = 20260417

   integer :: i
   integer :: t_pts(n_tpts)

   call sgrnd(seed)
   call build_log_time_grid(t_pts, n_tpts, T_steps)

   print '(A)', '==== TCRW Fig 1(c): MSD(t) ensemble run ===='
   print '(A,ES9.2,A,I10,A,I6,A,I4)',                      &
        '  D_r = ', D_r, '   T = ', T_steps,               &
        '   N_traj = ', N_traj, '   n_tpts = ', n_tpts
   print '(A,I10)', '  seed = ', seed
   print '(A)', ''

   do i = 1, size(omega_list)
      call run_msd(omega_list(i), D_r, T_steps, N_traj, t_pts, n_tpts)
   end do

   print '(A)', ''
   print '(A)', 'Wrote 4 files:  tcrw_fig1c_msd_w{0.5,0.7,0.9,1.0}.txt'
   print '(A)', 'Plot with:      bash run.sh tcrw-plot-fig1c'

contains

   !------------------------------------------------------------------
   ! Build a log-spaced integer time grid from 1 to T_max.
   ! Guards against collisions at small t (e.g. 10^0 = 1 = 10^{0.04} rounded).
   !------------------------------------------------------------------
   subroutine build_log_time_grid(t_pts, n_pts, T_max)
      integer, intent(out) :: t_pts(:)
      integer, intent(in)  :: n_pts, T_max
      integer  :: k, t_tmp
      real(dp) :: log_Tmax, frac

      log_Tmax = log10(real(T_max, dp))
      t_pts(1) = 1
      do k = 2, n_pts
         frac   = real(k - 1, dp) / real(n_pts - 1, dp)
         t_tmp  = int( 10.0_dp ** (log_Tmax * frac) )
         if (t_tmp <= t_pts(k - 1)) t_tmp = t_pts(k - 1) + 1
         if (t_tmp >  T_max)        t_tmp = T_max
         t_pts(k) = t_tmp
      end do
   end subroutine build_log_time_grid


   !------------------------------------------------------------------
   ! Ensemble MSD for a single ω.
   ! Inner loop is "all T steps for one walker, then next walker" —
   ! keeps (x, y, d) hot in registers; accumulate |r|² at checkpoints
   ! into shared sum_r2(k), divide by N at the end.
   !------------------------------------------------------------------
   subroutine run_msd(omega, D_r_in, T_in, N_in, t_pts, n_pts)
      real(dp), intent(in) :: omega, D_r_in
      integer,  intent(in) :: T_in, N_in, n_pts
      integer,  intent(in) :: t_pts(:)

      real(dp) :: grnd
      real(dp) :: sum_r2(n_pts), msd_k
      integer  :: j, it, k, u, x, y, d
      character(len=64) :: fname
      real(dp) :: t0, t1

      sum_r2 = 0.0_dp
      call cpu_time(t0)

      do j = 1, N_in
         ! ---- initialize one walker ----
         x = 0; y = 0
         d = int(4.0_dp * grnd())
         if (d == 4) d = 3                     ! guard grnd() ≈ 1.0

         ! ---- walk T_in steps, sampling at log-spaced checkpoints ----
         k = 1                                  ! pointer into t_pts
         do it = 1, T_in
            call tcrw_step_unbounded(x, y, d, omega, D_r_in)
            if (k <= n_pts) then
               if (it == t_pts(k)) then
                  sum_r2(k) = sum_r2(k) + real(x, dp)**2 + real(y, dp)**2
                  k = k + 1
               end if
            end if
         end do
      end do

      call cpu_time(t1)

      ! ---- write to file ----
      write(fname, '(A,F3.1,A)') 'tcrw_fig1c_msd_w', omega, '.txt'
      open(newunit=u, file=trim(fname), status='replace', action='write')
      write(u, '(A,F5.2,A,ES10.3,A,I10,A,I6)') &
           '# TCRW Fig 1(c) MSD | omega = ', omega, &
           ' | D_r = ', D_r_in,  &
           ' | T = ',   T_in,    &
           ' | N_traj = ', N_in
      write(u, '(A)') '# columns:  t   <|r|^2>'
      do k = 1, n_pts
         msd_k = sum_r2(k) / real(N_in, dp)
         write(u, '(I10, 1X, ES15.6)') t_pts(k), msd_k
      end do
      close(u)

      print '(A,F4.1,A,F7.1,A)', &
           '   ω = ', omega, ': done in ', t1 - t0, ' s'
   end subroutine run_msd


   !------------------------------------------------------------------
   ! Shared walker-step kernels (tcrw_step_unbounded / tcrw_step_obc).
   !------------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig1c

include 'mt.f90'
