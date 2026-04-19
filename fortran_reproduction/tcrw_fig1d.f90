!=====================================================================
! tcrw_fig1d.f90 — Fig 1(d): D(ω) for two D_r values (paper match)
!
! Paper: Osat, Meyberg, Metson, Speck,
!        "Topological chiral random walker"
!        arXiv:2602.12020v1 [cond-mat.stat-mech], 12 Feb 2026.
! Panel: Fig. 1(d) — long-time diffusion constant D(ω). The paper
!        shows three D_r values ∈ {10⁻⁴, 10⁻³, 10⁻²}; here we use
!        the inner two {10⁻³, 10⁻²} to demonstrate the collapse while
!        keeping runtime modest and staying in the regime of cleanest
!        statistics.
!        Caption: "The diffusion coefficient decreases linearly with
!                  chirality ω independent of the value of D_r."
!
! Expected physics
! ----------------
!   - Both curves collapse onto the same linear "tent" on linear-y
!     axis: D(ω) ≈ D_max · (1 - 2|ω - 0.5|).
!   - D_max ≈ (1 - D_r)/4 ≈ 0.25 (regardless of D_r), because at
!     ω = 0.5 the walker is effectively a 4-direction lattice random
!     walker with step probability (1 - D_r) per step.
!   - At ω = 0 or 1 the walker is caged in a tight CCW/CW unit square
!     and escapes only via noise events on the timescale 1/D_r. On
!     linear-y (0, 0.3) axis these residual D values look like zero.
!
! Parameters
! ----------
!   D_r values ∈ {1e-3, 1e-2}                   ! paper Fig 1(d) inner
!                                                 two values
!   T_steps    = 1_000_000                      ! same as Figs 1(b)(c)
!   N_traj     = 500                            ! per (D_r, ω) pair
!   n_tpts     = 50                             ! log-spaced checkpoints
!   ω grid     = {0.0, 0.1, ..., 1.0}           ! 11 points
!   t_fit_min  = 10/D_r   (per D_r)             ! ensures fit window is
!                                               ! past the caging
!                                               ! transient for every D_r
!
! Fit-window logic
! ----------------
!   The caging timescale is τ = 1/D_r. The diffusive regime only starts
!   for t ≫ τ. We therefore scale the fit window with D_r:
!       D_r = 1e-2  →  τ=100  →  t_fit_min = 1e3   (3 decades of fit)
!       D_r = 1e-3  →  τ=1e3  →  t_fit_min = 1e4   (2 decades of fit)
!   Both cases have plenty of diffusive data in T = 10⁶ steps.
!
! Cost
! ----
!   2 × 11 × 500 × 10⁶ = 1.1 × 10¹⁰ RNG-driven steps.
!   Scaling from Fig 1(c) (~15 ns/step): expect 3–7 min on a modern
!   Mac.
!
! Output
! ------
!   tcrw_fig1d_msd_dr{3,2}_w{0.0,...,1.0}.txt   (22 MSD files)
!   tcrw_fig1d_D_dr{3,2}.txt                    (2 summary files —
!       the Fig 1(d) datafiles)
!       columns:  omega   D_fit    D_err    slope    slope_err    D_end
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig1d.f90 -o tcrw_fig1d
! Run   :  ./tcrw_fig1d
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig1d
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)

   ! ---- Fig 1(d) parameters (single source of truth) ----
   real(dp), parameter :: D_r_list(2) = (/ 1.0e-3_dp, 1.0e-2_dp /)
   integer,  parameter :: T_steps     = 1000000
   integer,  parameter :: N_traj      = 500
   integer,  parameter :: n_tpts      = 50
   real(dp), parameter :: omega_list(11) = (/ 0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, &
                                              0.4_dp, 0.5_dp, 0.6_dp, 0.7_dp, &
                                              0.8_dp, 0.9_dp, 1.0_dp /)
   integer,  parameter :: seed = 20260418

   integer :: i, id, u_sum, dr_exp
   integer :: t_pts(n_tpts)
   real(dp) :: msd(n_tpts)
   real(dp) :: D_fit, D_err, slope, slope_err, D_end
   real(dp) :: D_r_cur, t_fit_min
   character(len=64) :: fname_sum

   call sgrnd(seed)
   call build_log_time_grid(t_pts, n_tpts, T_steps)

   print '(A)', '==== TCRW Fig 1(d): D(ω, D_r) — paper match ===='
   print '(A,I10,A,I6,A,I4)',                         &
        '  T = ', T_steps, '   N_traj = ', N_traj,    &
        '   n_tpts = ', n_tpts
   print '(A,I10)', '  seed = ', seed
   print '(A)', ''

   ! ---- outer loop over D_r ----
   do id = 1, size(D_r_list)
      D_r_cur   = D_r_list(id)
      t_fit_min = 10.0_dp / D_r_cur                   ! scales with noise τ
      dr_exp    = nint(-log10(D_r_cur))               ! 4, 3, or 2

      write(fname_sum, '(A,I0,A)') 'tcrw_fig1d_D_dr', dr_exp, '.txt'
      open(newunit=u_sum, file=trim(fname_sum), status='replace', action='write')
      write(u_sum, '(A,ES10.3,A,I10,A,I6)')          &
           '# TCRW Fig 1(d) | D_r = ', D_r_cur,      &
           '  T = ', T_steps, '  N_traj = ', N_traj
      write(u_sum, '(A,ES10.3)') '# fit window: t >= ', t_fit_min
      write(u_sum, '(A)') '# columns:  omega   D_fit   D_err   slope   slope_err   D_end'

      print '(A,ES9.2,A,ES9.2)', &
           '  --- D_r = ', D_r_cur, '   fit window t >= ', t_fit_min
      print '(A)', '     ω       D_fit        D_err        slope    s_err       D_end        cpu[s]'
      print '(A)', '  -------  -----------  -----------  --------  --------  -----------  --------'

      ! ---- inner loop over ω ----
      do i = 1, size(omega_list)
         call run_and_fit(omega_list(i), D_r_cur, dr_exp, T_steps, N_traj, &
                          t_pts, n_tpts, msd, t_fit_min, &
                          D_fit, D_err, slope, slope_err, D_end)
         write(u_sum, '(F6.2, 5(1X,ES13.5))') &
              omega_list(i), D_fit, D_err, slope, slope_err, D_end
      end do

      close(u_sum)
      print '(A,A)', '  wrote summary  -> ', trim(fname_sum)
      print '(A)', ''
   end do

   print '(A)', 'Wrote 2 summary files:  tcrw_fig1d_D_dr{3,2}.txt'
   print '(A)', 'Wrote 22 MSD files:     tcrw_fig1d_msd_dr{3,2}_w{0.0,...,1.0}.txt'
   print '(A)', 'Plot with:              bash run.sh tcrw-plot-fig1d'

contains

   !------------------------------------------------------------------
   ! Log-spaced integer time grid, same logic as Fig 1(c).
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
   ! Run ensemble MSD for one (D_r, ω), write per-pair MSD file,
   ! fit diffusive tail.
   !------------------------------------------------------------------
   subroutine run_and_fit(omega, D_r_in, dr_exp, T_in, N_in, t_pts, n_pts, msd, &
                          t_fit_min, D_fit, D_err, slope, slope_err, D_end)
      real(dp), intent(in)  :: omega, D_r_in, t_fit_min
      integer,  intent(in)  :: dr_exp, T_in, N_in, n_pts
      integer,  intent(in)  :: t_pts(:)
      real(dp), intent(out) :: msd(:)
      real(dp), intent(out) :: D_fit, D_err, slope, slope_err, D_end

      real(dp) :: grnd
      real(dp) :: sum_r2(n_pts)
      integer  :: j, it, k, u, x, y, d
      character(len=64) :: fname
      real(dp) :: t0, t1

      sum_r2 = 0.0_dp
      call cpu_time(t0)

      do j = 1, N_in
         x = 0; y = 0
         d = int(4.0_dp * grnd())
         if (d == 4) d = 3

         k = 1
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

      ! ---- MSD ----
      do k = 1, n_pts
         msd(k) = sum_r2(k) / real(N_in, dp)
      end do

      ! ---- write per-(D_r, ω) MSD file ----
      write(fname, '(A,I0,A,F3.1,A)') &
           'tcrw_fig1d_msd_dr', dr_exp, '_w', omega, '.txt'
      open(newunit=u, file=trim(fname), status='replace', action='write')
      write(u, '(A,F5.2,A,ES10.3,A,I10,A,I6)') &
           '# TCRW Fig 1(d) MSD | omega = ', omega, &
           ' | D_r = ', D_r_in,  &
           ' | T = ',   T_in,    &
           ' | N_traj = ', N_in
      write(u, '(A)') '# columns:  t   <|r|^2>'
      do k = 1, n_pts
         write(u, '(I10, 1X, ES15.6)') t_pts(k), msd(k)
      end do
      close(u)

      ! ---- fit log10(msd) vs log10(t) over t >= t_fit_min ----
      call loglog_fit(t_pts, msd, n_pts, t_fit_min, slope, slope_err, D_fit, D_err)

      ! ---- endpoint (slope-1) estimate ----
      D_end = msd(n_pts) / (4.0_dp * real(t_pts(n_pts), dp))

      print '(F6.2, 5(1X,ES11.3), 2X,F7.1)', &
           omega, D_fit, D_err, slope, slope_err, D_end, t1 - t0
   end subroutine run_and_fit


   !------------------------------------------------------------------
   ! OLS regression in log–log space (y = a + b x; y = log10 msd,
   !                                                x = log10 t)
   ! Over t >= t_fit_min.
   !------------------------------------------------------------------
   subroutine loglog_fit(t_pts, msd, n_pts, t_fit_min, slope, slope_err, D_fit, D_err)
      integer,  intent(in)  :: t_pts(:), n_pts
      real(dp), intent(in)  :: msd(:), t_fit_min
      real(dp), intent(out) :: slope, slope_err, D_fit, D_err

      integer  :: k, m
      real(dp) :: x, y
      real(dp) :: sx, sy, sxx, sxy, denom
      real(dp) :: a_int, resid_sum, s2, xbar, ybar

      m = 0
      sx = 0.0_dp; sy = 0.0_dp; sxx = 0.0_dp; sxy = 0.0_dp

      do k = 1, n_pts
         if (real(t_pts(k), dp) < t_fit_min)  cycle
         if (msd(k) <= 0.0_dp)                cycle
         x = log10(real(t_pts(k), dp))
         y = log10(msd(k))
         m = m + 1
         sx  = sx  + x
         sy  = sy  + y
         sxx = sxx + x*x
         sxy = sxy + x*y
      end do

      if (m < 3) then
         slope = -1.0_dp; slope_err = -1.0_dp
         D_fit = -1.0_dp; D_err     = -1.0_dp
         return
      end if

      xbar  = sx / real(m, dp)
      ybar  = sy / real(m, dp)
      denom = sxx - real(m, dp) * xbar * xbar
      slope = (sxy - real(m, dp) * xbar * ybar) / denom
      a_int = ybar - slope * xbar

      resid_sum = 0.0_dp
      do k = 1, n_pts
         if (real(t_pts(k), dp) < t_fit_min)  cycle
         if (msd(k) <= 0.0_dp)                cycle
         x = log10(real(t_pts(k), dp))
         y = log10(msd(k))
         resid_sum = resid_sum + (y - a_int - slope * x)**2
      end do
      s2        = resid_sum / real(m - 2, dp)
      slope_err = sqrt(s2 / denom)

      D_fit = 10.0_dp**a_int / 4.0_dp
      D_err = D_fit * log(10.0_dp) * sqrt( s2 * (1.0_dp/real(m,dp) + xbar*xbar/denom) )
   end subroutine loglog_fit


   !------------------------------------------------------------------
   ! Shared walker-step kernels (tcrw_step_unbounded / tcrw_step_obc).
   !------------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig1d

include 'mt.f90'
