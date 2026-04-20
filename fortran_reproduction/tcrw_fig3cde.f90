!=====================================================================
! tcrw_fig3cde.f90 — Fig 3(c), 3(d), 3(e):  left-wall current vector
!                    fields (c) J_Dr  (d) J_ω and (e) θ_JDr vs D_r.
!
! Paper: Osat, Meyberg, Metson, Speck,
!        "Topological chiral random walker"
!        arXiv:2602.12020v1 [cond-mat.stat-mech], 12 Feb 2026.
!
! Panels covered (one Fortran run feeds all three gnuplot files)
! --------------------------------------------------------------
!   Fig 3(c):  Quiver of J_Dr(y) on the left wall for each D_r.
!              x-axis = D_r (log),  y-axis = edge site index (1..L-2
!              shown; corners excluded in plot).  At each cell draw
!              an arrow  (⟨Δx⟩, ⟨Δy⟩)  with x_before = 0.
!              Paper claim: arrow points DOWN at small D_r (−π/2),
!              rotates toward the bulk (θ → 0) as D_r → 1.
!   Fig 3(d):  Same, but with J_ω instead of J_Dr.
!              Paper claim: arrow CONSTANT direction, θ = +π/4 for
!              all D_r (chiral-continuation current is rigidly set
!              by the CW orbit geometry; its magnitude changes with
!              D_r but not its angle).
!   Fig 3(e):  θ_JDr(D_r) — the angle of the TOTAL left-wall J_Dr
!              vector (summed over all wall sites) as a function of
!              D_r.  Monotone from −π/2 at D_r → 0 to ~0 at D_r → 1.
!
! Observable per (D_r, y) cell
! ----------------------------
!   For each D_r in the grid and each left-wall site y ∈ {0..L-1}:
!       Jx_Dr(y) = Σ over steps (step_type = 1 AND x_before = 0 AND
!                                y_before = y AND prev_noise) of DX(d_before)
!       Jy_Dr(y) = same, with DY(d_before)
!       Jx_om(y), Jy_om(y) = same but for prev_noise = .false.
!
!   We additionally output per-site magnitudes and angles for plotting
!   convenience:
!       |J_Dr|(y)  = sqrt(Jx_Dr^2 + Jy_Dr^2)
!       θ_JDr(y)   = atan2(Jy_Dr, Jx_Dr)          (radians, -π..π)
!       |J_ω|(y), θ_Jω(y) similarly.
!
!   For Fig 3(e) we aggregate:
!       Jx_Dr_tot = Σ_y Jx_Dr(y),     Jy_Dr_tot = Σ_y Jy_Dr(y)
!       θ_JDr_tot = atan2(Jy_Dr_tot, Jx_Dr_tot)
!       θ_Jω_tot  = atan2(Jy_om_tot,  Jx_om_tot)    (for reference in Fig 3j)
!
! Parameters  (paper: L = 10, ω = 1, D_r sweep)
! ---------------------------------------------
!   L             = 10
!   ω             = 1.0
!   D_r grid      = 25 log-spaced pts  10^-4 → 10^0  (match Fig 3a,b)
!   T_floor       = 10^8
!   N_burn_floor  = 10^7
!   K_meas        = 30                         ! reduced vs Fig 3(b)
!   K_burn        = 5
!   seed          = 20260423
!
!   Why K_meas = 30 instead of 100?  Fig 3(c)/(d) are QUIVER plots —
!   we need directions and relative magnitudes to be right, not
!   per-site magnitudes to 1% accuracy.  At small D_r the per-cell
!   budget is 30 · 10^8 = 3 · 10^9 steps; summed over L = 10 wall
!   sites and averaged, that's plenty of SNR for the arrow pattern.
!   Total cost ≈ 15–20 min on a modern Mac.
!
! Current decomposition (same rule as Fig 2, Fig 3(b), Fig 3(g))
! --------------------------------------------------------------
!   tcrw_step_mask returns step_type ∈ {0, 1, 2}.
!   prev_noise ← (step_type == 0).
!   At step_type == 1 AND x_before == 0:
!       prev_noise = .true.  → contribute (DX, DY) to J_Dr(y_before)
!       prev_noise = .false. → contribute (DX, DY) to J_ω (y_before)
!
! Output
! ------
!   tcrw_fig3cde_summary.txt   — per-site currents, columns:
!     # iD  D_r  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om  |J_Dr|  th_Dr  |J_om|  th_om
!     (L × n_Dr rows  =  10 × 25  =  250)
!
!   tcrw_fig3e_summary.txt     — aggregated angles for Fig 3(e), cols:
!     # iD  D_r  Jx_Dr_tot  Jy_Dr_tot  Jx_om_tot  Jy_om_tot  th_Dr_tot  th_om_tot
!     (n_Dr rows  =  25)
!
! Sanity checks (printed to stdout)
! ---------------------------------
!   - θ_JDr_tot at D_r = 10^-4:  should be ≈ −π/2  = −1.5708
!   - θ_JDr_tot at D_r = 10^0:   should be ≈  0
!   - θ_Jω_tot  for all D_r:     should stay near +π/4 = +0.7854
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig3cde.f90 -o tcrw_fig3cde
! Run   :  ./tcrw_fig3cde
! Plot  :  bash run.sh tcrw-plot-fig3c   (quiver J_Dr)
!          bash run.sh tcrw-plot-fig3d   (quiver J_ω)
!          bash run.sh tcrw-plot-fig3e   (θ_JDr vs D_r)
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig3cde
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- panel parameters ----
   integer,  parameter :: L_cur     = 10
   real(dp), parameter :: omega     = 1.0_dp
   integer,  parameter :: n_Dr      = 25
   real(dp), parameter :: log_Dr_min = -4.0_dp
   real(dp), parameter :: log_Dr_max =  0.0_dp
   integer(i8), parameter :: T_floor      = 100000000_i8   ! 10^8
   integer(i8), parameter :: N_burn_floor =  10000000_i8   ! 10^7
   real(dp),    parameter :: K_meas       = 30.0_dp
   real(dp),    parameter :: K_burn       =  5.0_dp
   integer,     parameter :: seed         = 20260423

   ! ---- locals ----
   integer  :: iD, y, u_field, u_angle
   real(dp) :: D_r_cur
   real(dp) :: D_r_values(n_Dr)
   real(dp) :: t0, t1, t_run

   ! Per-site wall currents for the current D_r.  Allocated for L = L_cur.
   integer(i8), allocatable :: Jx_Dr(:), Jy_Dr(:), Jx_om(:), Jy_om(:)
   real(dp)    :: abs_JDr_y, abs_Jom_y, th_JDr_y, th_Jom_y
   integer(i8) :: Jx_Dr_tot, Jy_Dr_tot, Jx_om_tot, Jy_om_tot
   real(dp)    :: th_JDr_tot, th_Jom_tot

   call sgrnd(seed)

   ! ---- build log-spaced D_r grid (matches Fig 3a, b) ----
   call build_log_grid(D_r_values, n_Dr, log_Dr_min, log_Dr_max)

   print '(A)',         '==== TCRW Fig 3(c)/(d)/(e): left-wall J_Dr, J_ω vs D_r (L = 10, ω = 1) ===='
   print '(A,I0,A,ES11.4)', '  L = ', L_cur, '     ω = ', omega
   print '(A,I12,A,I12)',   '  T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   print '(A,F6.1,A,F6.1,A)', &
        '  T_use      = max(T_floor,      ', K_meas, ' * max(L^2, 1/D_r)/D_r)    ' // &
        '  N_burn_use = max(N_burn_floor, ', K_burn, ' * max(L^2, 1/D_r)/D_r)'
   print '(A,I0,A,ES9.2,A,ES9.2)', &
        '  D_r grid  (', n_Dr, ' log-spaced):  ', D_r_values(1), ' ... ', D_r_values(n_Dr)
   print '(A,I0)',         '  seed  = ', seed
   print '(A)',            ''

   allocate(Jx_Dr(0:L_cur-1), Jy_Dr(0:L_cur-1), Jx_om(0:L_cur-1), Jy_om(0:L_cur-1))

   ! ---- open two output files (field + aggregate) ----
   open(newunit=u_field, file='tcrw_fig3cde_summary.txt', status='replace', action='write')
   write(u_field, '(A)') '# TCRW Fig 3(c)/(d): per-site left-wall currents  (L = 10, ω = 1)'
   write(u_field, '(A,F6.1,A,F6.1,A)') &
        '# T_use = max(T_floor, ', K_meas, '*max(L^2,1/D_r)/D_r);  ' // &
        'N_burn_use = max(N_burn_floor, ', K_burn, '*max(L^2,1/D_r)/D_r)'
   write(u_field, '(A,I0,A,I0)') '# L = ', L_cur, '   seed = ', seed
   write(u_field, '(A)') '# columns:  iD  D_r  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om  ' // &
                         '|J_Dr|  th_Dr  |J_om|  th_om'

   open(newunit=u_angle, file='tcrw_fig3e_summary.txt', status='replace', action='write')
   write(u_angle, '(A)') '# TCRW Fig 3(e):  θ_JDr and θ_Jω of the TOTAL left-wall current  (L = 10, ω = 1)'
   write(u_angle, '(A,I0,A,I0)') '# L = ', L_cur, '   seed = ', seed
   write(u_angle, '(A)') '# columns:  iD  D_r  Jx_Dr_tot  Jy_Dr_tot  Jx_om_tot  Jy_om_tot  ' // &
                         'th_Dr_tot  th_om_tot'

   ! ---- sweep D_r ----
   print '(A)', '     iD       D_r          Jy_Dr_tot    Jy_om_tot      th_Dr_tot(rad)   th_om_tot(rad)   cpu[s]'
   print '(A)', '   ----  -----------      ----------   ----------     ---------------  ---------------   ------'

   do iD = 1, n_Dr
      D_r_cur = D_r_values(iD)
      call cpu_time(t0)
      call run_one(omega, D_r_cur, L_cur, T_floor, N_burn_floor, &
                   Jx_Dr, Jy_Dr, Jx_om, Jy_om)
      call cpu_time(t1)
      t_run = t1 - t0

      ! ---- write per-site rows to Fig 3(c)/(d) file ----
      do y = 0, L_cur - 1
         abs_JDr_y = sqrt( real(Jx_Dr(y), dp)**2 + real(Jy_Dr(y), dp)**2 )
         abs_Jom_y = sqrt( real(Jx_om(y), dp)**2 + real(Jy_om(y), dp)**2 )
         if (abs_JDr_y > 0.0_dp) then
            th_JDr_y = atan2( real(Jy_Dr(y), dp), real(Jx_Dr(y), dp) )
         else
            th_JDr_y = 0.0_dp
         end if
         if (abs_Jom_y > 0.0_dp) then
            th_Jom_y = atan2( real(Jy_om(y), dp), real(Jx_om(y), dp) )
         else
            th_Jom_y = 0.0_dp
         end if

         write(u_field, '(I4, 1X, ES13.5, 1X, I3, 4(1X, I14), 1X, ES13.5, 1X, F10.5, 1X, ES13.5, 1X, F10.5)') &
              iD, D_r_cur, y, &
              Jx_Dr(y), Jy_Dr(y), Jx_om(y), Jy_om(y), &
              abs_JDr_y, th_JDr_y, abs_Jom_y, th_Jom_y
      end do
      write(u_field, '(A)') ''     ! blank between D_r blocks (nice for gnuplot index)

      ! ---- aggregate totals over wall for Fig 3(e) ----
      Jx_Dr_tot = sum(Jx_Dr)
      Jy_Dr_tot = sum(Jy_Dr)
      Jx_om_tot = sum(Jx_om)
      Jy_om_tot = sum(Jy_om)

      th_JDr_tot = 0.0_dp
      if (Jx_Dr_tot /= 0_i8 .or. Jy_Dr_tot /= 0_i8) then
         th_JDr_tot = atan2( real(Jy_Dr_tot, dp), real(Jx_Dr_tot, dp) )
      end if
      th_Jom_tot = 0.0_dp
      if (Jx_om_tot /= 0_i8 .or. Jy_om_tot /= 0_i8) then
         th_Jom_tot = atan2( real(Jy_om_tot, dp), real(Jx_om_tot, dp) )
      end if

      write(u_angle, '(I4, 1X, ES13.5, 4(1X, I14), 2(1X, F10.5))') &
           iD, D_r_cur, Jx_Dr_tot, Jy_Dr_tot, Jx_om_tot, Jy_om_tot, &
           th_JDr_tot, th_Jom_tot

      print '(2X, I4, 2X, ES12.4, 2(2X, I12), 2(2X, F14.5), 2X, F7.2)', &
           iD, D_r_cur, Jy_Dr_tot, Jy_om_tot, th_JDr_tot, th_Jom_tot, t_run
   end do

   close(u_field)
   close(u_angle)

   print '(A)', ''
   print '(A,I0,A)', 'Wrote Fig 3(c)/(d) field -> tcrw_fig3cde_summary.txt  (', &
                     L_cur * n_Dr, ' site-rows)'
   print '(A,I0,A)', 'Wrote Fig 3(e) angles    -> tcrw_fig3e_summary.txt    (', &
                     n_Dr, ' rows)'
   print '(A)',      'Plot:  bash run.sh tcrw-plot-fig3c   # quiver J_Dr'
   print '(A)',      '       bash run.sh tcrw-plot-fig3d   # quiver J_ω'
   print '(A)',      '       bash run.sh tcrw-plot-fig3e   # θ_JDr(D_r)'

   deallocate(Jx_Dr, Jy_Dr, Jx_om, Jy_om)

contains

   !------------------------------------------------------------------
   ! Build n log-spaced points between 10^xmin and 10^xmax (inclusive).
   !------------------------------------------------------------------
   subroutine build_log_grid(x, n, xmin, xmax)
      real(dp), intent(out) :: x(:)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: xmin, xmax
      integer  :: k
      real(dp) :: frac
      do k = 1, n
         frac = real(k - 1, dp) / real(n - 1, dp)
         x(k) = 10.0_dp ** ( xmin + (xmax - xmin) * frac )
      end do
   end subroutine build_log_grid

   !------------------------------------------------------------------
   ! Run a single (ω, D_r, L) MC trajectory, accumulating PER-SITE
   ! left-wall currents decomposed into J_Dr and J_ω by the standard
   ! prev_noise rule.
   !------------------------------------------------------------------
   subroutine run_one(omega_in, D_r_in, L_in, T_floor_in, N_burn_floor_in, &
                      Jx_Dr, Jy_Dr, Jx_om, Jy_om)
      real(dp),    intent(in)    :: omega_in, D_r_in
      integer,     intent(in)    :: L_in
      integer(i8), intent(in)    :: T_floor_in, N_burn_floor_in
      integer(i8), intent(inout) :: Jx_Dr(0:), Jy_Dr(0:), Jx_om(0:), Jy_om(0:)

      real(dp)    :: grnd                         ! RNG from mt.f90
      integer     :: x, y, d, x_before, y_before, d_before
      integer(i8) :: it, T_use, N_burn_use
      integer     :: step_type
      real(dp)    :: tau_bulk, tau_wall, tau_relax
      logical     :: prev_noise
      logical, allocatable :: mask(:,:)
      integer, parameter :: DX(0:3) = (/ 0,  1,  0, -1 /)
      integer, parameter :: DY(0:3) = (/ 1,  0, -1,  0 /)

      ! ---- adaptive T and burn-in (same as Fig 3a, 3b) ----
      tau_bulk   = real(L_in, dp) ** 2  / D_r_in
      tau_wall   = 1.0_dp / (D_r_in * D_r_in)
      tau_relax  = max(tau_bulk, tau_wall)
      N_burn_use = max( N_burn_floor_in, int(K_burn * tau_relax, i8) )
      T_use      = max( T_floor_in,      int(K_meas * tau_relax, i8) )

      ! ---- all-true mask = plain OBC on L × L box ----
      allocate(mask(0:L_in-1, 0:L_in-1))
      mask = .true.

      ! ---- random initial site & direction ----
      x = int( real(L_in, dp) * grnd() )
      y = int( real(L_in, dp) * grnd() )
      if (x == L_in) x = L_in - 1
      if (y == L_in) y = L_in - 1
      d = int( 4.0_dp * grnd() )
      if (d == 4) d = 3

      ! ---- burn-in ----
      prev_noise = .false.
      do it = 1_i8, N_burn_use
         call tcrw_step_mask(x, y, d, mask, L_in, L_in, omega_in, D_r_in, step_type)
         prev_noise = (step_type == 0)
      end do

      ! ---- zero accumulators ----
      Jx_Dr = 0_i8 ;  Jy_Dr = 0_i8
      Jx_om = 0_i8 ;  Jy_om = 0_i8

      ! ---- measurement ----
      do it = 1_i8, T_use
         x_before = x
         y_before = y
         d_before = d
         call tcrw_step_mask(x, y, d, mask, L_in, L_in, omega_in, D_r_in, step_type)

         if (step_type == 1 .and. x_before == 0) then
            if (prev_noise) then
               Jx_Dr(y_before) = Jx_Dr(y_before) + int(DX(d_before), i8)
               Jy_Dr(y_before) = Jy_Dr(y_before) + int(DY(d_before), i8)
            else
               Jx_om(y_before) = Jx_om(y_before) + int(DX(d_before), i8)
               Jy_om(y_before) = Jy_om(y_before) + int(DY(d_before), i8)
            end if
         end if

         prev_noise = (step_type == 0)
      end do

      deallocate(mask)
   end subroutine run_one

   !------------------------------------------------------------------
   ! Shared walker-step kernels.
   !------------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig3cde

include 'mt.f90'
