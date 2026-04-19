!=====================================================================
! tcrw_fig3a.f90 — Fig 3(a): P_edge/P_bulk vs D_r for ω = 1, multiple L
!
! Paper: Osat, Meyberg, Metson, Speck,
!        "Topological chiral random walker"
!        arXiv:2602.12020v1 [cond-mat.stat-mech], 12 Feb 2026.
! Panel: Fig 3(a) — log-log plot of the per-site edge/bulk ratio of the
!        steady-state occupation probability vs. the rotational noise
!        rate D_r, at maximal chirality ω = 1, for system sizes
!        L ∈ {4, 9, 19, 49}.  The paper goes up to L = 99 (and higher)
!        using exact diagonalization; MC at L = 99 with the small-D_r
!        end of the range would need T ≳ 10^9 for a clean ratio, so we
!        stop at L = 49 here.
!        Caption: ratio collapses onto a master curve, independent of L
!                 once normalized per site → "edge localization is a
!                 bulk-band-topology effect, not a finite-size artefact".
!
! Convention (matches Python `tcrw_obc.py`)
! -----------------------------------------
!   - Lattice IS the L × L playground.  Walker occupies any site in
!     [0, L-1] × [0, L-1].  No surrounding wall ring.
!   - Edge sites: x == 0 .or. x == L-1 .or. y == 0 .or. y == L-1
!     → n_edge = 4L − 4
!   - Bulk sites: 1 ≤ x ≤ L-2 .and. 1 ≤ y ≤ L-2
!     → n_bulk = (L − 2)²
!   - Per-site ratio:   r = (P_edge / n_edge) / (P_bulk / n_bulk)
!   - We use the same OBC step kernel as Fig 1: `tcrw_step_obc`
!     (NOT the wall-ringed mask kernel used in Fig 2 defects).
!
! Expected physics
! ----------------
!   - At small D_r the ω = 1 walker spirals into a tight CW loop near
!     each wall and the ratio r grows as roughly 1 / D_r.  Master curve
!     is a power-law on the log-log plot.
!   - At D_r ≳ 0.1 the noise rate exceeds the chiral-loop rate, the
!     edge trap is destroyed, r → 1 (uniform distribution).
!   - All four L curves should overlap ABOVE D_r ~ 10^-3 because the
!     edge-localization length stays small compared to L.  At the
!     smallest D_r the L = 4 curve will saturate (the whole system IS
!     the edge layer once the edge-localization length exceeds L/2).
!
! Parameters
! ----------
!   ω             = 1.0_dp                         ! fixed (paper Fig 3a)
!   L grid        = (4, 9, 19, 49)                 ! paper subset
!   D_r grid      = 25 log-spaced points, 10^-4 → 10^0 (exact paper grid)
!   T_floor       = 10^8                           ! minimum measurement length
!   N_burn_floor  = 10^7                           ! minimum discarded burn-in
!   K_meas        = 100                            ! T_use    = max(T_floor,    K_meas · L^2/D_r)
!   K_burn        = 10                             ! N_burn   = max(N_burn_floor, K_burn · L^2/D_r)
!   seed          = 20260419
!
! Why adaptive T and burn-in?
! ---------------------------
!   Two correlation times govern the chain at ω = 1:
!
!     τ_bulk = L^2 / D_r          -- bulk random-walk mixing
!                                    (need L^2 noise events to diffuse
!                                     across the box; one noise per 1/D_r
!                                     MC steps)
!
!     τ_wall = 1 / D_r^2          -- wall-trap escape time
!                                    (from a wall-bound state, two
!                                     consecutive noise events are
!                                     needed to flip d into a bulk-
!                                     pointing direction; each noise
!                                     attempt waits ~1/D_r steps)
!
!   The slower one sets the autocorrelation of the edge/bulk indicator:
!
!     τ_relax = max(τ_bulk, τ_wall)
!
!   At small D_r (D_r ≲ 1/L^2) the wall trap dominates; at large D_r
!   (D_r ≳ 1/L^2) bulk diffusion dominates.  For our worst corner
!   (L = 49, D_r = 10^-4), 1/D_r^2 = 10^8 ≫ L^2/D_r = 2.4×10^7, so
!   τ_relax ≈ 10^8 — a hundred times larger than the naive L^2/D_r
!   estimate.  This is why a first pass that scaled T as L^2/D_r still
!   produced ratty curves at small D_r: only ~1 wall-escape event was
!   sampled at L = 4, 9, where the floor capped T at 10^8.
!
!   Scaling both burn-in and T with τ_relax = max(L^2, 1/D_r) / D_r:
!     - equilibrates each (L, D_r) for 10 τ before measurement
!     - averages over 100 τ (giving ~10% statistical error per point)
!
! Cost
! ----
!   Worst cells are now the (any L, D_r ≤ 3×10^-4) cells, each costing
!   ~10^10 steps.  At ~15 ns/step that's ~150 s per cell × 4 L = ~10 min
!   for the D_r = 10^-4 column alone.  Summing the geometric series over
!   the 6 smallest D_r values: ~20 min total.  All other cells run at the
!   floor cost and finish in seconds.
!
! Output
! ------
!   tcrw_fig3a_summary.txt
!     # one row per (L, D_r); 100 rows total
!     # columns:  L  D_r  ratio  P_edge_norm  P_bulk_norm  n_edge  n_bulk
!     where P_edge_norm = P_edge / n_edge   (per-site occupation, edge)
!           P_bulk_norm = P_bulk / n_bulk   (per-site occupation, bulk)
!           ratio       = P_edge_norm / P_bulk_norm
!
! Sanity checks (will be visible in stdout)
! -----------------------------------------
!   - At D_r = 1 every L curve should give ratio ≈ 1.0  (well-mixed).
!   - For L = 4 at D_r = 10^-4 ratio should saturate (small-L cap).
!   - For L = 49 at D_r = 10^-4 ratio should be the largest of the four.
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig3a.f90 -o tcrw_fig3a
! Run   :  ./tcrw_fig3a
! Plot  :  bash run.sh tcrw-plot-fig3a
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig3a
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- Fig 3(a) parameters (single source of truth) ----
   real(dp), parameter :: omega   = 1.0_dp
   integer,  parameter :: L_list(4) = (/ 4, 9, 19, 49 /)
   integer,  parameter :: n_Dr   = 25
   real(dp), parameter :: log_Dr_min = -4.0_dp
   real(dp), parameter :: log_Dr_max =  0.0_dp
   integer(i8), parameter :: T_floor      = 100000000_i8   ! 10^8  (floor for T_steps)
   integer(i8), parameter :: N_burn_floor =  10000000_i8   ! 10^7  (floor for burn-in)
   real(dp),    parameter :: K_meas       = 100.0_dp       ! measurement = K_meas × τ_relax
   real(dp),    parameter :: K_burn       =  10.0_dp       ! burn-in     = K_burn × τ_relax
   integer,     parameter :: seed         = 20260419
   ! Both T and burn-in are set per-(L, D_r) inside run_one as
   !     τ_relax    = L^2 / D_r                              (mixing time)
   !     N_burn_use = max( N_burn_floor, K_burn · τ_relax )  (10 τ of equilibration)
   !     T_use      = max( T_floor,      K_meas · τ_relax )  (100 τ of measurement)
   ! At the (L=49, D_r=10^-4) corner τ_relax ≈ 2.4×10^7, so T_use ≈ 2.4×10^9
   ! (≈24× the floor) and N_burn_use ≈ 2.4×10^8 (≈24× the floor).  For
   ! (L, D_r) cells where L^2/D_r is small, both floors win and the cell
   ! costs the same as before.  Total runtime: ~6–8 min (vs ~2–3 min with
   ! fixed T = 10^8).

   ! ---- locals ----
   integer  :: iL, iD, u_sum
   integer  :: L_cur, n_edge, n_bulk
   real(dp) :: D_r_cur, P_edge_norm, P_bulk_norm, ratio
   real(dp) :: D_r_values(n_Dr)
   real(dp) :: t0, t1, t_run

   call sgrnd(seed)

   ! ---- build log-spaced D_r grid ----
   call build_log_grid(D_r_values, n_Dr, log_Dr_min, log_Dr_max)

   print '(A)',         '==== TCRW Fig 3(a): P_edge/P_bulk vs D_r (ω = 1) ===='
   print '(A,I12,A,I12)','  T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   print '(A,F6.1,A,F6.1)', &
        '  T_use      = max(T_floor,      ', K_meas, ' * L^2/D_r)    ' // &
        '  N_burn_use = max(N_burn_floor, ', K_burn, ' * L^2/D_r)'
   print '(A,I0,A,4(1X,I0))', '  L grid (', size(L_list), ') :', L_list
   print '(A,I0,A,ES9.2,A,ES9.2)', &
        '  D_r grid (', n_Dr, ' log-spaced):  ', D_r_values(1), '  ...  ', D_r_values(n_Dr)
   print '(A,I0)', '  seed  = ', seed
   print '(A)',    ''

   open(newunit=u_sum, file='tcrw_fig3a_summary.txt', status='replace', action='write')
   write(u_sum, '(A)') '# TCRW Fig 3(a)  P_edge/P_bulk vs D_r  (ω = 1, OBC)'
   write(u_sum, '(A,I12,A,I12)') '# T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   write(u_sum, '(A,F6.1,A,F6.1,A)') &
        '# T_use = max(T_floor, ', K_meas, '*L^2/D_r);  ' // &
        '  N_burn_use = max(N_burn_floor, ', K_burn, '*L^2/D_r)'
   write(u_sum, '(A,I0)')        '# seed    = ', seed
   write(u_sum, '(A)') '# columns:  L   D_r   ratio   P_edge_norm   P_bulk_norm   n_edge   n_bulk'

   ! ---- outer loop over L ----
   do iL = 1, size(L_list)
      L_cur  = L_list(iL)
      n_edge = 4 * L_cur - 4
      n_bulk = (L_cur - 2) ** 2

      print '(A,I3,A,I0,A,I0,A)', '  --- L = ', L_cur, &
            '   (n_edge = ', n_edge, ', n_bulk = ', n_bulk, ') ---'
      print '(A)', '       D_r        P_edge_norm    P_bulk_norm        ratio       cpu[s]'
      print '(A)', '   -----------  -------------  -------------  -------------  -----------'

      ! ---- inner loop over D_r ----
      do iD = 1, n_Dr
         D_r_cur = D_r_values(iD)
         call cpu_time(t0)
         call run_one(omega, D_r_cur, L_cur, T_floor, N_burn_floor, &
                      n_edge, n_bulk, P_edge_norm, P_bulk_norm, ratio)
         call cpu_time(t1)
         t_run = t1 - t0

         print '(2X,ES12.4, 2(2X,ES13.5), 2X,ES13.5, 2X,F9.2)', &
              D_r_cur, P_edge_norm, P_bulk_norm, ratio, t_run

         write(u_sum, '(I3, 1X, ES13.5, 1X, ES13.5, 1X, ES13.5, 1X, ES13.5, 1X, I8, 1X, I10)') &
              L_cur, D_r_cur, ratio, P_edge_norm, P_bulk_norm, n_edge, n_bulk
      end do
      print '(A)', ''
   end do

   close(u_sum)
   print '(A,I0,A)', 'Wrote summary -> tcrw_fig3a_summary.txt   (', &
                     size(L_list) * n_Dr, ' rows)'
   print '(A)',      'Plot with:    bash run.sh tcrw-plot-fig3a'

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
   ! Run a single (ω, D_r, L) MC trajectory of length T_steps after a
   ! burn-in of N_burn discarded steps.  Walker uses the OBC kernel
   ! (no wall ring; lattice IS the L × L playground).
   !
   ! Returns:
   !   P_edge_norm = (visits to any edge site)  / (n_edge × T_meas)
   !   P_bulk_norm = (visits to any bulk site)  / (n_bulk × T_meas)
   !   ratio       = P_edge_norm / P_bulk_norm
   ! These are *per-site* occupations, so the ratio is L-independent
   ! once the system is large enough (matches Python `compute_edge_bulk_ratio`
   ! with `per_site=True`).
   !------------------------------------------------------------------
   subroutine run_one(omega_in, D_r_in, L_in, T_floor_in, N_burn_floor_in, &
                      n_edge_in, n_bulk_in, P_edge_norm, P_bulk_norm, ratio)
      real(dp),    intent(in)  :: omega_in, D_r_in
      integer,     intent(in)  :: L_in, n_edge_in, n_bulk_in
      integer(i8), intent(in)  :: T_floor_in, N_burn_floor_in
      real(dp),    intent(out) :: P_edge_norm, P_bulk_norm, ratio

      real(dp)    :: grnd                           ! external RNG (mt.f90)
      integer     :: x, y, d
      integer(i8) :: it, edge_cnt, bulk_cnt, T_meas, T_use, N_burn_use
      real(dp)    :: tau_bulk, tau_wall, tau_relax
      logical     :: is_edge

      ! ---- adaptive burn-in and measurement length ----
      ! Two competing correlation times govern the chain at ω = 1:
      !
      !   τ_bulk = L^2 / D_r          -- bulk random-walk mixing
      !                                  (noise events random-walk the walker
      !                                   across the L × L box in ~L^2 noise
      !                                   events, i.e. L^2/D_r MC steps)
      !
      !   τ_wall = 1 / D_r^2          -- ω=1 wall-escape time
      !                                  (from d=3 on left wall:
      !                                     wait 1/D_r steps for noise → d=2,
      !                                     then P[next step is noise → d=1] = D_r,
      !                                     else walker slides along wall and
      !                                     restarts from d=3.  Mean time to
      !                                     escape to bulk = 1/D_r^2.)
      !
      ! The slower of the two sets the autocorrelation time of the
      ! edge/bulk indicator:
      !
      !   τ_relax = max(τ_bulk, τ_wall) = max(L^2, 1/D_r) / D_r
      !
      ! At D_r = 10^-4 with L ≤ 49 we have 1/D_r = 10^4 > L^2 for all L,
      ! so τ_wall wins and τ_relax ≈ 10^8 — 100x longer than the naive
      ! L^2/D_r estimate.  This is why the first adaptive-T pass still
      ! produced ratty curves at small D_r (only ~1 wall-escape event
      ! sampled at L = 4, 9).
      tau_bulk   = real(L_in, dp) ** 2  / D_r_in
      tau_wall   = 1.0_dp / (D_r_in * D_r_in)
      tau_relax  = max(tau_bulk, tau_wall)
      N_burn_use = max( N_burn_floor_in, int(K_burn * tau_relax, i8) )
      T_use      = max( T_floor_in,      int(K_meas * tau_relax, i8) )

      ! ---- random initial site & direction ----
      x = int( real(L_in, dp) * grnd() )
      y = int( real(L_in, dp) * grnd() )
      if (x == L_in) x = L_in - 1
      if (y == L_in) y = L_in - 1
      d = int( 4.0_dp * grnd() )
      if (d == 4) d = 3

      ! ---- burn-in (discard) ----
      do it = 1_i8, N_burn_use
         call tcrw_step_obc(x, y, d, L_in, omega_in, D_r_in)
      end do

      ! ---- measurement ----
      edge_cnt = 0_i8
      bulk_cnt = 0_i8
      do it = 1_i8, T_use
         call tcrw_step_obc(x, y, d, L_in, omega_in, D_r_in)
         is_edge = (x == 0 .or. x == L_in - 1 .or. y == 0 .or. y == L_in - 1)
         if (is_edge) then
            edge_cnt = edge_cnt + 1_i8
         else
            bulk_cnt = bulk_cnt + 1_i8
         end if
      end do
      T_meas = edge_cnt + bulk_cnt              ! == T_use by construction

      P_edge_norm = real(edge_cnt, dp) / ( real(n_edge_in, dp) * real(T_meas, dp) )
      P_bulk_norm = real(bulk_cnt, dp) / ( real(n_bulk_in, dp) * real(T_meas, dp) )

      if (P_bulk_norm > 0.0_dp) then
         ratio = P_edge_norm / P_bulk_norm
      else
         ratio = huge(1.0_dp)                   ! walker never visited bulk
      end if
   end subroutine run_one


   !------------------------------------------------------------------
   ! Shared walker-step kernels (tcrw_step_unbounded / tcrw_step_obc /
   ! tcrw_step_mask).
   !------------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig3a

include 'mt.f90'
