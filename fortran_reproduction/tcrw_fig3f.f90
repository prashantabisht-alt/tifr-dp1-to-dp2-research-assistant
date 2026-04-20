!=====================================================================
! tcrw_fig3f.f90 — Fig 3(f): P_edge/P_bulk vs ω at fixed D_r = 10^-3
!
! Paper: Osat, Meyberg, Metson, Speck,
!        "Topological chiral random walker"
!        arXiv:2602.12020v1 [cond-mat.stat-mech], 12 Feb 2026.
!
! Panel: Fig 3(f) — plot of the per-site edge/bulk ratio vs the noise
!        chirality ω, at fixed rotational noise rate D_r = 10^-3, for
!        system sizes L ∈ {10, 19, 49}.
!
!        Paper caption punchline: the ratio is approximately INDEPENDENT
!        of ω.  Even at ω = 0.5 (no net chirality), edge localization
!        persists at the same value as at ω = 1 (full chirality).  This
!        is the central evidence that edge localization is a TOPOLOGICAL
!        BAND-STRUCTURE effect (non-trivial Chern number of the 4-band
!        Bloch Hamiltonian), NOT a chirality effect.  Chirality sets the
!        direction of the edge current (Fig 3(b), (g)) but not its
!        magnitude or the existence of edge localization.
!
! Convention (matches Python `tcrw_obc.py:fig3f_Pedge_vs_omega`)
! -------------------------------------------------------------
!   - Lattice IS the L × L playground (OBC, no wall ring).
!   - Edge sites: x == 0 .or. x == L-1 .or. y == 0 .or. y == L-1
!     → n_edge = 4L − 4
!   - Bulk sites: 1 ≤ x ≤ L-2 .and. 1 ≤ y ≤ L-2
!     → n_bulk = (L − 2)^2
!   - Per-site ratio: r = (P_edge / n_edge) / (P_bulk / n_bulk)
!   - Uses tcrw_step_obc (same kernel as Fig 3(a)).
!
! Expected physics
! ----------------
!   - r ≈ 700 for ALL ω ∈ [0, 1] at D_r = 10^-3 (matches the Fig 3(a)
!     master curve at that D_r).  Plot should be nearly a horizontal line.
!   - Mirror symmetry: r(ω) = r(1-ω) within MC noise (chain is invariant
!     under CW ↔ CCW relabelling).
!   - All three L curves should overlap — lack of L-dependence.
!   - At ω = 1, the point should agree with Fig 3(a) at D_r = 10^-3 to
!     within MC error (~10%).  Cross-check between the two panels.
!
! Parameters
! ----------
!   D_r           = 10^-3                      ! fixed (paper Fig 3f)
!   L grid        = (10, 19, 49)               ! paper subset
!   ω grid        = 21 linearly-spaced pts  0.0 → 1.0
!   T_floor       = 3·10^8                     ! bumped from 10^8 for ω∈{0,1} cleanup
!   N_burn_floor  = 3·10^7                     ! bumped from 10^7
!   K_meas        = 300                        ! bumped from 100 (2026-04-20)
!   K_burn        = 30                         ! bumped from 10
!   seed          = 20260420                   ! distinct from Fig 3(a)'s 20260419
!
!   Why the bump?  Original run left ω = 0 (L=19) and ω = 1 (L=10) as
!   visible outliers (ratio ~ 450 vs plateau ~ 260).  Root cause: at
!   ω ∈ {0, 1} the noise step is a deterministic 90° rotation — the
!   director cycles {↑→↓←} (or reverse) exactly every 4 steps.  This
!   kills Markov-chain randomness for the angular coordinate, so the
!   lattice-position mixing has to come purely from chiral-move
!   scatter.  That's much slower than generic ω, so we throw 3× more
!   steps at it and 3× longer burn-in.
!
! Why adaptive T and burn-in?
! ---------------------------
!   Two correlation times compete (same reasoning as Fig 3(a)):
!
!     τ_bulk = L^2 / D_r           -- bulk random-walk mixing
!     τ_wall = 1 / D_r^2           -- ω=0 or ω=1 wall-escape time
!
!   τ_relax = max(τ_bulk, τ_wall)
!
!   At D_r = 10^-3:  τ_bulk = {10^5, 3.6×10^5, 2.4×10^6} for L = {10, 19, 49};
!                   τ_wall = 10^6  at ω = 0, 1  (less at intermediate ω,
!                   but we take the conservative max for uniform scaling).
!   Worst τ_relax = 2.4×10^6 at L = 49, so T_use = 2.4×10^8 there.
!   For L = 10, 19 the floor T_use = 10^8 wins.
!
! Cost
! ----
!   ~ 10^10 total RNG-driven OBC steps (plus ~10× less burn-in on worst
!   cells).  At ~15 ns/step: ~3–5 min measurement + 1–2 min burn-in =
!   ~5–10 min total on a modern Mac.  Much lighter than Fig 3(a) because
!   D_r is in a mild regime.
!
! Output
! ------
!   tcrw_fig3f_summary.txt
!     # one row per (L, ω); 63 rows total
!     # columns:  L  ω  ratio  P_edge_norm  P_bulk_norm  n_edge  n_bulk
!     where P_edge_norm = P_edge / n_edge   (per-site occupation, edge)
!           P_bulk_norm = P_bulk / n_bulk   (per-site occupation, bulk)
!           ratio       = P_edge_norm / P_bulk_norm
!
! Sanity checks (printed to stdout)
! ---------------------------------
!   - ω = 0 and ω = 1 points: should give same ratio within MC noise.
!   - All three L curves: should overlap (no finite-size splitting).
!   - At ω = 1 cross-check: compare to tcrw_fig3a_summary.txt row
!     (L = 19, D_r ≈ 10^-3).  Currently different L sets though — Fig 3(a)
!     uses L = 4, 9, 19, 49; Fig 3(f) uses L = 10, 19, 49.  The L = 19
!     row is shared and should match.
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig3f.f90 -o tcrw_fig3f
! Run   :  ./tcrw_fig3f
! Plot  :  bash run.sh tcrw-plot-fig3f
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig3f
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- Fig 3(f) parameters (single source of truth) ----
   real(dp), parameter :: D_r_fixed = 1.0e-3_dp
   integer,  parameter :: L_list(3) = (/ 10, 19, 49 /)
   integer,  parameter :: n_omega   = 21
   real(dp), parameter :: omega_min = 0.0_dp
   real(dp), parameter :: omega_max = 1.0_dp
   integer(i8), parameter :: T_floor      = 300000000_i8   ! 3·10^8  (bumped from 10^8)
   integer(i8), parameter :: N_burn_floor =  30000000_i8   ! 3·10^7  (bumped from 10^7)
   real(dp),    parameter :: K_meas       = 300.0_dp       ! bumped 100 → 300 for ω∈{0,1} cleanup
   real(dp),    parameter :: K_burn       =  30.0_dp       ! bumped  10 →  30
   integer,     parameter :: seed         = 20260420

   ! ---- locals ----
   integer  :: iL, iW, u_sum
   integer  :: L_cur, n_edge, n_bulk
   real(dp) :: omega_cur, P_edge_norm, P_bulk_norm, ratio
   real(dp) :: omega_values(n_omega)
   real(dp) :: t0, t1, t_run

   call sgrnd(seed)

   ! ---- build linearly-spaced ω grid ----
   call build_linear_grid(omega_values, n_omega, omega_min, omega_max)

   print '(A)',         '==== TCRW Fig 3(f): P_edge/P_bulk vs ω  (D_r = 10^-3, OBC) ===='
   print '(A,ES11.4)',  '  D_r (fixed)  = ', D_r_fixed
   print '(A,I12,A,I12)', '  T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   print '(A,F6.1,A,F6.1,A)', &
        '  T_use      = max(T_floor,      ', K_meas, ' * max(L^2, 1/D_r)/D_r)    ' // &
        '  N_burn_use = max(N_burn_floor, ', K_burn, ' * max(L^2, 1/D_r)/D_r)'
   print '(A,I0,A,3(1X,I0))', '  L grid  (', size(L_list), ') :', L_list
   print '(A,I0,A,F4.2,A,F4.2)', &
        '  ω grid  (', n_omega, ' linearly-spaced) :  ', &
        omega_values(1), ' ... ', omega_values(n_omega)
   print '(A,I0)', '  seed  = ', seed
   print '(A)',    ''

   open(newunit=u_sum, file='tcrw_fig3f_summary.txt', status='replace', action='write')
   write(u_sum, '(A)') '# TCRW Fig 3(f)  P_edge/P_bulk vs ω  (D_r = 10^-3, OBC)'
   write(u_sum, '(A,ES11.4)') '# D_r_fixed = ', D_r_fixed
   write(u_sum, '(A,I12,A,I12)') '# T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   write(u_sum, '(A,F6.1,A,F6.1,A)') &
        '# T_use = max(T_floor, ', K_meas, '*max(L^2,1/D_r)/D_r);  ' // &
        'N_burn_use = max(N_burn_floor, ', K_burn, '*max(L^2,1/D_r)/D_r)'
   write(u_sum, '(A,I0)') '# seed    = ', seed
   write(u_sum, '(A)') '# columns:  L   ω   ratio   P_edge_norm   P_bulk_norm   n_edge   n_bulk'

   ! ---- outer loop over L ----
   do iL = 1, size(L_list)
      L_cur  = L_list(iL)
      n_edge = 4 * L_cur - 4
      n_bulk = (L_cur - 2) ** 2

      print '(A,I3,A,I0,A,I0,A)', '  --- L = ', L_cur, &
            '   (n_edge = ', n_edge, ', n_bulk = ', n_bulk, ') ---'
      print '(A)', '       ω        P_edge_norm    P_bulk_norm        ratio       cpu[s]'
      print '(A)', '   --------   -------------  -------------  -------------  -----------'

      ! ---- inner loop over ω ----
      do iW = 1, n_omega
         omega_cur = omega_values(iW)
         call cpu_time(t0)
         call run_one(omega_cur, D_r_fixed, L_cur, T_floor, N_burn_floor, &
                      n_edge, n_bulk, P_edge_norm, P_bulk_norm, ratio)
         call cpu_time(t1)
         t_run = t1 - t0

         print '(2X,F8.4, 2(2X,ES13.5), 2X,ES13.5, 2X,F9.2)', &
              omega_cur, P_edge_norm, P_bulk_norm, ratio, t_run

         write(u_sum, '(I3, 1X, F8.4, 1X, ES13.5, 1X, ES13.5, 1X, ES13.5, 1X, I8, 1X, I10)') &
              L_cur, omega_cur, ratio, P_edge_norm, P_bulk_norm, n_edge, n_bulk
      end do
      print '(A)', ''
   end do

   close(u_sum)
   print '(A,I0,A)', 'Wrote summary -> tcrw_fig3f_summary.txt   (', &
                     size(L_list) * n_omega, ' rows)'
   print '(A)',      'Plot with:    bash run.sh tcrw-plot-fig3f'

contains

   !------------------------------------------------------------------
   ! Build n linearly-spaced points between xmin and xmax (inclusive).
   !------------------------------------------------------------------
   subroutine build_linear_grid(x, n, xmin, xmax)
      real(dp), intent(out) :: x(:)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: xmin, xmax
      integer  :: k
      real(dp) :: frac
      do k = 1, n
         frac = real(k - 1, dp) / real(n - 1, dp)
         x(k) = xmin + (xmax - xmin) * frac
      end do
   end subroutine build_linear_grid

   !------------------------------------------------------------------
   ! Run a single (ω, D_r, L) MC trajectory of length T_use after a
   ! burn-in of N_burn_use discarded steps, where
   !   τ_relax    = max(L^2, 1/D_r) / D_r
   !   N_burn_use = max(N_burn_floor, K_burn · τ_relax)
   !   T_use      = max(T_floor,      K_meas · τ_relax)
   !
   ! Returns:
   !   P_edge_norm = (visits to any edge site)  / (n_edge × T_meas)
   !   P_bulk_norm = (visits to any bulk site)  / (n_bulk × T_meas)
   !   ratio       = P_edge_norm / P_bulk_norm
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
      ! τ_bulk = L^2 / D_r            bulk random-walk mixing
      ! τ_wall = 1 / D_r^2            ω → 0, 1 wall-escape time
      ! τ_relax = max(τ_bulk, τ_wall) — conservative over intermediate ω
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

end program tcrw_fig3f

include 'mt.f90'
