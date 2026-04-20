!=====================================================================
! tcrw_fig3g.f90 — Fig 3(g): |J_Dr|_wall / |J_ω|_wall  vs  ω
!                            (D_r = 10^-3, multiple L)
!
! Paper: Osat, Meyberg, Metson, Speck,
!        "Topological chiral random walker"
!        arXiv:2602.12020v1 [cond-mat.stat-mech], 12 Feb 2026.
!
! Panel: Fig 3(g) — linear–log plot of the ratio of the TOTAL
!        magnitudes of the two left-wall current components, at
!        fixed D_r = 10^-3, as a function of the chirality ω ∈ [0, 1],
!        for system sizes L ∈ {10, 19, 49}.
!
!        Paper claim (page 3):
!          "|J_Dr| decreases several orders of magnitude when going
!           from a fully chiral ω = 0(1) to an achiral ω = 0.5 case.
!           However, |J_ω| is almost constant. A symmetric behavior
!           is observed where the edge current decreases as we span
!           the two extreme CRW chiralities ω = 0 and ω = 1.  At
!           ω = 0.5 the walker behaves almost like a normal random
!           walker..."
!        So the RATIO |J_Dr|/|J_ω| shows a deep symmetric DIP at
!        ω = 0.5.  It's the ω-sweep twin of Fig 3(b)'s D_r-sweep.
!
! Current decomposition (identical to Fig 2, Fig 3(b))
! ----------------------------------------------------
!   step_type ∈ {0, 1, 2}  from tcrw_step_mask:
!       0 = noise step         (rotate only; no translation)
!       1 = chiral success     (translate by (DX(d_before), DY(d_before))
!                               then rotate)
!       2 = chiral blocked     (target outside L × L box; skipped)
!   prev_noise ← (step_type == 0)    carried across all steps.
!
!   At step_type == 1 the translation is attributed to:
!       J_Dr  if prev_noise  (noise-then-chiral: post-scatter step)
!       J_ω   otherwise      (chiral-then-chiral: in-orbit step)
!
!   We sum signed (Δx, Δy) over LEFT-WALL sites only (x_before == 0)
!   into four int64 scalar accumulators.
!
! Expected sign pattern  (derived from orbit structure, matches paper)
! --------------------------------------------------------------------
!   Rotation encoding:  d = 0↑, 1→, 2↓, 3←.
!
!   At ω = 1 (chiral = CW director, spatial orbit CCW):
!       walker stuck at (0, y) facing ←; noise rotates CCW to ↓; next
!       chiral step translates DOWN.  Each noise event drifts walker
!       DOWN by one → Jy_w_Dr < 0   (θ_JDr = -π/2,  paper Fig 3(e))
!       orbit-repeat at x=0 departs UP → Jy_w_om > 0
!                                      (θ_Jω = +π/4,  paper Fig 3(j))
!
!   At ω = 0 (chiral = CCW director, spatial orbit CW):
!       stuck at (0, y) facing ←; noise rotates CW (noise uses 1-ω rule
!       at ω=0, which is CW) to ↑; next chiral translates UP.  Each
!       noise event drifts walker UP by one → Jy_w_Dr > 0  (+π/2).
!       orbit-repeat at x=0 departs DOWN → Jy_w_om < 0  (-π/4).
!
!   At ω = 0.5 (achiral): walker is an ordinary persistent RW.  Signed
!       sums of (Δx, Δy) over wall visits have ZERO MEAN; |J_Dr| and
!       |J_ω| both scale as √(#visits), MC noise dominates the ratio.
!
!   Therefore:
!       sign(Jy_w_om) expected:    +1  for ω > 0.5,   -1  for ω < 0.5
!       sign(Jy_w_Dr) expected:    -1  for ω > 0.5,   +1  for ω < 0.5
!       both near zero and sign-noisy in a window around ω = 0.5.
!
! Convention / kernel
! -------------------
!   - Lattice IS the L × L box, OBC, no wall ring.
!   - tcrw_step_mask with an all-TRUE mask (exposes step_type).
!   - Adaptive τ_relax = max(L^2/D_r, 1/D_r^2), same as Fig 3(a,b,f).
!
! Parameters
! ----------
!   D_r           = 10^-3                      ! fixed (paper Fig 3g)
!   L grid        = (10, 19, 49)               ! paper subset (matches Fig 3f)
!   ω grid        = 21 linearly-spaced pts 0.0 → 1.0
!   T_floor       = 10^8                       ! minimum measurement length
!   N_burn_floor  = 10^7                       ! minimum discarded burn-in
!   K_meas        = 100                        ! T_use  = max(T_floor, K_meas · τ_relax)
!   K_burn        = 10                         ! N_burn = max(N_burn_floor, K_burn · τ_relax)
!   seed          = 20260422                   ! distinct from Fig 3(a/b/f)
!
! Cost
! ----
!   At D_r = 10^-3, τ_wall = 1/D_r^2 = 10^6, τ_bulk up to L^2/D_r.
!   For L = 49: τ_bulk = 2.4×10^6, τ_relax = 2.4×10^6, so T_use = 2.4×10^8.
!   For L = 10, 19: T_use = 10^8 (floor wins).
!   Worst burn-in: 2.4×10^7 (L = 49).
!   Total: 63 cells × ~10^8 avg ≈ 6×10^9 RNG-driven steps ≈ 10–12 min
!   on a modern Mac.
!
! Output
! ------
!   tcrw_fig3g_summary.txt
!     # columns:  L   ω   ratio   |J_Dr|_wall   |J_ω|_wall   Jx_w_Dr   Jy_w_Dr   Jx_w_om   Jy_w_om
!     63 rows total (3 L × 21 ω).
!
! Sanity checks (printed to stdout)
! ---------------------------------
!   - sign_Jy_om:   +1 for ω > 0.5,   -1 for ω < 0.5   (may flicker near 0.5)
!   - sign_Jy_Dr:   -1 for ω > 0.5,   +1 for ω < 0.5   (may flicker near 0.5)
!   - ratio(ω) ≈ ratio(1-ω) within MC noise  (CW↔CCW relabelling symmetry)
!   - Cross-check:  ratio at ω = 1 should equal Fig 3(b) ratio at
!     D_r = 10^-3 (same L) within MC error — both panels measure the
!     same quantity at their common corner.
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig3g.f90 -o tcrw_fig3g
! Run   :  ./tcrw_fig3g
! Plot  :  bash run.sh tcrw-plot-fig3g
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig3g
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- Fig 3(g) parameters ----
   real(dp), parameter :: D_r_fixed = 1.0e-3_dp
   integer,  parameter :: L_list(3) = (/ 9, 19, 49 /)   ! paper Fig 3 legend (authors' convention); ⇒ 10×10, 20×20, 50×50 grids
   integer,  parameter :: n_omega   = 21
   real(dp), parameter :: omega_min = 0.0_dp
   real(dp), parameter :: omega_max = 1.0_dp
   integer(i8), parameter :: T_floor      = 100000000_i8   ! 10^8
   integer(i8), parameter :: N_burn_floor =  10000000_i8   ! 10^7
   real(dp),    parameter :: K_meas       = 100.0_dp
   real(dp),    parameter :: K_burn       =  10.0_dp
   integer,     parameter :: seed         = 20260422

   ! ---- locals ----
   integer     :: iL, iW, u_sum
   ! Authors' lattice convention: label L means (L+1)×(L+1) sites (indices 0..L).
   ! L_paper = paper legend label (written to output); L_cur = actual sites per side.
   integer     :: L_paper, L_cur
   real(dp)    :: omega_cur, abs_JDr, abs_Jom, ratio
   real(dp)    :: omega_values(n_omega)
   real(dp)    :: t0, t1, t_run
   integer(i8) :: Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om

   call sgrnd(seed)

   ! ---- build linearly-spaced ω grid ----
   call build_linear_grid(omega_values, n_omega, omega_min, omega_max)

   print '(A)',         '==== TCRW Fig 3(g): |J_Dr|_wall / |J_ω|_wall vs ω  (D_r = 10^-3, OBC) ===='
   print '(A,ES11.4)',  '  D_r (fixed)  = ', D_r_fixed
   print '(A,I12,A,I12)', '  T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   print '(A,F6.1,A,F6.1,A)', &
        '  T_use      = max(T_floor,      ', K_meas, ' * max(L^2, 1/D_r)/D_r)    ' // &
        '  N_burn_use = max(N_burn_floor, ', K_burn, ' * max(L^2, 1/D_r)/D_r)'
   print '(A,I0,A,3(1X,I0))', '  L grid  (', size(L_list), ') :', L_list
   print '(A,I0,A,F4.2,A,F4.2)', &
        '  ω grid  (', n_omega, ' linearly-spaced) :  ', &
        omega_values(1), ' ... ', omega_values(n_omega)
   print '(A,I0)',              '  seed  = ', seed
   print '(A)',                 ''

   open(newunit=u_sum, file='tcrw_fig3g_summary.txt', status='replace', action='write')
   write(u_sum, '(A)') '# TCRW Fig 3(g)  |J_Dr|_wall / |J_ω|_wall vs ω  (D_r = 10^-3, OBC)'
   write(u_sum, '(A,ES11.4)') '# D_r_fixed = ', D_r_fixed
   write(u_sum, '(A,I12,A,I12)') '# T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   write(u_sum, '(A,F6.1,A,F6.1,A)') &
        '# T_use = max(T_floor, ', K_meas, '*max(L^2,1/D_r)/D_r);  ' // &
        'N_burn_use = max(N_burn_floor, ', K_burn, '*max(L^2,1/D_r)/D_r)'
   write(u_sum, '(A,I0)') '# seed    = ', seed
   write(u_sum, '(A)') '# columns: L   ω   ratio   |J_Dr|_wall   |J_ω|_wall   ' // &
                        'Jx_w_Dr   Jy_w_Dr   Jx_w_om   Jy_w_om'

   ! ---- outer loop over L ----
   do iL = 1, size(L_list)
      L_paper = L_list(iL)
      L_cur   = L_paper + 1                ! authors' convention: L=N ⇒ (N+1)×(N+1) sites

      print '(A,I3,A,I0,A)', '  --- L = ', L_paper, ' (sites = ', L_cur, ') ---'
      print '(A)', '        ω         |J_Dr|_wall     |J_ω|_wall         ratio       sJyom  sJyDr   cpu[s]'
      print '(A)', '   ----------   -------------  -------------  -------------   -----  -----   -------'

      ! ---- inner loop over ω ----
      do iW = 1, n_omega
         omega_cur = omega_values(iW)
         call cpu_time(t0)
         call run_one(omega_cur, D_r_fixed, L_cur, T_floor, N_burn_floor, &
                      Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om,                  &
                      abs_JDr, abs_Jom, ratio)
         call cpu_time(t1)
         t_run = t1 - t0

         print '(2X,F10.4, 3(2X,ES13.5), 2X,I4,3X,I4,4X,F7.2)', &
              omega_cur, abs_JDr, abs_Jom, ratio, &
              sign_of(Jy_w_om), sign_of(Jy_w_Dr), t_run

         write(u_sum, '(I3, 1X, F8.4, 1X, ES13.5, 1X, ES13.5, 1X, ES13.5, 4(1X, I14))') &
              L_paper, omega_cur, ratio, abs_JDr, abs_Jom, &
              Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om
      end do
      print '(A)', ''
   end do

   close(u_sum)
   print '(A,I0,A)', 'Wrote summary -> tcrw_fig3g_summary.txt   (', &
                     size(L_list) * n_omega, ' rows)'
   print '(A)',      'Plot with:    bash run.sh tcrw-plot-fig3g'

contains

   !------------------------------------------------------------------
   ! sign of an int64  (-1, 0, +1)  — used to verify the paper's
   ! predicted sign pattern on the left wall (see header).
   !------------------------------------------------------------------
   integer function sign_of(n) result(s)
      integer(i8), intent(in) :: n
      if      (n > 0_i8) then ; s =  1
      else if (n < 0_i8) then ; s = -1
      else                    ; s =  0
      end if
   end function sign_of

   !------------------------------------------------------------------
   ! n linearly-spaced points in [xmin, xmax] (inclusive).
   !------------------------------------------------------------------
   subroutine build_linear_grid(x, n, xmin, xmax)
      real(dp), intent(out) :: x(:)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: xmin, xmax
      integer :: k
      do k = 1, n
         x(k) = xmin + (xmax - xmin) * real(k - 1, dp) / real(n - 1, dp)
      end do
   end subroutine build_linear_grid

   !------------------------------------------------------------------
   ! Run a single (ω, D_r, L) MC trajectory of length T_use steps
   ! after N_burn_use burn-in steps.  Accumulate the two left-wall
   ! currents J_Dr and J_ω decomposed by the "prev_noise" rule.
   !
   !   τ_bulk    = L^2 / D_r
   !   τ_wall    = 1 / D_r^2               (ω-independent upper bound)
   !   τ_relax   = max(τ_bulk, τ_wall)
   !   N_burn_use = max(N_burn_floor_in, K_burn · τ_relax)
   !   T_use      = max(T_floor_in,      K_meas · τ_relax)
   !
   ! Note: τ_wall = 1/D_r^2 is the ω = 0 or ω = 1 edge-trap time.  At
   ! intermediate ω the actual wall-escape time is shorter, but we use
   ! the worst-case bound for uniform scaling across the ω grid.  This
   ! costs us a constant factor ~(1 + D_r²·L²) in runtime at ω = 0.5
   ! but avoids under-sampling the chiral extremes.
   !------------------------------------------------------------------
   subroutine run_one(omega_in, D_r_in, L_in, T_floor_in, N_burn_floor_in, &
                      Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om,                   &
                      abs_JDr, abs_Jom, ratio)
      real(dp),    intent(in)  :: omega_in, D_r_in
      integer,     intent(in)  :: L_in
      integer(i8), intent(in)  :: T_floor_in, N_burn_floor_in
      integer(i8), intent(out) :: Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om
      real(dp),    intent(out) :: abs_JDr, abs_Jom, ratio

      real(dp)    :: grnd                     ! external RNG (mt.f90)
      integer     :: x, y, d, x_before, d_before
      integer(i8) :: it, T_use, N_burn_use
      integer     :: step_type
      real(dp)    :: tau_bulk, tau_wall, tau_relax
      logical     :: prev_noise
      logical, allocatable :: mask(:,:)
      integer, parameter :: DX(0:3) = (/ 0,  1,  0, -1 /)
      integer, parameter :: DY(0:3) = (/ 1,  0, -1,  0 /)

      ! ---- adaptive burn-in and measurement length ----
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
         if (step_type /= 2) prev_noise = (step_type == 0)   ! authors' rule: blocked chiral leaves flag unchanged
      end do

      ! ---- measurement ----
      Jx_w_Dr = 0_i8
      Jy_w_Dr = 0_i8
      Jx_w_om = 0_i8
      Jy_w_om = 0_i8

      do it = 1_i8, T_use
         x_before = x
         d_before = d
         call tcrw_step_mask(x, y, d, mask, L_in, L_in, omega_in, D_r_in, step_type)

         if (step_type == 1 .and. x_before == 0) then
            if (prev_noise) then
               Jx_w_Dr = Jx_w_Dr + int(DX(d_before), i8)
               Jy_w_Dr = Jy_w_Dr + int(DY(d_before), i8)
            else
               Jx_w_om = Jx_w_om + int(DX(d_before), i8)
               Jy_w_om = Jy_w_om + int(DY(d_before), i8)
            end if
         end if

         if (step_type /= 2) prev_noise = (step_type == 0)   ! authors' rule: blocked chiral leaves flag unchanged
      end do

      deallocate(mask)

      abs_JDr = sqrt( real(Jx_w_Dr, dp)**2 + real(Jy_w_Dr, dp)**2 )
      abs_Jom = sqrt( real(Jx_w_om, dp)**2 + real(Jy_w_om, dp)**2 )

      if (abs_Jom > 0.0_dp) then
         ratio = abs_JDr / abs_Jom
      else
         ratio = huge(1.0_dp)
      end if
   end subroutine run_one

   !------------------------------------------------------------------
   ! Shared walker-step kernels.
   !------------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig3g

include 'mt.f90'
