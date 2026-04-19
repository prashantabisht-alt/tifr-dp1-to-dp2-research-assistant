!=====================================================================
! tcrw_fig3b.f90 — Fig 3(b): |J_Dr|_wall / |J_ω|_wall  vs  D_r
!                            (ω = 1, multiple L)
!
! Paper: Osat, Meyberg, Metson, Speck,
!        "Topological chiral random walker"
!        arXiv:2602.12020v1 [cond-mat.stat-mech], 12 Feb 2026.
!
! Panel: Fig 3(b) — log-log plot of the ratio of the TOTAL magnitudes
!        of the two current components on the LEFT WALL of the clean
!        OBC box, at maximal chirality ω = 1, for system sizes
!        L ∈ {4, 9, 19, 49}.
!
!        Paper claim (page 3):
!          "The ratio is constant for small values of D_r and starts
!           to diverge for higher values of D_r.  This seems to
!           contradict the fact that higher D_r values increase the
!           scattering of the walker from the wall.  However, at higher
!           D_r many of the steps are counted as steps after rotational
!           noise, hence, the ratio |J_Dr|/|J_ω| diverges despite the
!           lack of a strong edge current."
!        In short: the D_r → 1 divergence is a BOOKKEEPING effect, not
!        a physical edge-current effect.  At D_r near 1 almost every
!        chiral translation is preceded by a noise event and is thus
!        attributed to J_Dr by construction.
!
! Current decomposition (same rule as Fig 2)
! ------------------------------------------
!   At each MC step, tcrw_step_mask reports step_type ∈ {0, 1, 2}:
!       0 = noise step         (rotate only; no translation)
!       1 = chiral success     (translate by (DX(d_before), DY(d_before))
!                               then rotate)
!       2 = chiral blocked     (target was outside the L × L box;
!                               step skipped)
!
!   We carry a flag prev_noise across steps:
!       prev_noise ← (step_type == 0)
!
!   For a step with step_type == 1 (a successful translation) we bin
!   the translation into one of two current fields at the DEPARTING
!   site (x_before, y_before):
!       prev_noise  → contributes to  J_Dr  (noise-triggered current)
!       otherwise   → contributes to  J_ω   (chiral-continuation current)
!
!   This matches the paper's "it contributes to J_Dr if the previous
!   step was noise, otherwise J_ω" convention and the already-tested
!   tcrw_fig2_clean pipeline.
!
! Observable (Fig 3(b))
! ---------------------
!   Sum J_Dr and J_ω over the LEFT-WALL sites only.  "Left wall" =
!   { (x = 0, y = 0, …, L-1) }.  Because we only need the SUM (not
!   the 2D field), we accumulate four int64 scalar counters:
!       Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om
!   incremented only when x_before == 0 and step_type == 1.
!
!   The reported ratio is
!       |J_Dr|_wall   = sqrt( Jx_w_Dr^2 + Jy_w_Dr^2 )
!       |J_ω|_wall    = sqrt( Jx_w_om^2 + Jy_w_om^2 )
!       ratio         = |J_Dr|_wall / |J_ω|_wall
!   (units: number of translations, so the ratio is dimensionless.)
!
! Expected physics
! ----------------
!   - Left wall carries a CCW chiral edge current at ω = 1; on the
!     left wall CCW means walker slides DOWNWARD, so the dominant
!     contribution to J_ω on the left wall is in the -y direction
!     (J_y_om < 0, |J_y_om| large).
!   - J_Dr on the left wall is the post-noise current.  At small D_r
!     it's weak (noise is rare); at large D_r every other step is a
!     noise followed by a scattering step.
!   - Ratio at small D_r:  order 1, D_r-independent plateau.
!   - Ratio at D_r → 1:   diverges as roughly  D_r / (1 - D_r)
!     (trivial counting limit — probability of "chiral-then-chiral" is
!      (1-D_r)^2, of "noise-then-chiral" is D_r(1-D_r), so the RATIO
!      of counts is D_r / (1-D_r), IF both are unbiased in direction.
!      The wall constraint makes the actual ratio deviate below this,
!      but the log-log slope at large D_r should still approach 1).
!   - All four L curves should overlap in the plateau region (same
!     master curve, no finite-size correction at small D_r); the
!     divergence at D_r → 1 is universal too.
!
! Convention / kernel
! -------------------
!   - Lattice IS the L × L box, OBC, no wall ring.
!   - We use tcrw_step_mask (the kernel from Fig 2 that reports
!     step_type) with an all-TRUE mask, which reduces exactly to
!     the OBC dynamics of tcrw_step_obc while exposing step_type.
!   - Same adaptive τ_relax = max(L^2/D_r, 1/D_r^2) as Fig 3(a).
!
! Parameters
! ----------
!   ω             = 1.0                        ! fixed (paper Fig 3b)
!   L grid        = (4, 9, 19, 49)             ! paper subset, matches Fig 3a
!   D_r grid      = 25 log-spaced points, 10^-4 → 10^0 (exact paper grid)
!   T_floor       = 10^8                       ! minimum measurement length
!   N_burn_floor  = 10^7                       ! minimum discarded burn-in
!   K_meas        = 100                        ! T_use   = max(T_floor,    K_meas · τ_relax)
!   K_burn        = 10                         ! N_burn  = max(N_burn_floor, K_burn · τ_relax)
!   seed          = 20260421                   ! distinct from Fig 3(a), 3(f)
!
! Cost
! ----
!   Same budget as Fig 3(a): worst cell (L = 49, D_r = 10^-4) has
!   τ_relax = 1/D_r^2 = 10^8, so T_use = 10^10 there, and N_burn = 10^9.
!   Full sweep ~ 25 min on a modern Mac (same ballpark as Fig 3a).
!
! Output
! ------
!   tcrw_fig3b_summary.txt
!     # columns:  L   D_r   ratio   |J_Dr|_wall   |J_ω|_wall   Jx_w_Dr   Jy_w_Dr   Jx_w_om   Jy_w_om
!     one row per (L, D_r); 100 rows total for the default L grid.
!
! Sanity checks (printed to stdout)
! ---------------------------------
!   - At D_r = 1 the ratio should go to infinity (|J_ω|_wall → 0 since
!     all translations are post-noise).  Numerically we'll see a huge
!     value; the log-log plot should climb steeply.
!   - At D_r = 10^-3 the plateau value should be O(1) and the same for
!     all four L.
!   - Jy_w_om should be NEGATIVE and dominate |J_ω|_wall (CCW edge
!     current flows down along the left wall).
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig3b.f90 -o tcrw_fig3b
! Run   :  ./tcrw_fig3b
! Plot  :  bash run.sh tcrw-plot-fig3b
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig3b
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- Fig 3(b) parameters ----
   real(dp), parameter :: omega    = 1.0_dp
   integer,  parameter :: L_list(4) = (/ 4, 9, 19, 49 /)
   integer,  parameter :: n_Dr    = 25
   real(dp), parameter :: log_Dr_min = -4.0_dp
   real(dp), parameter :: log_Dr_max =  0.0_dp
   integer(i8), parameter :: T_floor      = 100000000_i8   ! 10^8
   integer(i8), parameter :: N_burn_floor =  10000000_i8   ! 10^7
   real(dp),    parameter :: K_meas       = 100.0_dp
   real(dp),    parameter :: K_burn       =  10.0_dp
   integer,     parameter :: seed         = 20260421

   ! ---- locals ----
   integer     :: iL, iD, u_sum
   integer     :: L_cur
   real(dp)    :: D_r_cur, abs_JDr, abs_Jom, ratio
   real(dp)    :: D_r_values(n_Dr)
   real(dp)    :: t0, t1, t_run
   integer(i8) :: Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om

   call sgrnd(seed)

   ! ---- build log-spaced D_r grid (same as Fig 3a) ----
   call build_log_grid(D_r_values, n_Dr, log_Dr_min, log_Dr_max)

   print '(A)',         '==== TCRW Fig 3(b): |J_Dr|_wall / |J_ω|_wall vs D_r (ω = 1, OBC) ===='
   print '(A,I12,A,I12)',       '  T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   print '(A,F6.1,A,F6.1,A)', &
        '  T_use      = max(T_floor,      ', K_meas, ' * max(L^2, 1/D_r)/D_r)    ' // &
        '  N_burn_use = max(N_burn_floor, ', K_burn, ' * max(L^2, 1/D_r)/D_r)'
   print '(A,I0,A,4(1X,I0))',   '  L grid (', size(L_list), ') :', L_list
   print '(A,I0,A,ES9.2,A,ES9.2)', &
        '  D_r grid (', n_Dr, ' log-spaced):  ', D_r_values(1), ' ... ', D_r_values(n_Dr)
   print '(A,I0)',              '  seed  = ', seed
   print '(A)',                 ''

   open(newunit=u_sum, file='tcrw_fig3b_summary.txt', status='replace', action='write')
   write(u_sum, '(A)') '# TCRW Fig 3(b)  |J_Dr|_wall / |J_ω|_wall vs D_r  (ω = 1, OBC)'
   write(u_sum, '(A,I12,A,I12)') '# T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   write(u_sum, '(A,F6.1,A,F6.1,A)') &
        '# T_use = max(T_floor, ', K_meas, '*max(L^2,1/D_r)/D_r);  ' // &
        'N_burn_use = max(N_burn_floor, ', K_burn, '*max(L^2,1/D_r)/D_r)'
   write(u_sum, '(A,I0)') '# seed    = ', seed
   write(u_sum, '(A)') '# columns: L   D_r   ratio   |J_Dr|_wall   |J_ω|_wall   ' // &
                        'Jx_w_Dr   Jy_w_Dr   Jx_w_om   Jy_w_om'

   ! ---- outer loop over L ----
   do iL = 1, size(L_list)
      L_cur = L_list(iL)

      print '(A,I3,A)', '  --- L = ', L_cur, ' ---'
      print '(A)', '       D_r        |J_Dr|_wall     |J_ω|_wall         ratio      Jy_om_sign  cpu[s]'
      print '(A)', '   -----------  --------------  --------------  --------------  ----------  ------'

      ! ---- inner loop over D_r ----
      do iD = 1, n_Dr
         D_r_cur = D_r_values(iD)
         call cpu_time(t0)
         call run_one(omega, D_r_cur, L_cur, T_floor, N_burn_floor, &
                      Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om,           &
                      abs_JDr, abs_Jom, ratio)
         call cpu_time(t1)
         t_run = t1 - t0

         print '(2X,ES12.4, 3(2X,ES13.5), 2X,I6,6X,F8.2)', &
              D_r_cur, abs_JDr, abs_Jom, ratio, sign_of(Jy_w_om), t_run

         write(u_sum, '(I3, 1X, ES13.5, 1X, ES13.5, 1X, ES13.5, 1X, ES13.5, 4(1X, I14))') &
              L_cur, D_r_cur, ratio, abs_JDr, abs_Jom, &
              Jx_w_Dr, Jy_w_Dr, Jx_w_om, Jy_w_om
      end do
      print '(A)', ''
   end do

   close(u_sum)
   print '(A,I0,A)', 'Wrote summary -> tcrw_fig3b_summary.txt   (', &
                     size(L_list) * n_Dr, ' rows)'
   print '(A)',      'Plot with:    bash run.sh tcrw-plot-fig3b'

contains

   !------------------------------------------------------------------
   ! Small helper: sign of an int64 (-1, 0, +1).  Used to verify that
   ! J_y on the left wall is NEGATIVE (CCW edge current slides walker
   ! downward on the left wall) at ω = 1.
   !------------------------------------------------------------------
   integer function sign_of(n) result(s)
      integer(i8), intent(in) :: n
      if      (n > 0_i8) then ; s =  1
      else if (n < 0_i8) then ; s = -1
      else                    ; s =  0
      end if
   end function sign_of

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
   ! Run a single (ω, D_r, L) MC trajectory of length T_use steps
   ! (after N_burn_use discarded burn-in steps), accumulating the
   ! TWO left-wall currents J_Dr and J_ω decomposed by the
   ! "prev_was_noise" rule.
   !
   ! Timescale:
   !   τ_bulk    = L^2 / D_r                bulk random-walk mixing
   !   τ_wall    = 1 / D_r^2                ω=1 wall-escape trap
   !   τ_relax   = max(τ_bulk, τ_wall)
   !   N_burn_use = max(N_burn_floor_in, K_burn · τ_relax)
   !   T_use      = max(T_floor_in,      K_meas · τ_relax)
   !
   ! Uses tcrw_step_mask with an all-true mask (equivalent to plain
   ! OBC but exposes step_type so we can track prev_noise).
   !
   ! Returns:
   !   Jx_w_Dr, Jy_w_Dr = sum over left-wall (x=0) sites of the
   !                       D_r-attributed translations' (Δx, Δy)
   !   Jx_w_om, Jy_w_om = same, for ω-attributed translations
   !   abs_JDr          = sqrt(Jx_w_Dr^2 + Jy_w_Dr^2)
   !   abs_Jom          = sqrt(Jx_w_om^2 + Jy_w_om^2)
   !   ratio            = abs_JDr / abs_Jom       (huge if abs_Jom == 0)
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
         prev_noise = (step_type == 0)
      end do

      ! ---- measurement ----
      Jx_w_Dr = 0_i8
      Jy_w_Dr = 0_i8
      Jx_w_om = 0_i8
      Jy_w_om = 0_i8

      do it = 1_i8, T_use
         x_before = x          ! y_before not needed (we only filter on x_before == 0)
         d_before = d
         call tcrw_step_mask(x, y, d, mask, L_in, L_in, omega_in, D_r_in, step_type)

         if (step_type == 1 .and. x_before == 0) then
            ! A successful chiral translation departing from the left wall.
            if (prev_noise) then
               Jx_w_Dr = Jx_w_Dr + int(DX(d_before), i8)
               Jy_w_Dr = Jy_w_Dr + int(DY(d_before), i8)
            else
               Jx_w_om = Jx_w_om + int(DX(d_before), i8)
               Jy_w_om = Jy_w_om + int(DY(d_before), i8)
            end if
         end if

         prev_noise = (step_type == 0)
      end do

      deallocate(mask)

      abs_JDr = sqrt( real(Jx_w_Dr, dp)**2 + real(Jy_w_Dr, dp)**2 )
      abs_Jom = sqrt( real(Jx_w_om, dp)**2 + real(Jy_w_om, dp)**2 )

      if (abs_Jom > 0.0_dp) then
         ratio = abs_JDr / abs_Jom
      else
         ratio = huge(1.0_dp)        ! no ω-attributed translations: ratio infinite
      end if
   end subroutine run_one

   !------------------------------------------------------------------
   ! Shared walker-step kernels.
   !------------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig3b

include 'mt.f90'
