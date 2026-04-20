!=====================================================================
! tcrw_fig3hij.f90 — Fig 3(h), 3(i), 3(j):  left-wall current vector
!                    fields (h) J_Dr  (i) J_ω and (j) θ_JDr, θ_Jω of
!                    the TOTAL left-wall current — ALL vs ω at fixed
!                    D_r = 10^-3, L = 10.
!
! Paper: Osat, Meyberg, Metson, Speck,
!        "Topological chiral random walker"
!        arXiv:2602.12020v1 [cond-mat.stat-mech], 12 Feb 2026.
!
! Panels covered  (one Fortran run feeds all three gnuplot files)
! ---------------------------------------------------------------
!   Fig 3(h):  Quiver of J_Dr(y) on the left wall for each ω.
!              x-axis = ω (linear),  y-axis = edge site index.
!              Paper claim: arrow ROTATES continuously from θ = +π/2
!              (straight UP) at ω = 0, through θ = 0 at ω = 0.5
!              (achirality), to θ = -π/2 (straight DOWN) at ω = 1.
!              This is the ω-sweep twin of Fig 3(c) (which fixed ω = 1
!              and swept D_r).
!   Fig 3(i):  Same, but with J_ω instead of J_Dr.
!              Paper claim: arrow ROTATES from θ = -π/4 (down-right)
!              at ω = 0, through θ = 0 at ω = 0.5, to θ = +π/4
!              (up-right) at ω = 1.  Rigid rotation — magnitude
!              stays roughly constant across ω (unlike |J_Dr|).
!   Fig 3(j):  θ_JDr(ω) and θ_Jω(ω) — angles of the TOTAL left-wall
!              current vectors as functions of ω.
!              θ_JDr:  sigmoid +π/2 → 0 → -π/2  (antisymmetric in ω-½).
!              θ_Jω:   sigmoid -π/4 → 0 → +π/4  (antisymmetric in ω-½).
!
! Why the angle signs flip between ω = 0 and ω = 1
! ------------------------------------------------
!   Rotation encoding  d = 0↑, 1→, 2↓, 3←.
!   Noise rotates CCW with prob (1-ω), CW with prob ω.
!   Chiral rotation is the complement (CW w.p. 1-ω, CCW w.p. ω).
!
!   ω = 1:   noise rotates CW (director CW),  chiral step rotates CCW.
!            Walker jammed at (0, y) facing ←; noise CW → ↓; chiral
!            step translates DOWN.  Each noise event drifts walker
!            DOWN → Jy_Dr < 0   (θ_JDr ≈ -π/2).
!            Orbit departure from x=0 goes  → then ↑ → UP-RIGHT,
!            so J_ω ≈ (+, +)    (θ_Jω ≈ +π/4).
!
!   ω = 0:   noise rotates CCW, chiral step rotates CW.
!            Walker jammed at (0, y) facing ←; noise CCW → ↑; chiral
!            step translates UP.  Each noise event drifts walker
!            UP → Jy_Dr > 0     (θ_JDr ≈ +π/2).
!            Orbit departure from x=0 goes  → then ↓ → DOWN-RIGHT,
!            so J_ω ≈ (+, -)    (θ_Jω ≈ -π/4).
!
!   ω = 0.5: noise and chiral steps are unbiased (50/50).  Signed
!            Jy accumulators have zero mean → θ ≈ 0 (with MC noise
!            that shows up as a small angular wobble near achirality).
!
!   Note: even at ω = 0 and ω = 1, θ_JDr does not reach exactly ±π/2:
!   every now and then a noise event lands the walker into the bulk
!   (so the next chiral step goes sideways rather than along the wall),
!   leaving a small positive Jx_Dr baseline that limits |θ_JDr| to
!   about atan(Jy/Jx) ≈ ±1.46 rad instead of ±π/2 = ±1.571.
!
! Observable per (ω, y) cell   [same definitions as Fig 3(c/d/e)]
! ---------------------------------------------------------------
!   For each ω in the grid and each left-wall site y ∈ {0..L-1}:
!       Jx_Dr(y) = Σ over steps (step_type = 1 AND x_before = 0 AND
!                                y_before = y AND prev_noise) of DX(d_before)
!       Jy_Dr(y) = same, with DY(d_before)
!       Jx_om(y), Jy_om(y) = same but for prev_noise = .false.
!
!       |J_Dr|(y) = √(Jx² + Jy²),    θ_JDr(y)  = atan2(Jy_Dr, Jx_Dr)
!       |J_ω|(y),                    θ_Jω(y)   similarly
!
!   For Fig 3(j) we aggregate over the wall:
!       Jx_Dr_tot = Σ_y Jx_Dr(y), etc.
!       θ_JDr_tot = atan2(Jy_Dr_tot, Jx_Dr_tot)
!       θ_Jω_tot  = atan2(Jy_om_tot,  Jx_om_tot)
!
! Parameters  (paper Fig 3(h/i/j): L = 10, D_r = 10^-3, ω sweep)
! -------------------------------------------------------------
!   L             = 10
!   D_r           = 10^-3                        ! fixed
!   ω grid        = 21 linearly-spaced pts  0.0 → 1.0
!   T_floor       = 3·10^8                       ! bumped from 10^8 for ω=0.5 cleanup
!   N_burn_floor  = 3·10^7                       ! bumped from 10^7
!   K_meas        = 300                          ! bumped from 100 (2026-04-20)
!   K_burn        = 30                           ! bumped from 10
!   seed          = 20260424                     ! distinct from fig3cde
!
! Cost  (updated for bumped parameters)
! -------------------------------------
!   D_r = 10^-3  →  τ_bulk = L²/D_r = 10^5,  τ_wall = 1/D_r² = 10^6,
!                   τ_relax = 10^6.
!   T_use  = max(3·10^8, 300·10^6) = 3·10^8 steps / ω-point.
!   N_burn = max(3·10^7,  30·10^6) = 3·10^7 steps / ω-point.
!   Total: 21 × (3·10^7 + 3·10^8) ≈ 7 × 10^9 RNG-driven steps.  On a
!   modern Mac (~40–60 M steps/s) that's roughly 2–3 min.
!
!   Why the bump?  Original run showed ω = 0.5 column arrows still
!   tilted slightly DOWN in Fig 3(h) (should be ≈ →, since Jy → 0 by
!   symmetry and only noise-scatter Jx survives).  3× more statistics
!   shrinks the residual Jy sample noise enough that atan2(Jy, Jx)
!   lands closer to 0.
!
! Current decomposition (same rule as Fig 2, Fig 3(b), Fig 3(c/d/e/g))
! --------------------------------------------------------------------
!   tcrw_step_mask returns step_type ∈ {0, 1, 2}.
!   prev_noise ← (step_type == 0).
!   At step_type == 1 AND x_before == 0:
!       prev_noise = .true.  → contribute (DX, DY) to J_Dr(y_before)
!       prev_noise = .false. → contribute (DX, DY) to J_ω (y_before)
!
! Output
! ------
!   tcrw_fig3hij_summary.txt   — per-site currents (for (h) and (i)):
!     # iW  ω  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om  |J_Dr|  th_Dr  |J_om|  th_om
!     (L × n_ω rows  =  10 × 21  =  210)
!
!   tcrw_fig3j_summary.txt     — aggregated angles (for (j)):
!     # iW  ω  Jx_Dr_tot  Jy_Dr_tot  Jx_om_tot  Jy_om_tot  th_Dr_tot  th_om_tot
!     (n_ω rows  =  21)
!
! Sanity checks (printed to stdout)
! ---------------------------------
!   - θ_JDr_tot at ω = 0:      should be ≈ +π/2 (modulo the Jx baseline)
!   - θ_JDr_tot at ω = 1:      should be ≈ -π/2 (modulo the Jx baseline)
!   - θ_JDr_tot at ω = 0.5:    should be ≈ 0 (noise dominates Jy)
!   - θ_Jω_tot  at ω = 0:      should be ≈ -π/4
!   - θ_Jω_tot  at ω = 1:      should be ≈ +π/4
!   - θ_Jω_tot  at ω = 0.5:    should be ≈ 0
!   - antisymmetry:  θ(ω) + θ(1-ω) ≈ 0  for both components.
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig3hij.f90 -o tcrw_fig3hij
! Run   :  ./tcrw_fig3hij
! Plot  :  bash run.sh tcrw-plot-fig3h   (quiver J_Dr  vs ω)
!          bash run.sh tcrw-plot-fig3i   (quiver J_ω   vs ω)
!          bash run.sh tcrw-plot-fig3j   (θ vs ω scalar plot)
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig3hij
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- panel parameters ----
   ! Authors' lattice convention: label L means (L+1)×(L+1) sites
   ! (indices 0..L inclusive).  Paper's Fig 3 uses L = 10 ⇒ 11×11 grid.
   ! We keep two variables: L_paper for output labels, L_cur for the sim.
   integer,  parameter :: L_paper   = 10                  ! paper's legend label
   integer,  parameter :: L_cur     = L_paper + 1         ! # sites per side (authors' convention)
   real(dp), parameter :: D_r_fixed = 1.0e-3_dp
   integer,  parameter :: n_omega   = 21
   real(dp), parameter :: omega_min = 0.0_dp
   real(dp), parameter :: omega_max = 1.0_dp
   integer(i8), parameter :: T_floor      = 300000000_i8   ! 3·10^8  (bumped from 10^8)
   integer(i8), parameter :: N_burn_floor =  30000000_i8   ! 3·10^7  (bumped from 10^7)
   real(dp),    parameter :: K_meas       = 300.0_dp       ! bumped 100 → 300 for ω=0.5 column cleanup
   real(dp),    parameter :: K_burn       =  30.0_dp       ! bumped  10 →  30
   integer,     parameter :: seed         = 20260424

   ! ---- locals ----
   integer  :: iW, y, u_field, u_angle
   real(dp) :: omega_cur
   real(dp) :: omega_values(n_omega)
   real(dp) :: t0, t1, t_run

   ! Per-site wall currents for the current ω.  Allocated for L = L_cur.
   integer(i8), allocatable :: Jx_Dr(:), Jy_Dr(:), Jx_om(:), Jy_om(:)
   real(dp)    :: abs_JDr_y, abs_Jom_y, th_JDr_y, th_Jom_y
   integer(i8) :: Jx_Dr_tot, Jy_Dr_tot, Jx_om_tot, Jy_om_tot
   real(dp)    :: th_JDr_tot, th_Jom_tot

   call sgrnd(seed)

   ! ---- build linearly-spaced ω grid ----
   call build_linear_grid(omega_values, n_omega, omega_min, omega_max)

   print '(A)',         '==== TCRW Fig 3(h)/(i)/(j): left-wall J_Dr, J_ω vs ω  (L = 10, D_r = 10^-3) ===='
   print '(A,I0,A,I0,A,ES11.4)', '  L = ', L_paper, ' (sites per side: ', L_cur, ')     D_r (fixed) = ', D_r_fixed
   print '(A,I12,A,I12)',   '  T_floor = ', T_floor, '   N_burn_floor = ', N_burn_floor
   print '(A,F6.1,A,F6.1,A)', &
        '  T_use      = max(T_floor,      ', K_meas, ' * max(L^2, 1/D_r)/D_r)    ' // &
        '  N_burn_use = max(N_burn_floor, ', K_burn, ' * max(L^2, 1/D_r)/D_r)'
   print '(A,I0,A,F4.2,A,F4.2)', &
        '  ω grid  (', n_omega, ' linearly-spaced):  ', &
        omega_values(1), ' ... ', omega_values(n_omega)
   print '(A,I0)',         '  seed  = ', seed
   print '(A)',            ''

   allocate(Jx_Dr(0:L_cur-1), Jy_Dr(0:L_cur-1), Jx_om(0:L_cur-1), Jy_om(0:L_cur-1))

   ! ---- open two output files (per-site field + aggregate angle) ----
   open(newunit=u_field, file='tcrw_fig3hij_summary.txt', status='replace', action='write')
   write(u_field, '(A)') '# TCRW Fig 3(h)/(i): per-site left-wall currents  (L = 10, D_r = 10^-3)'
   write(u_field, '(A,F6.1,A,F6.1,A)') &
        '# T_use = max(T_floor, ', K_meas, '*max(L^2,1/D_r)/D_r);  ' // &
        'N_burn_use = max(N_burn_floor, ', K_burn, '*max(L^2,1/D_r)/D_r)'
   write(u_field, '(A,I0,A,I0,A,ES11.4,A,I0)') '# L = ', L_paper, ' (sites = ', L_cur, ')   D_r = ', D_r_fixed, '   seed = ', seed
   write(u_field, '(A)') '# columns:  iW  ω  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om  ' // &
                         '|J_Dr|  th_Dr  |J_om|  th_om'

   open(newunit=u_angle, file='tcrw_fig3j_summary.txt', status='replace', action='write')
   write(u_angle, '(A)') '# TCRW Fig 3(j): θ_JDr and θ_Jω of the TOTAL left-wall current  (L = 10, D_r = 10^-3)'
   write(u_angle, '(A,I0,A,I0,A,ES11.4,A,I0)') '# L = ', L_paper, ' (sites = ', L_cur, ')   D_r = ', D_r_fixed, '   seed = ', seed
   write(u_angle, '(A)') '# columns:  iW  ω  Jx_Dr_tot  Jy_Dr_tot  Jx_om_tot  Jy_om_tot  ' // &
                         'th_Dr_tot  th_om_tot'

   ! ---- sweep ω ----
   print '(A)', '     iW     ω         Jy_Dr_tot     Jy_om_tot       th_Dr_tot(rad)   th_om_tot(rad)   cpu[s]'
   print '(A)', '   ----  ------      ----------    ----------      ---------------  ---------------   ------'

   do iW = 1, n_omega
      omega_cur = omega_values(iW)
      call cpu_time(t0)
      call run_one(omega_cur, D_r_fixed, L_cur, T_floor, N_burn_floor, &
                   Jx_Dr, Jy_Dr, Jx_om, Jy_om)
      call cpu_time(t1)
      t_run = t1 - t0

      ! ---- write per-site rows to Fig 3(h)/(i) file ----
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

         write(u_field, '(I4, 1X, F8.4, 1X, I3, 4(1X, I14), 1X, ES13.5, 1X, F10.5, 1X, ES13.5, 1X, F10.5)') &
              iW, omega_cur, y, &
              Jx_Dr(y), Jy_Dr(y), Jx_om(y), Jy_om(y), &
              abs_JDr_y, th_JDr_y, abs_Jom_y, th_Jom_y
      end do
      write(u_field, '(A)') ''     ! blank between ω blocks (nice for gnuplot index)

      ! ---- aggregate totals over wall for Fig 3(j) ----
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

      write(u_angle, '(I4, 1X, F8.4, 4(1X, I14), 2(1X, F10.5))') &
           iW, omega_cur, Jx_Dr_tot, Jy_Dr_tot, Jx_om_tot, Jy_om_tot, &
           th_JDr_tot, th_Jom_tot

      print '(2X, I4, 2X, F6.3, 2(2X, I12), 2(2X, F14.5), 2X, F7.2)', &
           iW, omega_cur, Jy_Dr_tot, Jy_om_tot, th_JDr_tot, th_Jom_tot, t_run
   end do

   close(u_field)
   close(u_angle)

   print '(A)', ''
   print '(A,I0,A)', 'Wrote Fig 3(h)/(i) field -> tcrw_fig3hij_summary.txt  (', &
                     L_cur * n_omega, ' site-rows)'
   print '(A,I0,A)', 'Wrote Fig 3(j) angles    -> tcrw_fig3j_summary.txt    (', &
                     n_omega, ' rows)'
   print '(A)',      'Plot:  bash run.sh tcrw-plot-fig3h   # quiver J_Dr vs ω'
   print '(A)',      '       bash run.sh tcrw-plot-fig3i   # quiver J_ω  vs ω'
   print '(A)',      '       bash run.sh tcrw-plot-fig3j   # θ vs ω'

   deallocate(Jx_Dr, Jy_Dr, Jx_om, Jy_om)

contains

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
   ! Run a single (ω, D_r, L) MC trajectory, accumulating PER-SITE
   ! left-wall currents decomposed into J_Dr and J_ω by the standard
   ! prev_noise rule.   (Identical body to tcrw_fig3cde.f90's run_one
   ! with the only difference being which variable is the sweep axis
   ! in the caller.)
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

      ! ---- adaptive T and burn-in (same as Fig 3a/b/cde/g) ----
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

         if (step_type /= 2) prev_noise = (step_type == 0)   ! authors' rule: blocked chiral leaves flag unchanged
      end do

      deallocate(mask)
   end subroutine run_one

   !------------------------------------------------------------------
   ! Shared walker-step kernels.
   !------------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig3hij

include 'mt.f90'
