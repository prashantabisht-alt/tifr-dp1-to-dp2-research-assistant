!=====================================================================
! tcrw_fig2_clean.f90 — Fig 2 rows (a)–(e) and (f)–(j): clean OBC box
!
! Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020v1, 12 Feb 2026.
! Panels reproduced:
!   Row 1 (a–e): ω = 0.0  (fully chiral CCW)  — localization + CCW edge current
!   Row 2 (f–j): ω = 0.5  (achiral)          — localization, no directed current
!   Row 3       : ω = 1.0  (fully chiral CW)  — mirror of ω=0, symmetry check
!                                              (NOT a paper panel, added as
!                                               a clean-chirality sanity test)
!   (Row with internal defects is a separate driver.)
!
! System
! ------
!   L × L lattice, L = 10.  OBC: the outer ring of lattice points
!   (x=0, x=L-1, y=0, y=L-1) are walls — the walker lives on the inner
!   (L-2) × (L-2) = 8 × 8 playground.  When a chiral step targets a
!   wall site the attempt is skipped (director unchanged).  This matches
!   the paper's definition: "an edge or a boundary is defined as
!   lattice points along the boundary not allowed to be occupied by
!   the walker".
!
! Why one driver for two ω
! ------------------------
!   Both ω ∈ {0, 0.5} share the identical mask and step kernel — the
!   physics difference is purely in the step rule.  Packaging both runs
!   here keeps the parameter block a single source of truth and makes
!   the ω=0 vs ω=0.5 comparison trivially reproducible.
!
! Accumulators
! ------------
!   occ  (L,L)             visits per site.  P(X,Y) = occ / T_steps.
!   f_tot(4,L,L)           translation count per direction per site (total)
!   f_om (4,L,L)           translations whose PREVIOUS step was chiral
!   f_dr (4,L,L)           translations whose PREVIOUS step was noise
!
!   Net current at (X,Y):
!      J_x(X,Y) = (f(right,X,Y) - f(left,X,Y))  / T_steps
!      J_y(X,Y) = (f(up,X,Y)    - f(down,X,Y))  / T_steps
!   A translation leaves the source site (x_before, y_before) in
!   direction d_before (since chiral = translate-then-rotate).
!
! J = J_ω + J_Dr decomposition  (paper eqn (III)-ish)
! ---------------------------------------------------
!   "If a translational move occurs immediately after a noise step then
!    it contributes to J_Dr, otherwise it is part of J_ω."
!   We classify each translation by the boolean `prev_noise`:
!      prev_noise = .true.  → flow goes into f_dr
!      prev_noise = .false. → flow goes into f_om
!   `prev_noise` is updated to .true. iff the just-completed step was a
!   noise step (step_type == 0); any chiral attempt (success OR wall-
!   blocked) sets it to .false. — that's consistent with the paper's
!   "otherwise" clause.
!
! Cost
! ----
!   3 ω × 10¹⁰ steps × ~30 ns/step  ≈  15 minutes total on a modern Mac.
!   This matches the paper's T = 10¹⁰.  At 10⁹ we previously saw the
!   mirror-symmetry check (ω=0 vs ω=1) only agree at the 10% level, with
!   residuals dominated by the quasi-cancellation J_tot = J_ω + J_Dr.
!   At 10¹⁰ the statistical noise on each component drops by √10 ≈ 3.2,
!   so the residual in J_tot tightens by the same factor, and |ΔP(x,y)|
!   between nominally-mirror runs drops into the 10⁻⁵ range.
!
! Output per ω  ∈ {0.0, 0.5, 1.0}   (5 files × 3 = 15 total)
! -------------------------------
!   tcrw_fig2_occ_w{0.0,0.5,1.0}.txt       : x  y  P(x,y)
!   tcrw_fig2_traj_w{0.0,0.5,1.0}.txt      : t  x  y           (first 10⁶ steps)
!   tcrw_fig2_Jtot_w{0.0,0.5,1.0}.txt      : x  y  Jx  Jy  |J|
!   tcrw_fig2_Jomega_w{0.0,0.5,1.0}.txt    : x  y  Jx  Jy  |J|
!   tcrw_fig2_JDr_w{0.0,0.5,1.0}.txt       : x  y  Jx  Jy  |J|
!
!   Occupancy file writes ALL L × L sites (wall sites have P = 0) so
!   gnuplot draws the full box.  Current files write ONLY allowed
!   (inner) sites so the vector field is clean.
!
! Build : gfortran -O2 -fno-range-check -ffree-line-length-none \
!                  tcrw_fig2_clean.f90 -o tcrw_fig2_clean
! Run   : ./tcrw_fig2_clean
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig2_clean
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- Fig 2 parameters (single source of truth) ----
   integer,     parameter :: L         = 10
   integer(i8), parameter :: T_steps   = 10000000000_i8      ! 10^10 (paper)
   real(dp),    parameter :: D_r       = 1.0e-3_dp
   real(dp),    parameter :: omega_list(3) = (/ 0.0_dp, 0.5_dp, 1.0_dp /)
   integer,     parameter :: seed      = 20260418
   integer,     parameter :: traj_len  = 1000000             ! first 10^6 steps
   integer,     parameter :: n_prog    = 10                  ! print every 10%

   integer, parameter :: DX(0:3) = (/ 0,  1,  0, -1 /)
   integer, parameter :: DY(0:3) = (/ 1,  0, -1,  0 /)

   ! ---- state ----
   integer     :: x, y, d
   integer     :: x_before, y_before, d_before
   integer(i8) :: it
   integer     :: step_type
   logical     :: prev_noise
   logical     :: mask(0:L-1, 0:L-1)

   ! ---- accumulators (reset per ω) ----
   integer(i8) :: occ(0:L-1, 0:L-1)
   integer(i8) :: f_tot(0:3, 0:L-1, 0:L-1)
   integer(i8) :: f_om (0:3, 0:L-1, 0:L-1)
   integer(i8) :: f_dr (0:3, 0:L-1, 0:L-1)

   integer :: traj_x(traj_len), traj_y(traj_len)

   ! ---- misc ----
   real(dp) :: grnd
   integer  :: iw, ix, iy, u
   real(dp) :: omega
   character(len=80) :: fname
   real(dp) :: t0, t1
   integer(i8) :: prog_every
   real(dp) :: P_edge_sum, P_bulk_sum
   integer  :: n_edge, n_bulk

   ! ---- build clean-box mask ----
   mask = .false.
   do iy = 1, L - 2
      do ix = 1, L - 2
         mask(ix, iy) = .true.
      end do
   end do

   call sgrnd(seed)

   print '(A)',                               '==== TCRW Fig 2 (clean OBC box) ===='
   print '(A,I0,A,I0,A)',                     '  grid: ', L, ' x ', L, &
                                              '  (inner 8 x 8 = 64 allowed sites)'
   print '(A,I14,A,ES9.2,A,I0)',              '  T_steps = ', T_steps, &
                                              '   D_r = ', D_r, '   seed = ', seed
   print '(A,I0,A)',                          '  first ', traj_len, ' steps saved as trajectory'
   print '(A)', ''

   prog_every = T_steps / int(n_prog, i8)

   ! =================================================================
   ! outer loop over ω
   ! =================================================================
   do iw = 1, size(omega_list)
      omega = omega_list(iw)

      ! reset accumulators
      occ   = 0_i8
      f_tot = 0_i8
      f_om  = 0_i8
      f_dr  = 0_i8

      ! initial state: center of the inner box, random director
      x = L / 2
      y = L / 2
      d = int(4.0_dp * grnd())
      if (d == 4) d = 3
      prev_noise = .false.

      print '(A,F4.2,A)', '--- ω = ', omega, ' ---'
      call cpu_time(t0)

      ! -------------- hot loop (10^9 steps) --------------
      do it = 1, T_steps
         x_before = x
         y_before = y
         d_before = d

         call tcrw_step_mask(x, y, d, mask, L, L, omega, D_r, step_type)

         ! occupancy: always counted
         occ(x, y) = occ(x, y) + 1_i8

         ! flow: only on successful translation
         if (step_type == 1) then
            f_tot(d_before, x_before, y_before) = &
                 f_tot(d_before, x_before, y_before) + 1_i8
            if (prev_noise) then
               f_dr (d_before, x_before, y_before) = &
                    f_dr(d_before, x_before, y_before) + 1_i8
            else
               f_om (d_before, x_before, y_before) = &
                    f_om(d_before, x_before, y_before) + 1_i8
            end if
         end if

         ! trajectory buffer (first 10^6 steps)
         if (it <= int(traj_len, i8)) then
            traj_x(int(it)) = x
            traj_y(int(it)) = y
         end if

         prev_noise = (step_type == 0)

         if (mod(it, prog_every) == 0_i8) then
            call cpu_time(t1)
            print '(A,I4,A,F8.1,A)', '    progress ', &
                 int(100_i8 * it / T_steps), ' %   cpu = ', t1 - t0, ' s'
         end if
      end do

      call cpu_time(t1)
      print '(A,F8.1,A)', '    loop done    cpu = ', t1 - t0, ' s'

      ! -------- write outputs for this ω --------
      call write_occupancy(omega)
      call write_trajectory(omega)
      call write_current(omega, 'Jtot',   f_tot)
      call write_current(omega, 'Jomega', f_om )
      call write_current(omega, 'JDr',    f_dr )

      ! -------- cheap sanity print --------
      ! edge = allowed sites adjacent to wall (perimeter of inner box)
      ! bulk = strictly interior of inner box
      P_edge_sum = 0.0_dp;  n_edge = 0
      P_bulk_sum = 0.0_dp;  n_bulk = 0
      do iy = 1, L - 2
         do ix = 1, L - 2
            if (ix == 1 .or. ix == L-2 .or. iy == 1 .or. iy == L-2) then
               P_edge_sum = P_edge_sum + real(occ(ix, iy), dp)
               n_edge = n_edge + 1
            else
               P_bulk_sum = P_bulk_sum + real(occ(ix, iy), dp)
               n_bulk = n_bulk + 1
            end if
         end do
      end do
      P_edge_sum = P_edge_sum / (real(T_steps, dp) * real(n_edge, dp))
      P_bulk_sum = P_bulk_sum / (real(T_steps, dp) * real(n_bulk, dp))
      print '(A,ES11.3,A,ES11.3,A,F7.2)', &
           '    <P>_edge = ', P_edge_sum, &
           '   <P>_bulk = ', P_bulk_sum, &
           '   ratio = ',    P_edge_sum / max(P_bulk_sum, 1.0e-30_dp)
      print '(A)', ''
   end do

   print '(A)', 'Fig 2 clean-box: wrote 10 output files (5 per ω).'
   print '(A)', 'Plot with:'
   print '(A)', '  bash run.sh tcrw-plot-fig2-occ        # P(X,Y) heatmaps'
   print '(A)', '  bash run.sh tcrw-plot-fig2-traj       # rainbow trajectories'
   print '(A)', '  bash run.sh tcrw-plot-fig2-currents   # 2x3 vector fields'

contains

   !------------------------------------------------------------------
   ! write P(X,Y) on the FULL L × L grid  (wall sites → P = 0)
   !------------------------------------------------------------------
   subroutine write_occupancy(omega)
      real(dp), intent(in) :: omega
      integer :: ix_, iy_, u_

      write(fname, '(A,F3.1,A)') 'tcrw_fig2_occ_w', omega, '.txt'
      open(newunit=u_, file=trim(fname), status='replace', action='write')
      write(u_, '(A)') '# TCRW Fig 2 occupancy P(X,Y)'
      write(u_, '(A,F4.2,A,ES10.3,A,I14)') &
         '# omega = ', omega, '   D_r = ', D_r, '   T = ', T_steps
      write(u_, '(A)') '# columns: x   y   P(x,y)'
      do ix_ = 0, L - 1
         do iy_ = 0, L - 1
            write(u_, '(I4,1X,I4,1X,ES15.6)') &
                 ix_, iy_, real(occ(ix_, iy_), dp) / real(T_steps, dp)
         end do
         write(u_, '(A)') ''          ! blank line between x-columns → pm3d-friendly
      end do
      close(u_)
      print '(A,A)', '    wrote ', trim(fname)
   end subroutine write_occupancy


   !------------------------------------------------------------------
   ! write trajectory (t, x, y) for the first traj_len steps
   !------------------------------------------------------------------
   subroutine write_trajectory(omega)
      real(dp), intent(in) :: omega
      integer :: k, u_

      write(fname, '(A,F3.1,A)') 'tcrw_fig2_traj_w', omega, '.txt'
      open(newunit=u_, file=trim(fname), status='replace', action='write')
      write(u_, '(A,F4.2)') '# TCRW Fig 2 trajectory for omega = ', omega
      write(u_, '(A)')      '# columns: t   x   y'
      do k = 1, traj_len
         write(u_, '(I10,1X,I4,1X,I4)') k, traj_x(k), traj_y(k)
      end do
      close(u_)
      print '(A,A)', '    wrote ', trim(fname)
   end subroutine write_trajectory


   !------------------------------------------------------------------
   ! write a vector current field on ALLOWED sites only.
   ! tag ∈ {'Jtot', 'Jomega', 'JDr'}
   !------------------------------------------------------------------
   subroutine write_current(omega, tag, fcount)
      real(dp),     intent(in) :: omega
      character(*), intent(in) :: tag
      integer(i8),  intent(in) :: fcount(0:, 0:, 0:)

      integer :: ix_, iy_, u_
      real(dp) :: Jx, Jy, Jmag, T_inv

      write(fname, '(A,A,A,F3.1,A)') 'tcrw_fig2_', tag, '_w', omega, '.txt'
      open(newunit=u_, file=trim(fname), status='replace', action='write')
      write(u_, '(A,A,A,F4.2)') &
           '# TCRW Fig 2 current ', tag, ' for omega = ', omega
      write(u_, '(A,ES10.3,A,I14)') '# D_r = ', D_r, '   T = ', T_steps
      write(u_, '(A)') '# columns:  x  y  Jx  Jy  |J|   (allowed sites only; counts/T_steps)'

      T_inv = 1.0_dp / real(T_steps, dp)
      do ix_ = 0, L - 1
         do iy_ = 0, L - 1
            if (.not. mask(ix_, iy_)) cycle
            Jx = real(fcount(1, ix_, iy_) - fcount(3, ix_, iy_), dp) * T_inv
            Jy = real(fcount(0, ix_, iy_) - fcount(2, ix_, iy_), dp) * T_inv
            Jmag = sqrt(Jx*Jx + Jy*Jy)
            write(u_, '(I4,1X,I4,3(1X,ES14.6))') ix_, iy_, Jx, Jy, Jmag
         end do
      end do
      close(u_)
      print '(A,A)', '    wrote ', trim(fname)
   end subroutine write_current


   ! ----------------------------------------------------------------
   ! Shared walker-step kernels (tcrw_step_unbounded / tcrw_step_obc /
   ! tcrw_step_mask).
   ! ----------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig2_clean

include 'mt.f90'
