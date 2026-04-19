!=====================================================================
! tcrw_fig2_defects.f90 — Fig 2 rows (k)-(o): clean OBC box with an
!                          internal plus-sign defect
!
! Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020v1, 12 Feb 2026.
!
! Goal
! ----
!   The clean box (tcrw_fig2_clean) showed that in a chiral walker:
!     - density is ω-invariant (noise-driven localization)
!     - J_ω and J_Dr nearly cancel on the outer edge
!     - J_Dr has NO normal component on a straight edge (symmetry)
!
!   Adding an internal defect (a cluster of wall cells inside the bulk)
!   does two things the clean box could not show:
!     (i)  The walker hugs the DEFECT walls too — P(x,y) develops a
!          second localization ring around the defect.
!     (ii) Interior walls have a well-defined "normal", so J_Dr now
!          picks up a genuine perpendicular component near the defect.
!          This is what the paper's text actually describes by
!          "J_Dr perpendicular to the edge".
!
! Geometry
! --------
!   L × L lattice, L = 10.  Outer ring = walls (same as clean box).
!   Plus-sign defect centered at (4, 5), total 5 cells:
!
!       y=9  W W W W W W W W W W
!       y=8  W . . . . . . . . W
!       y=7  W . . . . . . . . W
!       y=6  W . . . . D . . . W     ← top arm  (4, 6)
!       y=5  W . . . D D D . . W     ← center + L/R arms (3,5), (4,5), (5,5)
!       y=4  W . . . . D . . . W     ← bot arm  (4, 4)
!       y=3  W . . . . . . . . W
!       y=2  W . . . . . . . . W
!       y=1  W . . . . . . . . W
!       y=0  W W W W W W W W W W
!            0 1 2 3 4 5 6 7 8 9
!
!   Allowed sites: 64 (inner 8x8) − 5 (defect) = 59.
!
! Step kernel
! -----------
!   Unchanged — uses tcrw_step_mask from tcrw_step.f90.  The mask just
!   has more FALSE cells now; the walker rebounds off interior walls
!   exactly the same way as outer walls.
!
! ω
! -
!   Only ω = 0.0 (fully chiral CCW).  The defect-induced J_Dr normal
!   component is a chiral effect — at ω = 0.5 it would be zero by
!   symmetry, so the achiral run adds no information.
!
! Cost
! ----
!   1 ω × 10^10 steps × ~30 ns/step  ≈  5 minutes total on a modern Mac.
!
! Output (5 files)
! ----------------
!   tcrw_fig2_occ_defects.txt       : x  y  P(x,y)         (full L×L)
!   tcrw_fig2_traj_defects.txt      : t  x  y              (first 10^6)
!   tcrw_fig2_Jtot_defects.txt      : x  y  Jx  Jy  |J|    (allowed sites)
!   tcrw_fig2_Jomega_defects.txt    : x  y  Jx  Jy  |J|    (allowed sites)
!   tcrw_fig2_JDr_defects.txt       : x  y  Jx  Jy  |J|    (allowed sites)
!   tcrw_fig2_defects_layout.txt    : x  y                 (defect cells,
!                                                           for gnuplot
!                                                           overlay)
!
! Build : gfortran -O2 -fno-range-check -ffree-line-length-none \
!                  tcrw_fig2_defects.f90 -o tcrw_fig2_defects
! Run   : ./tcrw_fig2_defects
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig2_defects
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- parameters (single source of truth) ----
   integer,     parameter :: L         = 10
   integer(i8), parameter :: T_steps   = 10000000000_i8     ! 10^10 (paper)
   real(dp),    parameter :: D_r       = 1.0e-3_dp
   real(dp),    parameter :: omega     = 0.0_dp             ! fully chiral
   integer,     parameter :: seed      = 20260418
   integer,     parameter :: traj_len  = 1000000            ! first 10^6 steps
   integer,     parameter :: n_prog    = 10                 ! print every 10%

   integer, parameter :: DX(0:3) = (/ 0,  1,  0, -1 /)
   integer, parameter :: DY(0:3) = (/ 1,  0, -1,  0 /)

   ! ---- plus-sign defect (5 cells) ----
   integer, parameter :: n_defects = 5
   integer, parameter :: defect_x(n_defects) = (/ 4, 3, 5, 4, 4 /)
   integer, parameter :: defect_y(n_defects) = (/ 5, 5, 5, 4, 6 /)

   ! ---- walker state ----
   integer     :: x, y, d
   integer     :: x_before, y_before, d_before
   integer(i8) :: it
   integer     :: step_type
   logical     :: prev_noise
   logical     :: mask(0:L-1, 0:L-1)

   ! ---- accumulators ----
   integer(i8) :: occ(0:L-1, 0:L-1)
   integer(i8) :: f_tot(0:3, 0:L-1, 0:L-1)
   integer(i8) :: f_om (0:3, 0:L-1, 0:L-1)
   integer(i8) :: f_dr (0:3, 0:L-1, 0:L-1)

   integer :: traj_x(traj_len), traj_y(traj_len)

   ! ---- misc ----
   real(dp) :: grnd
   integer  :: ix, iy, k, u, n_allowed
   real(dp) :: t0, t1
   integer(i8) :: prog_every
   real(dp) :: P_outer_sum, P_defect_ring_sum, P_bulk_sum
   integer  :: n_outer, n_defect_ring, n_bulk

   ! -----------------------------------------------------------------
   ! Build mask: outer ring walls + plus-sign defect
   ! -----------------------------------------------------------------
   mask = .false.
   do iy = 1, L - 2
      do ix = 1, L - 2
         mask(ix, iy) = .true.
      end do
   end do
   do k = 1, n_defects
      mask(defect_x(k), defect_y(k)) = .false.
   end do

   n_allowed = count(mask)

   call sgrnd(seed)

   print '(A)', '==== TCRW Fig 2 defects (plus-sign at (4,5)) ===='
   print '(A,I0,A,I0,A,I0,A)', '  grid: ', L, ' x ', L, &
                               '  allowed sites = ', n_allowed, &
                               '  (inner 8x8 minus 5 defect cells)'
   print '(A,I14,A,ES9.2,A,F4.2)', '  T_steps = ', T_steps, &
                                   '   D_r = ', D_r, &
                                   '   ω = ', omega
   print '(A,I0,A)', '  first ', traj_len, ' steps saved as trajectory'
   print '(A)', '  defect cells (x, y):'
   do k = 1, n_defects
      print '(A,I0,A,I0,A)', '    (', defect_x(k), ', ', defect_y(k), ')'
   end do
   print '(A)', ''

   prog_every = T_steps / int(n_prog, i8)

   ! -----------------------------------------------------------------
   ! Write defect layout (for gnuplot overlays)
   ! -----------------------------------------------------------------
   call write_defect_layout()

   ! -----------------------------------------------------------------
   ! Initialize accumulators + walker
   ! -----------------------------------------------------------------
   occ   = 0_i8
   f_tot = 0_i8
   f_om  = 0_i8
   f_dr  = 0_i8

   ! initial position: (2, 2) — guaranteed allowed (not on defect, not on wall)
   x = 2
   y = 2
   d = int(4.0_dp * grnd())
   if (d == 4) d = 3
   prev_noise = .false.

   ! sanity: if somehow initial position is forbidden, bail loudly
   if (.not. mask(x, y)) then
      print '(A,I0,A,I0,A)', 'ERROR: initial position (', x, ', ', y, ') is a wall/defect. Aborting.'
      stop 1
   end if

   print '(A,I0,A,I0,A,I0,A)', 'Initial state: (x, y, d) = (', x, ', ', y, ', ', d, ')'
   print '(A)', ''

   call cpu_time(t0)

   ! -----------------------------------------------------------------
   ! HOT LOOP — 10^10 steps
   ! -----------------------------------------------------------------
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
            f_dr(d_before, x_before, y_before) = &
                 f_dr(d_before, x_before, y_before) + 1_i8
         else
            f_om(d_before, x_before, y_before) = &
                 f_om(d_before, x_before, y_before) + 1_i8
         end if
      end if

      ! trajectory buffer
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

   ! -----------------------------------------------------------------
   ! Write outputs
   ! -----------------------------------------------------------------
   call write_occupancy()
   call write_trajectory()
   call write_current('Jtot',   f_tot)
   call write_current('Jomega', f_om )
   call write_current('JDr',    f_dr )

   ! -----------------------------------------------------------------
   ! Sanity prints:
   !   outer_edge     = allowed sites adjacent to outer wall
   !   defect_ring    = allowed sites adjacent to a defect cell
   !   bulk           = everything else
   ! -----------------------------------------------------------------
   P_outer_sum       = 0.0_dp;  n_outer       = 0
   P_defect_ring_sum = 0.0_dp;  n_defect_ring = 0
   P_bulk_sum        = 0.0_dp;  n_bulk        = 0

   do iy = 1, L - 2
      do ix = 1, L - 2
         if (.not. mask(ix, iy)) cycle                ! skip defect cells
         if (ix == 1 .or. ix == L-2 .or. iy == 1 .or. iy == L-2) then
            P_outer_sum = P_outer_sum + real(occ(ix, iy), dp)
            n_outer = n_outer + 1
         else if (is_adjacent_to_defect(ix, iy)) then
            P_defect_ring_sum = P_defect_ring_sum + real(occ(ix, iy), dp)
            n_defect_ring = n_defect_ring + 1
         else
            P_bulk_sum = P_bulk_sum + real(occ(ix, iy), dp)
            n_bulk = n_bulk + 1
         end if
      end do
   end do

   P_outer_sum       = P_outer_sum       / (real(T_steps, dp) * real(n_outer, dp))
   P_defect_ring_sum = P_defect_ring_sum / (real(T_steps, dp) * real(max(n_defect_ring, 1), dp))
   P_bulk_sum        = P_bulk_sum        / (real(T_steps, dp) * real(max(n_bulk, 1), dp))

   print '(A)', ''
   print '(A,I3,A,ES11.3)', '  <P>_outer_edge  (n=', n_outer,       ') = ', P_outer_sum
   print '(A,I3,A,ES11.3)', '  <P>_defect_ring (n=', n_defect_ring, ') = ', P_defect_ring_sum
   print '(A,I3,A,ES11.3)', '  <P>_bulk        (n=', n_bulk,        ') = ', P_bulk_sum
   if (P_bulk_sum > 0.0_dp) then
      print '(A,F7.2)', '  outer/bulk       ratio = ', P_outer_sum / P_bulk_sum
      print '(A,F7.2)', '  defect_ring/bulk ratio = ', P_defect_ring_sum / P_bulk_sum
   end if
   print '(A)', ''
   print '(A)', 'Fig 2 defects: wrote 6 output files (5 data + 1 layout).'
   print '(A)', 'Plot with:'
   print '(A)', '  bash run.sh tcrw-plot-fig2-defects-occ'
   print '(A)', '  bash run.sh tcrw-plot-fig2-defects-traj'
   print '(A)', '  bash run.sh tcrw-plot-fig2-defects-currents'

contains

   !------------------------------------------------------------------
   ! is_adjacent_to_defect: returns .true. if (ix, iy) is an ALLOWED
   ! site with at least one defect cell as an NN neighbor.
   !------------------------------------------------------------------
   function is_adjacent_to_defect(ix, iy) result(adj)
      integer, intent(in) :: ix, iy
      logical :: adj
      integer :: k_, nx, ny

      adj = .false.
      do k_ = 0, 3
         nx = ix + DX(k_)
         ny = iy + DY(k_)
         if (nx >= 0 .and. nx <= L-1 .and. ny >= 0 .and. ny <= L-1) then
            if (is_defect(nx, ny)) then
               adj = .true.
               return
            end if
         end if
      end do
   end function is_adjacent_to_defect


   function is_defect(ix, iy) result(isd)
      integer, intent(in) :: ix, iy
      logical :: isd
      integer :: k_

      isd = .false.
      do k_ = 1, n_defects
         if (ix == defect_x(k_) .and. iy == defect_y(k_)) then
            isd = .true.
            return
         end if
      end do
   end function is_defect


   !------------------------------------------------------------------
   ! write_defect_layout: dumps (x, y) of each defect cell so gnuplot
   ! can overlay filled squares showing the defect geometry.
   !------------------------------------------------------------------
   subroutine write_defect_layout()
      integer :: u_, k_

      open(newunit=u_, file='tcrw_fig2_defects_layout.txt', &
           status='replace', action='write')
      write(u_, '(A)') '# TCRW Fig 2 defect cell layout'
      write(u_, '(A,I0)') '# n_defects = ', n_defects
      write(u_, '(A)') '# columns: x  y'
      do k_ = 1, n_defects
         write(u_, '(I4,1X,I4)') defect_x(k_), defect_y(k_)
      end do
      close(u_)
      print '(A)', '    wrote tcrw_fig2_defects_layout.txt'
   end subroutine write_defect_layout


   !------------------------------------------------------------------
   ! write_occupancy: P(x,y) on the FULL L × L grid.  Outer walls
   ! AND defect cells have P = 0 (walker never there).
   !------------------------------------------------------------------
   subroutine write_occupancy()
      integer :: ix_, iy_, u_

      open(newunit=u_, file='tcrw_fig2_occ_defects.txt', &
           status='replace', action='write')
      write(u_, '(A)') '# TCRW Fig 2 occupancy P(X,Y) with plus-sign defect'
      write(u_, '(A,F4.2,A,ES10.3,A,I14)') &
         '# omega = ', omega, '   D_r = ', D_r, '   T = ', T_steps
      write(u_, '(A)') '# columns: x   y   P(x,y)'
      do ix_ = 0, L - 1
         do iy_ = 0, L - 1
            write(u_, '(I4,1X,I4,1X,ES15.6)') &
                 ix_, iy_, real(occ(ix_, iy_), dp) / real(T_steps, dp)
         end do
         write(u_, '(A)') ''          ! blank line between x-columns (pm3d-friendly)
      end do
      close(u_)
      print '(A)', '    wrote tcrw_fig2_occ_defects.txt'
   end subroutine write_occupancy


   !------------------------------------------------------------------
   ! write_trajectory: (t, x, y) for first traj_len steps.
   !------------------------------------------------------------------
   subroutine write_trajectory()
      integer :: k_, u_

      open(newunit=u_, file='tcrw_fig2_traj_defects.txt', &
           status='replace', action='write')
      write(u_, '(A,F4.2)') '# TCRW Fig 2 trajectory with defects, omega = ', omega
      write(u_, '(A)')      '# columns: t   x   y'
      do k_ = 1, traj_len
         write(u_, '(I10,1X,I4,1X,I4)') k_, traj_x(k_), traj_y(k_)
      end do
      close(u_)
      print '(A)', '    wrote tcrw_fig2_traj_defects.txt'
   end subroutine write_trajectory


   !------------------------------------------------------------------
   ! write_current: vector field on ALLOWED sites only.
   ! tag ∈ {'Jtot', 'Jomega', 'JDr'}
   !------------------------------------------------------------------
   subroutine write_current(tag, fcount)
      character(*), intent(in) :: tag
      integer(i8),  intent(in) :: fcount(0:, 0:, 0:)

      character(len=80) :: fname
      integer  :: ix_, iy_, u_
      real(dp) :: Jx, Jy, Jmag, T_inv

      write(fname, '(A,A,A)') 'tcrw_fig2_', trim(tag), '_defects.txt'
      open(newunit=u_, file=trim(fname), status='replace', action='write')
      write(u_, '(A,A,A,F4.2)') &
           '# TCRW Fig 2 current ', trim(tag), ' defects, omega = ', omega
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
   ! Shared walker-step kernel
   ! ----------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig2_defects

include 'mt.f90'
