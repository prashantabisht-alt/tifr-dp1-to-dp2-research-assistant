!=====================================================================
! tcrw_fig1b.f90 — Fig 1(b) only: sample trajectories
!
! Paper: Osat, Meyberg, Metson, Speck — "Topological Chiral Random Walker"
! Panel 1(b): single long unwrapped trajectory in (x,y), for 4 chiralities.
!
! Parameters (match the Python reference tcrw_phase1_pbc.py):
!   D_r = 1.0e-3                         ! noise/chiral split
!   T   = 1_000_000                      ! total steps per trajectory
!   skip = 100                           ! subsample stride for output
!   ω   ∈ { 0.00, 0.50, 0.75, 1.00 }    ! the four panels overlaid
!
! Output :  tcrw_fig1b_traj_wX.XX.txt   (one file per ω)
! Columns:  step   x   y                 (integer unwrapped positions)
!
! Direction convention (see tcrw_step.f90):
!   d = 0 → ↑   d = 1 → →   d = 2 → ↓   d = 3 → ←
!
! Build :  gfortran -O2 -fno-range-check -ffree-line-length-none \
!                   tcrw_fig1b.f90 -o tcrw_fig1b
! Run   :  ./tcrw_fig1b
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program tcrw_fig1b
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)

   ! ---- Fig 1(b) parameters (single source of truth) ----
   real(dp), parameter :: D_r      = 1.0e-3_dp
   integer,  parameter :: T_steps  = 1000000
   integer,  parameter :: skip     = 1           ! every step — full resolution
                                                 ! ~33 MB per ω file
                                                 ! resolves 4-step CCW/CW loops
                                                 ! paper-facing Python panel.
                                                 ! The dynamics are still run
                                                 ! for all T_steps.
   real(dp), parameter :: omega_list(4) = (/ 0.0_dp, 0.5_dp, 0.75_dp, 1.0_dp /)
   integer,  parameter :: seed     = 123         ! align with Python reference
   integer :: i
   call sgrnd(seed)

   print '(A)', '==== TCRW Fig 1(b): sample trajectories ===='
   print '(A,ES9.2,A,I10,A,I6,A,I8)', &
        '  D_r = ', D_r, '   T = ', T_steps, '   skip = ', skip, '   seed = ', seed

   do i = 1, size(omega_list)
      call run_sample_traj(omega_list(i), D_r, T_steps, skip)
   end do

   print '(A)', ''
   print '(A)', 'Wrote 4 files:  tcrw_fig1b_traj_w{0.00,0.50,0.75,1.00}.txt'
   print '(A)', 'Plot with:      bash run.sh tcrw-plot-fig1b'
   print '(A)', 'PNG preview:    bash run.sh tcrw-plot-fig1b-png'

contains

   !------------------------------------------------------------------
   ! Single long unwrapped trajectory.
   ! Writes (step, x, y) every `skip` steps to tcrw_fig1b_traj_wX.XX.txt.
   !------------------------------------------------------------------
   subroutine run_sample_traj(omega, D_r_in, T_in, skip_in)
      real(dp), intent(in) :: omega, D_r_in
      integer,  intent(in) :: T_in, skip_in
      integer  :: x, y, d, it, u
      real(dp) :: grnd
      character(len=64) :: fname

      x = 0
      y = 0
      d = int(4.0_dp * grnd())               ! uniform in {0,1,2,3}
      if (d == 4) d = 3                       ! guard grnd() ≈ 1.0

      write(fname, '(A,F4.2,A)') 'tcrw_fig1b_traj_w', omega, '.txt'
      open(newunit=u, file=trim(fname), status='replace', action='write')
      write(u, '(A,F5.2,A,ES10.3,A,I10,A,I6)') &
           '# TCRW Fig 1(b) trajectory | omega = ', omega, &
           ' | D_r = ', D_r_in, &
           ' | T = ', T_in,     &
           ' | skip = ', skip_in
      write(u, '(A)') '# columns:  step   x   y'
      write(u, '(I10,2(1X,I10))') 0, x, y

      do it = 1, T_in
         call tcrw_step_unbounded(x, y, d, omega, D_r_in)
         if (mod(it, skip_in) == 0) then
            write(u, '(I10,2(1X,I10))') it, x, y
         end if
      end do

      close(u)
      print '(A,F5.2,A,I10,A,I10,A)',                &
           '   ω = ', omega, ': final (x,y) = (',   &
           x, ',', y, ')'
   end subroutine run_sample_traj


   !------------------------------------------------------------------
   ! Shared walker-step kernels (tcrw_step_unbounded / tcrw_step_obc).
   ! Included inside `contains` so they see `dp` and the external grnd.
   !------------------------------------------------------------------
   include 'tcrw_step.f90'

end program tcrw_fig1b

include 'mt.f90'
