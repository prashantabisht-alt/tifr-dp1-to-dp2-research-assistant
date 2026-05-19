!=====================================================================
! kmc_chiral_rtw.f90  —  Kinetic Monte Carlo for the chiral RTW on
! the triangular lattice.
!
! Model (Phase 1 of DP2):
!   Position: triangular-lattice site in axial coordinates (n1, n2).
!   Director: m in {0..5}, labelling the 6 nearest-neighbour directions.
!   Events (continuous-time Gillespie):
!     run:      hop one step along director m,   rate v
!     CCW:      m -> (m+1) mod 6,                rate gamma_plus
!     CW:       m -> (m-1) mod 6,                rate gamma_minus
!     reversal: m -> (m+3) mod 6,                rate gamma_r
!
!   Tumble rates:  gamma_plus  = gamma * (1 + b) / 2
!                  gamma_minus = gamma * (1 - b) / 2
!   where gamma is total adjacent-tumble rate, b is chirality bias.
!
! Outputs:
!   kmc_chiral_msd.dat  — MSD(t) and displacement moments at log-spaced
!                          times, for cross-validation against the Bloch
!                          matrix theory.
!
! Cartesian coordinates:
!   x = n1 + n2/2,   y = sqrt(3) n2 / 2
!   All 6 lattice steps have |Delta|^2_cart = 1.
!
! Build:
!   gfortran -O3 -fno-range-check -ffree-line-length-none \
!            kmc_chiral_rtw.f90 -o kmc_chiral
!
! Author: Prashant Bisht, TIFR Hyderabad (DP2, Phase 1)
!=====================================================================
program kmc_chiral_rtw
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ================================================================
   ! Physical parameters — edit these and recompile
   ! ================================================================
   real(dp), parameter :: v_run     = 1.0_dp       ! run rate
   real(dp), parameter :: gamma_tot = 1.0_dp       ! total adjacent-tumble rate
   real(dp), parameter :: bias      = 0.5_dp       ! chirality bias in [-1,1]
   real(dp), parameter :: gamma_r   = 0.0_dp       ! reversal rate
   real(dp), parameter :: t_final   = 200.0_dp     ! simulation time per walker
   integer(i8), parameter :: n_walkers = 1000000_i8 ! number of independent walkers
   integer,  parameter :: seed      = 20260519

   ! ================================================================
   ! Derived rates
   ! ================================================================
   real(dp), parameter :: gamma_plus  = 0.5_dp * gamma_tot * (1.0_dp + bias)
   real(dp), parameter :: gamma_minus = 0.5_dp * gamma_tot * (1.0_dp - bias)
   real(dp), parameter :: rate_total  = v_run + gamma_tot + gamma_r
   !   (gamma_plus + gamma_minus = gamma_tot)

   ! Cumulative event thresholds (for uniform-draw event selection)
   real(dp), parameter :: cum_run = v_run / rate_total
   real(dp), parameter :: cum_ccw = cum_run + gamma_plus / rate_total
   real(dp), parameter :: cum_cw  = cum_ccw + gamma_minus / rate_total
   !   remainder [cum_cw, 1] -> reversal

   ! ================================================================
   ! NN lattice displacements in axial coords (n1, n2)
   ! ================================================================
   !   m=0: (+1, 0)   m=1: (0,+1)   m=2: (-1,+1)
   !   m=3: (-1, 0)   m=4: (0,-1)   m=5: (+1,-1)
   integer, parameter :: NN_n1(0:5) = (/ +1,  0, -1, -1,  0, +1 /)
   integer, parameter :: NN_n2(0:5) = (/  0, +1, +1,  0, -1, -1 /)

   ! ================================================================
   ! Measurement time grid (log-spaced)
   ! ================================================================
   integer, parameter :: n_times = 60
   real(dp) :: t_meas(n_times)

   ! Accumulators for displacement moments (Cartesian)
   real(dp) :: sum_x(n_times),  sum_y(n_times)
   real(dp) :: sum_x2(n_times), sum_y2(n_times), sum_xy(n_times)
   integer(i8) :: n_samples(n_times)

   ! ================================================================
   ! Walker state variables
   ! ================================================================
   integer  :: n1, n2, m
   real(dp) :: t, dt, u_time, u_event
   integer  :: i_time              ! index of next measurement time
   integer(i8) :: iw               ! walker counter
   real(dp) :: x_cart, y_cart      ! Cartesian position
   real(dp) :: cpu_start, cpu_end
   integer  :: k

   ! Cartesian conversion
   real(dp), parameter :: sqrt3h = 0.86602540378443864676_dp  ! sqrt(3)/2

   ! Theoretical exact D_even and D_odd for comparison
   real(dp) :: alpha_th, beta_th, lam1sq, D_even_th, D_odd_th
   real(dp) :: msd_last, D_even_kmc

   ! External RNG (Mersenne Twister, mt.f90)
   real(dp) :: grnd
   external :: grnd, sgrnd

   ! ================================================================
   ! Setup
   ! ================================================================

   ! Build log-spaced measurement times: t_min to t_final
   !   t_meas(k) = t_min * (t_final / t_min)^((k-1)/(n_times-1))
   do k = 1, n_times
      t_meas(k) = 0.01_dp * (t_final / 0.01_dp) &
                   ** (dble(k - 1) / dble(n_times - 1))
   end do

   ! Zero accumulators
   sum_x  = 0.0_dp
   sum_y  = 0.0_dp
   sum_x2 = 0.0_dp
   sum_y2 = 0.0_dp
   sum_xy = 0.0_dp
   n_samples = 0_i8

   ! Theoretical prediction
   alpha_th = gamma_tot / 2.0_dp + 2.0_dp * gamma_r
   beta_th  = sqrt(3.0_dp) * gamma_tot * bias / 2.0_dp
   lam1sq   = alpha_th**2 + beta_th**2
   D_even_th = v_run**2 * alpha_th / (2.0_dp * lam1sq) + v_run / 4.0_dp
   D_odd_th  = v_run**2 * beta_th  / (2.0_dp * lam1sq)

   call sgrnd(seed)

   write(*, '(A)')      '==== Chiral RTW — KMC on Triangular Lattice ===='
   write(*, '(A,F10.4)') '  v         = ', v_run
   write(*, '(A,F10.4)') '  gamma     = ', gamma_tot
   write(*, '(A,F10.4)') '  bias      = ', bias
   write(*, '(A,F10.4)') '  gamma_r   = ', gamma_r
   write(*, '(A,F10.4)') '  gamma+    = ', gamma_plus
   write(*, '(A,F10.4)') '  gamma-    = ', gamma_minus
   write(*, '(A,F10.4)') '  R_total   = ', rate_total
   write(*, '(A,ES12.4)')'  t_final   = ', t_final
   write(*, '(A,I12)')   '  walkers   = ', n_walkers
   write(*, '(A)')       ''
   write(*, '(A,F12.8)') '  D_even (theory) = ', D_even_th
   write(*, '(A,F12.8)') '  D_odd  (theory) = ', D_odd_th
   write(*, '(A)')       ''

   call cpu_time(cpu_start)

   ! ================================================================
   ! Main walker loop
   ! ================================================================
   do iw = 1_i8, n_walkers

      ! Initialise: origin, uniform random director
      n1 = 0
      n2 = 0
      m  = int(6.0_dp * grnd())
      if (m >= 6) m = 5
      t  = 0.0_dp
      i_time = 1          ! next measurement slot

      ! ---- Gillespie loop for one walker ----
      do
         ! Draw exponential waiting time
         u_time = grnd()
         dt = -log(u_time) / rate_total

         ! Record position at any measurement times we cross
         do while (i_time <= n_times .and. t + dt >= t_meas(i_time))
            x_cart = dble(n1) + 0.5_dp * dble(n2)
            y_cart = sqrt3h * dble(n2)
            sum_x(i_time)  = sum_x(i_time)  + x_cart
            sum_y(i_time)  = sum_y(i_time)  + y_cart
            sum_x2(i_time) = sum_x2(i_time) + x_cart * x_cart
            sum_y2(i_time) = sum_y2(i_time) + y_cart * y_cart
            sum_xy(i_time) = sum_xy(i_time) + x_cart * y_cart
            n_samples(i_time) = n_samples(i_time) + 1_i8
            i_time = i_time + 1
         end do

         ! Stop if we've passed t_final
         if (t + dt >= t_final) exit
         t = t + dt

         ! Draw event type
         u_event = grnd()

         if (u_event < cum_run) then
            ! ---- Run: hop along current director ----
            n1 = n1 + NN_n1(m)
            n2 = n2 + NN_n2(m)

         else if (u_event < cum_ccw) then
            ! ---- Tumble CCW: m -> m+1 ----
            m = mod(m + 1, 6)

         else if (u_event < cum_cw) then
            ! ---- Tumble CW: m -> m-1 ----
            m = mod(m + 5, 6)

         else
            ! ---- Reversal: m -> m+3 ----
            m = mod(m + 3, 6)
         end if
      end do

      ! Any measurement times beyond t_final get the final position
      do while (i_time <= n_times)
         x_cart = dble(n1) + 0.5_dp * dble(n2)
         y_cart = sqrt3h * dble(n2)
         sum_x(i_time)  = sum_x(i_time)  + x_cart
         sum_y(i_time)  = sum_y(i_time)  + y_cart
         sum_x2(i_time) = sum_x2(i_time) + x_cart * x_cart
         sum_y2(i_time) = sum_y2(i_time) + y_cart * y_cart
         sum_xy(i_time) = sum_xy(i_time) + x_cart * y_cart
         n_samples(i_time) = n_samples(i_time) + 1_i8
         i_time = i_time + 1
      end do

      ! Progress report every 10%
      if (mod(iw, max(n_walkers / 10, 1_i8)) == 0_i8) then
         call cpu_time(cpu_end)
         write(*, '(A,I12,A,F8.1,A)') &
              '  done: ', iw, '   cpu = ', cpu_end - cpu_start, ' s'
      end if
   end do

   call cpu_time(cpu_end)
   write(*, '(A,F10.1,A)') '  ALL DONE in ', cpu_end - cpu_start, ' s'
   write(*, '(A)') ''

   ! ================================================================
   ! Write output
   ! ================================================================
   open(unit=10, file='kmc_chiral_msd.dat', status='replace', action='write')
   write(10, '(A,F8.4,A,F8.4,A,F8.4,A,F8.4,A,I0)') &
        '# v=', v_run, ' gamma=', gamma_tot, ' bias=', bias, &
        ' gamma_r=', gamma_r, ' N=', n_walkers
   write(10, '(A,F12.8,A,F12.8)') &
        '# D_even_theory=', D_even_th, ' D_odd_theory=', D_odd_th
   write(10, '(A)') '# columns: t  <x>  <y>  <x^2>  <y^2>  <xy>  MSD  N_samples'
   do k = 1, n_times
      if (n_samples(k) > 0_i8) then
         write(10, '(ES14.6, 5(1X,ES18.10), 1X, ES18.10, 1X, I0)') &
              t_meas(k), &
              sum_x(k)  / dble(n_samples(k)), &
              sum_y(k)  / dble(n_samples(k)), &
              sum_x2(k) / dble(n_samples(k)), &
              sum_y2(k) / dble(n_samples(k)), &
              sum_xy(k) / dble(n_samples(k)), &
              (sum_x2(k) + sum_y2(k)) / dble(n_samples(k)), &
              n_samples(k)
      end if
   end do
   close(10)
   write(*, '(A)') 'Saved kmc_chiral_msd.dat'

   ! ================================================================
   ! Quick D_even extraction from last time point
   ! ================================================================
   if (n_samples(n_times) > 0_i8) then
      msd_last = (sum_x2(n_times) + sum_y2(n_times)) / dble(n_samples(n_times))
      D_even_kmc = msd_last / (4.0_dp * t_meas(n_times))
      write(*, '(A)')       ''
      write(*, '(A)')       '---- Quick D_even check (last time point) ----'
      write(*, '(A,ES14.6)')'  t_last        = ', t_meas(n_times)
      write(*, '(A,F14.6)') '  MSD(t_last)   = ', msd_last
      write(*, '(A,F12.8)') '  D_even (KMC)  = ', D_even_kmc
      write(*, '(A,F12.8)') '  D_even (exact) = ', D_even_th
      write(*, '(A,ES12.4)')'  rel error     = ', &
           abs(D_even_kmc - D_even_th) / D_even_th
   end if

end program kmc_chiral_rtw

include 'mt.f90'
