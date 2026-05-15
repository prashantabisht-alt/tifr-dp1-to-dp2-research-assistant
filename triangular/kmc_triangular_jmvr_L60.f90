!=====================================================================
! kmc_triangular_jmvr_L60.f90  —  Kinetic Monte Carlo for the original-style
! the JMVR-style active random walker on a triangular lattice.
!
! Implements master equations B2-B7 of the Confinement_enhanced_clustering
! draft (Mandal-Barma-Ramola 2021) directly, with NO Bloch matrix shortcuts.
! Output is the histogram of walker final positions at t = t_final on a larger
! L x L lattice torus.  L=60 avoids visible wrap-around at t=50 and reproduces
! the annular/hexagonal visual style of the Confinement draft Fig. 11.
!
! Used to verify the sign-error bug found in coeff3 of rtp_tl_2.nb.
!
! Model
! -----
!   Walker on triangular lattice with primitive vectors a1 = (2a, 0),
!   a2 = (a, b) -> uses lattice coordinates (n1, n2) with PBC mod L.
!   Internal state m in {0..5} corresponds to 6 NN directions:
!     0: +a1 = (+1, 0) in lattice
!     1: +a2 = (0, +1)
!     2: -a1 + a2 = (-1, +1)
!     3: -a1 = (-1, 0)
!     4: -a2 = (0, -1)
!     5: +a1 - a2 = (+1, -1)
!   Rates (continuous time):
!     translation rate to direction d when in state m:
!         forward (d=m):              1/6 + epsilon
!         backward (d=(m+3) mod 6):   1/6 - epsilon
!         sideways (other 4):         1/6
!     state-flip rate (rotation):
!         m -> m+1 mod 6:  gamma/2
!         m -> m-1 mod 6:  gamma/2
!   Total outgoing rate per walker = 1 + gamma.
!
! Build:
!   gfortran -O3 -fno-range-check -ffree-line-length-none \
!            kmc_triangular_jmvr.f90 -o kmc_triangular
!
! Run:
!   ./kmc_triangular_L60 > kmc_L60_run.log
!
! Output:
!   kmc_triangular_counts_L60.txt   columns:  n1  n2  count
!
! Author: Prashant Bisht, TIFR Hyderabad
!=====================================================================
program kmc_triangular_jmvr
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   integer, parameter :: i8 = selected_int_kind(18)

   ! ---- Physical parameters (match Confinement draft Fig 11) ----
   real(dp), parameter   :: gamma         = 0.01_dp
   real(dp), parameter   :: epsilon_param = 0.15_dp
   real(dp), parameter   :: t_final       = 50.0_dp
   integer,  parameter   :: L             = 60
   integer(i8), parameter :: n_walkers    = 20000000_i8       ! 20 M, enough for visual comparison
   integer,  parameter   :: seed          = 20260515

   ! ---- NN lattice displacements for each director ----
   integer, parameter :: NN_n1(0:5) = (/ +1,  0, -1, -1,  0, +1 /)
   integer, parameter :: NN_n2(0:5) = (/  0, +1, +1,  0, -1, -1 /)

   ! ---- Hop cumulative probabilities (precomputed per director) ----
   real(dp) :: hop_cum(0:5, 0:5)    ! hop_cum(m, d) = P(direction <= d | state m)

   ! ---- Walker state and counters ----
   integer     :: n1, n2, m
   integer(i8) :: i_walker
   integer(i8) :: counts(0:L-1, 0:L-1)
   real(dp)    :: t, dt, u, rate_total
   real(dp)    :: r_event, r_dir, r_rot
   integer     :: d_chosen, mm, dd
   real(dp)    :: r_rate
   real(dp)    :: cpu_start, cpu_end
   integer     :: i, j

   ! External RNG (Mersenne Twister from mt.f90)
   real(dp) :: grnd
   external :: grnd, sgrnd

   ! ---- Build hop_cum table ----
   ! hop_cum(m, d) = sum of rates for state m, directions 0..d, divided by 1
   ! (since total translation rate = 1 always)
   do mm = 0, 5
      do dd = 0, 5
         if (dd == mm) then
            r_rate = 1.0_dp/6 + epsilon_param      ! forward
         else if (dd == mod(mm + 3, 6)) then
            r_rate = 1.0_dp/6 - epsilon_param      ! backward
         else
            r_rate = 1.0_dp/6                       ! sideways
         end if
         if (dd == 0) then
            hop_cum(mm, dd) = r_rate
         else
            hop_cum(mm, dd) = hop_cum(mm, dd-1) + r_rate
         end if
      end do
   end do

   call sgrnd(seed)
   counts = 0_i8
   rate_total = 1.0_dp + gamma

   write(*, '(A)') '==== Triangular JMVR KMC ===='
   write(*, '(A,F8.4,A,F8.4)') '  gamma   = ', gamma, '   epsilon = ', epsilon_param
   write(*, '(A,F8.2,A,I4)')   '  t_final = ', t_final, '   L = ', L
   write(*, '(A,I12)')         '  walkers = ', n_walkers
   write(*, '(A,I0)')          '  seed    = ', seed
   write(*, '(A)')             ''

   call cpu_time(cpu_start)

   ! ──────────────────────────────────────────────────────────────────
   ! Main walker loop
   ! ──────────────────────────────────────────────────────────────────
   do i_walker = 1_i8, n_walkers
      ! Walker initial state: (0, 0) with uniform random director
      n1 = 0
      n2 = 0
      m  = int(6.0_dp * grnd())
      if (m == 6) m = 5
      t = 0.0_dp

      do
         ! Exponential time to next event
         u  = grnd()
         dt = -log(u) / rate_total
         if (t + dt >= t_final) exit
         t = t + dt

         ! Choose event: rotation vs translation
         r_event = grnd()
         if (r_event < gamma / rate_total) then
            ! ----- Rotation -----
            r_rot = grnd()
            if (r_rot < 0.5_dp) then
               m = mod(m + 1, 6)
            else
               m = mod(m + 5, 6)            ! m - 1 mod 6
            end if
         else
            ! ----- Translation: pick direction by cumulative probabilities -----
            r_dir = grnd()
            d_chosen = 5                     ! default last direction
            do dd = 0, 5
               if (r_dir <= hop_cum(m, dd)) then
                  d_chosen = dd
                  exit
               end if
            end do
            n1 = modulo(n1 + NN_n1(d_chosen), L)
            n2 = modulo(n2 + NN_n2(d_chosen), L)
         end if
      end do

      ! Bin the final position
      counts(n1, n2) = counts(n1, n2) + 1_i8

      ! Progress every 10%
      if (mod(i_walker, n_walkers / 10) == 0_i8) then
         call cpu_time(cpu_end)
         write(*, '(A,I12,A,F8.1,A)') &
              '  walkers done: ', i_walker, '   cpu = ', cpu_end - cpu_start, ' s'
      end if
   end do

   call cpu_time(cpu_end)
   write(*, '(A,F8.1,A)') '  ALL DONE in ', cpu_end - cpu_start, ' s'

   ! ──────────────────────────────────────────────────────────────────
   ! Write output
   ! ──────────────────────────────────────────────────────────────────
   open(unit=10, file='kmc_triangular_counts_L60.txt', status='replace', action='write')
   write(10, '(A,F8.4,A,F8.4,A,F8.2,A,I4,A,I12)') &
        '# gamma=', gamma, ' epsilon=', epsilon_param, ' t=', t_final, &
        ' L=', L, ' N=', n_walkers
   write(10, '(A)') '# columns:  n1  n2  count'
   do i = 0, L-1
      do j = 0, L-1
         write(10, '(I4,1X,I4,1X,I14)') i, j, counts(i, j)
      end do
   end do
   close(10)

   write(*, '(A)') 'Saved kmc_triangular_counts_L60.txt'

   ! Sanity check
   write(*, '(A,I14,A,I12)') '  Sum counts = ', sum(counts), '   (should equal n_walkers = ', n_walkers, ')'

end program kmc_triangular_jmvr

include 'mt.f90'
