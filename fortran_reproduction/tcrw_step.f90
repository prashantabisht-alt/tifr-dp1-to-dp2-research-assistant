!=====================================================================
! tcrw_step.f90  —  TCRW walker-update kernel (INCLUDE FILE)
!
! Included inside the `contains` block of any driver program.
! Provides:
!   subroutine tcrw_step_unbounded(x, y, d, omega, D_r)
!   subroutine tcrw_step_obc      (x, y, d, L, omega, D_r)
!
! Requires in host scope:
!   - parameter dp  = selected_real_kind(15, 300)
!   - external     :: grnd   (double-precision RNG from mt.f90)
!
! Direction encoding (matches the Python reference tcrw_core.py):
!   d = 0 → ↑   (Δx, Δy) = ( 0, +1)
!   d = 1 → →   (Δx, Δy) = (+1,  0)
!   d = 2 → ↓   (Δx, Δy) = ( 0, -1)
!   d = 3 → ←   (Δx, Δy) = (-1,  0)
!
! Rotation rule:
!   CCW (angle ↑): d → mod(d + 3, 4)     ! equivalent to d - 1 mod 4
!   CW  (angle ↓): d → mod(d + 1, 4)
!
! Step rule:
!   With probability D_r       — NOISE step (no translation):
!       prob ω   → CCW rotation
!       prob 1-ω → CW  rotation
!   With probability 1 - D_r   — CHIRAL step (translate THEN rotate):
!       prob ω   → CW  rotation
!       prob 1-ω → CCW rotation
!
! OBC rule: if the chiral translation would exit the box, the ENTIRE
! chiral step is skipped (direction unchanged). Noise step is unaffected.
!=====================================================================

subroutine tcrw_step_unbounded(x, y, d, omega, D_r)
   implicit none
   integer,  intent(inout) :: x, y, d
   real(dp), intent(in)    :: omega, D_r
   real(dp) :: grnd                               ! external RNG (mt.f90)
   real(dp) :: r_step, r_rot
   integer, parameter :: DX(0:3) = (/ 0,  1,  0, -1 /)
   integer, parameter :: DY(0:3) = (/ 1,  0, -1,  0 /)

   r_step = grnd()
   r_rot  = grnd()

   if (r_step < D_r) then
      ! ---- NOISE step: rotate only ----
      if (r_rot < omega) then
         d = mod(d + 3, 4)                         ! CCW
      else
         d = mod(d + 1, 4)                         ! CW
      end if
   else
      ! ---- CHIRAL step: translate then rotate (no boundary) ----
      x = x + DX(d)
      y = y + DY(d)
      if (r_rot < omega) then
         d = mod(d + 1, 4)                         ! CW
      else
         d = mod(d + 3, 4)                         ! CCW
      end if
   end if
end subroutine tcrw_step_unbounded


subroutine tcrw_step_obc(x, y, d, L, omega, D_r)
   implicit none
   integer,  intent(inout) :: x, y, d
   integer,  intent(in)    :: L
   real(dp), intent(in)    :: omega, D_r
   real(dp) :: grnd
   real(dp) :: r_step, r_rot
   integer  :: nx, ny
   integer, parameter :: DX(0:3) = (/ 0,  1,  0, -1 /)
   integer, parameter :: DY(0:3) = (/ 1,  0, -1,  0 /)

   r_step = grnd()
   r_rot  = grnd()

   if (r_step < D_r) then
      if (r_rot < omega) then
         d = mod(d + 3, 4)
      else
         d = mod(d + 1, 4)
      end if
   else
      nx = x + DX(d)
      ny = y + DY(d)
      if (nx >= 0 .and. nx < L .and. ny >= 0 .and. ny < L) then
         x = nx
         y = ny
         if (r_rot < omega) then
            d = mod(d + 1, 4)
         else
            d = mod(d + 3, 4)
         end if
      end if
      ! chiral translation blocked → entire step skipped (direction unchanged)
   end if
end subroutine tcrw_step_obc
