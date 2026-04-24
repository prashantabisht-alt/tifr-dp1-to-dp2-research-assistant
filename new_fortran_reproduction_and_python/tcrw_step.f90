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


!---------------------------------------------------------------------
! tcrw_step_mask —  OBC with a generic wall mask (supports defects)
!
!   mask(0:Lx-1, 0:Ly-1) : .true.  → site is allowed (walker may stand here)
!                         .false. → site is a wall  (walker cannot enter)
!
!   step_type (out):
!      0 = noise step         (director rotated, no translation attempt)
!      1 = chiral success     (walker translated by one lattice unit in
!                              the OLD direction, then director rotated)
!      2 = chiral blocked     (target off-grid or mask=.false.;
!                              entire step skipped, director unchanged)
!
! Intended use for Fig 2:  the caller saves d BEFORE the call; if
! step_type == 1 the direction of the translation is the saved d
! (since chiral-step rule is translate-then-rotate).
!
! The distinction between step_type 0 and 2 lets the caller maintain the
! "prev_was_noise" flag that drives the J = J_ω + J_Dr decomposition in
! Osat et al. Fig 2(c)-(e): a translation at step t is attributed to
!   J_Dr  if the last EFFECTIVE event was a noise step   (prev_was_noise = .true.)
!   J_ω   if the last EFFECTIVE event was a chiral move  (prev_was_noise = .false.)
!
! "Effective" means non-blocked.  A blocked chiral attempt (step_type == 2)
! does NOT reset prev_was_noise — the flag keeps whatever it was before
! the blocked attempt.  Only step_type == 0 sets it to .true.; only
! step_type == 1 sets it to .false.  This matches the reference Python
! implementation (ChiralWalker.step in TRW._original_code_by_paperauthors.py),
! which updates self.noise_step only inside the successful branches.
! Physical intuition: while bouncing on a wall the walker is still in the
! "just-scattered" regime; it has not performed a new chiral event yet.
!
! Caller loop pattern (required to reproduce the paper's J_Dr / J_ω curves):
!     select case (step_type)
!        case (0); prev_noise = .true.
!        case (1); prev_noise = .false.
!        case (2); continue                    ! blocked — leave flag alone
!     end select
! Equivalently:  if (step_type /= 2) prev_noise = (step_type == 0)
!---------------------------------------------------------------------
subroutine tcrw_step_mask(x, y, d, mask, Lx, Ly, omega, D_r, step_type)
   implicit none
   integer,  intent(inout) :: x, y, d
   logical,  intent(in)    :: mask(0:, 0:)
   integer,  intent(in)    :: Lx, Ly
   real(dp), intent(in)    :: omega, D_r
   integer,  intent(out)   :: step_type
   real(dp) :: grnd
   real(dp) :: r_step, r_rot
   integer  :: nx, ny
   integer, parameter :: DX(0:3) = (/ 0,  1,  0, -1 /)
   integer, parameter :: DY(0:3) = (/ 1,  0, -1,  0 /)

   r_step = grnd()
   r_rot  = grnd()

   if (r_step < D_r) then
      ! ---- NOISE step ----
      if (r_rot < omega) then
         d = mod(d + 3, 4)           ! CCW
      else
         d = mod(d + 1, 4)           ! CW
      end if
      step_type = 0
   else
      ! ---- CHIRAL step (translate then rotate, honouring mask) ----
      nx = x + DX(d)
      ny = y + DY(d)
      if (nx >= 0 .and. nx < Lx .and. ny >= 0 .and. ny < Ly) then
         if (mask(nx, ny)) then
            x = nx
            y = ny
            if (r_rot < omega) then
               d = mod(d + 1, 4)     ! CW
            else
               d = mod(d + 3, 4)     ! CCW
            end if
            step_type = 1
         else
            step_type = 2             ! wall (defect or boundary site)
         end if
      else
         step_type = 2                ! would leave the L × L grid
      end if
   end if
end subroutine tcrw_step_mask
