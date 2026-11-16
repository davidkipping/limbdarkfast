;+
; NAME:
;   LIMBDARKFAST_HYBRID_FLUX_OBLATE
;
; PURPOSE:
;   Normalised flux of a star with hybrid2, hybrid4 or hybrid5 limb
;   darkening, occulted by an oblate planet whose projected outline is an
;   ellipse, in closed form (the IDL version of limbdarkfast).
;
; CALLING SEQUENCE:
;   F = LIMBDARKFAST_HYBRID_FLUX_OBLATE(X, Y, r, f, theta, w [, law])
;
; INPUTS:
;   X, Y  - the planet's sky position in stellar radii, X along its
;           direction of motion (scalars or arrays of one length)
;   r     - area-equivalent radius sqrt(A B), with semi-axes
;           A = r / sqrt(1 - f) and B = A (1 - f)
;   f     - projected flattening 1 - B/A, in [0, 1)
;   theta - angle of the long axis from X, in RADIANS
;           (r, f and theta: scalars, or one per point)
;   w     - the law's shape weights (2, 4 or 5 of them)
;   law   - 'hybrid2', 'hybrid4' (the default) or 'hybrid5'
;
; OUTPUT:
;   The normalised flux, one value per point (double precision).
;
; RESTRICTIONS:
;   r <= (1 - f)^(3/2), so that the outline crosses the stellar limb at
;   most twice (0.59 at f = 0.3, 0.35 at f = 0.5).
;
; EXAMPLE:
;   X = dindgen(1001) / 500d - 1d
;   F = limbdarkfast_hybrid_flux_oblate(X, 0.3d, 0.1d, 0.2d, 0.52d, $
;                                       [0.3d, 0.2d, 0.1d, 0.05d], 'hybrid4')
;
; CITATION AND LICENSE:
;   If you use this code, please cite
;
;       Kipping, D. 2027, "Limb Darkened Transits of Oblate Planets as
;       Elementary Functions", AJ, submitted.
;
;   Copyright (c) 2026 David Kipping. Released under the MIT License: see
;   the LICENSE file, or https://opensource.org/licenses/MIT.
;-
function limbdarkfast_hybrid_flux_oblate, X, Y, r, f, theta, w, law
  compile_opt idl2
  limbdarkfast_core
  L = limbdarkfast_law(law)
  if n_elements(w) ne L.nw then $
    message, 'this law takes ' + strtrim(L.nw, 2) + ' weights; got ' + strtrim(n_elements(w), 2)
  n = n_elements(X) > n_elements(Y)
  xx = limbdarkfast_spread(X, n, 'X')
  yy = limbdarkfast_spread(Y, n, 'Y')
  rr = limbdarkfast_spread(r, n, 'r')
  ff = limbdarkfast_spread(f, n, 'f')
  th = limbdarkfast_spread(theta, n, 'theta')
  if min(ff) lt 0 || max(ff) ge 1 then message, 'need 0 <= f < 1'
  if min(rr) lt 0 || max(rr - (1d - ff)^1.5d * (1d + 1d-12)) gt 0 then $
    message, 'need 0 <= r <= (1 - f)^(3/2): the planet must cross the stellar limb at most twice'
  if n eq 1 then begin         ; never one element (see limbdarkfast_where in limbdarkfast_core)
    xx = [xx, xx] & yy = [yy, yy] & rr = [rr, rr] & ff = [ff, ff] & th = [th, th]
  endif
  x0 = xx * cos(th) + yy * sin(th)            ; principal frame
  y0 = -xx * sin(th) + yy * cos(th)
  flux = (limbdarkfast_flux_from_occ(limbdarkfast_oblate_occ(x0, y0, rr, ff, L), double(w), L))[0:n - 1]
  void = check_math()          ; branches never taken may raise flags; clear them
  big = n_elements(X) ge n_elements(Y) ? X : Y
  return, size(big, /n_dimensions) eq 0 ? flux[0] : reform(flux, size(big, /dimensions))
end
