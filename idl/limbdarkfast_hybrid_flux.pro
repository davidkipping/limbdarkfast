;+
; NAME:
;   LIMBDARKFAST_HYBRID_FLUX
;
; PURPOSE:
;   Normalised flux of a star with hybrid2, hybrid4 or hybrid5 limb
;   darkening, occulted by a spherical planet, in closed form (the IDL
;   version of limbdarkfast).
;
; CALLING SEQUENCE:
;   F = LIMBDARKFAST_HYBRID_FLUX(z, r, w [, law])
;
; INPUTS:
;   z   - sky separation of the planet's centre from the star's, in stellar
;         radii (scalar or array)
;   r   - planet radius in stellar radii (scalar, or one per z)
;   w   - the law's shape weights (2, 4 or 5 of them):
;         I(mu)/I(1) = 1 - sum_j w_j S_j(mu), each shape falling from 1 at
;         the limb to 0 at the centre
;   law - 'hybrid2', 'hybrid4' (the default) or 'hybrid5'
;
; OUTPUT:
;   The normalised flux, one value per z (double precision). Exact for
;   r < 1; for a companion larger than the star the even terms lose about
;   r^5 x 1e-16.
;
; EXAMPLE:
;   z = dindgen(1001) / 1000d * 1.2d
;   F = limbdarkfast_hybrid_flux(z, 0.1d, [0.3d, 0.2d, 0.1d, 0.05d], 'hybrid4')
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
function limbdarkfast_hybrid_flux, z, r, w, law
  compile_opt idl2
  limbdarkfast_core
  L = limbdarkfast_law(law)
  if n_elements(w) ne L.nw then $
    message, 'this law takes ' + strtrim(L.nw, 2) + ' weights; got ' + strtrim(n_elements(w), 2)
  if min(r) lt 0 then message, 'need r >= 0'
  n = n_elements(z)
  zz = limbdarkfast_spread(z, n, 'z')
  rr = limbdarkfast_spread(r, n, 'r')
  if n eq 1 then begin
    zz = [zz, zz]
    rr = [rr, rr]
  endif
  occ = limbdarkfast_sphere_occ(zz, rr, L)
  flux = (limbdarkfast_flux_from_occ(occ, double(w), L))[0:n - 1]
  void = check_math()          ; branches never taken may raise flags; clear them
  return, size(z, /n_dimensions) eq 0 ? flux[0] : reform(flux, size(z, /dimensions))
end
