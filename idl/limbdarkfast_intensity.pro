;+
; NAME:
;   LIMBDARKFAST_INTENSITY
;
; PURPOSE:
;   The hybrid law's intensity profile, I(mu)/I(1) = 1 - sum_j w_j S_j(mu).
;
; CALLING SEQUENCE:
;   I = LIMBDARKFAST_INTENSITY(mu, w [, law])
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
function limbdarkfast_intensity, mu, w, law
  compile_opt idl2
  limbdarkfast_core
  L = limbdarkfast_law(law)
  if n_elements(w) ne L.nw then $
    message, 'this law takes ' + strtrim(L.nw, 2) + ' weights; got ' + strtrim(n_elements(w), 2)
  n = n_elements(mu)
  m2 = limbdarkfast_spread(mu, n, 'mu')^2
  if n eq 1 then m2 = [m2, m2]   ; never one element (see limbdarkfast_where in limbdarkfast_core)
  I = 1d + 0d * m2
  for j = 0, L.nw - 1 do begin
    S = L.G[0, j] + L.G[1, j] * m2 + L.G[2, j] * m2 * m2
    for k = 0, L.K - 1 do S = S + L.G[3 + k, j] / (m2 + L.eps[k])^2
    I = I - w[j] * S
  endfor
  I = I[0:n - 1]
  return, size(mu, /n_dimensions) eq 0 ? I[0] : reform(I, size(mu, /dimensions))
end
