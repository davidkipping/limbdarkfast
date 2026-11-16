;+
; LIMBDARKFAST_CORE
;
; The internals of limbdarkfast's IDL version: the flux of a star with
; hybrid2, hybrid4 or hybrid5 limb darkening, occulted by a spherical or an
; oblate planet, in closed form and double precision. Users call
; LIMBDARKFAST_HYBRID_FLUX, LIMBDARKFAST_HYBRID_FLUX_OBLATE and
; LIMBDARKFAST_INTENSITY; each compiles this file on first use (by calling
; the empty procedure LIMBDARKFAST_CORE at its end).
;
; If you use this code, please cite
;
;     Kipping, D. 2027, "Limb Darkened Transits of Oblate Planets as
;     Elementary Functions", AJ, submitted.
;
; Copyright (c) 2026 David Kipping. Released under the MIT License: see
; the LICENSE file, or https://opensource.org/licenses/MIT.
;
; Method (python/limbdarkfast.py has the same structure, with more comments).
; Green's theorem with the radial potential h(s) = (1/2s) int_0^s I,
; s = x^2 + y^2, turns the occulted flux into h(1) times the limb arc inside
; the planet plus the integral of h(s) along the planet's outline inside the
; star. For a circle that reduces to the lens geometry; for an ellipse the
; even terms are finite Fourier sums and each pole a partial-fraction sum
; over the roots of a quartic, near-double roots being integrated as one
; quadratic factor.
;
; Everything is vectorised over points: arrays of length n, and n x 4 for
; the quartic's roots. Each regime and branch is evaluated on the points it
; owns (WHERE subsets), as EXOFAST's occultquad does.
;-

; ---------------------------------------------------------------------------
; constants and small helpers
; ---------------------------------------------------------------------------

; A law: nw weights and K poles eps; G[*, j] gives shape S_j as a combination
; of the generators [mu^0, mu^2, mu^4, (mu^2 + eps_k)^-2 ...]; norms[j] is the
; unocculted flux of S_j. From tools/law_constants.py.
function limbdarkfast_law, law
  compile_opt idl2, hidden
  if n_elements(law) eq 0 then law = 'hybrid4'
  case strlowcase(string(law)) of
    'hybrid2': return, {nw: 2, K: 1, $
        eps: [0.35999999999999999d, 0d, 0d], $
        G: reform([ $
          1d, -1d, 0d, 0d, 0d, 0d, $
          -0.075348837209302341d, 0d, 0d, 0.13936520930232557d, 0d, 0d, $
          0d, 0d, 0d, 0d, 0d, 0d, $
          0d, 0d, 0d, 0d, 0d, 0d, $
          0d, 0d, 0d, 0d, 0d, 0d], 6, 5), $
        norms: [1.5707963267948966d, 0.6575426484257707d, 0d, 0d, 0d]}
    'hybrid4': return, {nw: 4, K: 2, $
        eps: [0.028601207249486824d, 0.38582417584735224d, 0d], $
        G: reform([ $
          1d, -1d, 0d, 0d, 0d, 0d, $
          1d, -2d, 1d, 0d, 0d, 0d, $
          0.0031720957997626959d, -0.0028080275108228465d, 0.0007779295064964302d, $
          0.0008358238567779617d, -0.003710396155608333d, 0d, $
          -0.45048769921539672d, 0.51370715445091153d, -0.17569337890652659d, $
          -4.8355348357818505d-07, 0.21600802064223382d, 0d, $
          0d, 0d, 0d, 0d, 0d, 0d], 6, 5), $
        norms: [1.5707963267948966d, 1.0471975511965974d, 0.073823572890099914d, $
                0.47682011132930302d, 0d]}
    'hybrid5': return, {nw: 5, K: 3, $
        eps: [0.0015965451792942786d, 0.041275044772743019d, 0.38743471584956696d], $
        G: reform([ $
          1d, 0d, -1d, 0d, 0d, 0d, $
          1d, -2d, 1d, 0d, 0d, 0d, $
          -2.5408433670697262d-06, 0d, 0d, 2.5489629860270405d-06, 0d, 0d, $
          -0.0015737187012817282d, 0d, 0d, 0d, 0.0017063103543144378d, 0d, $
          -0.084572789617394117d, 0d, 0d, 0d, 0d, 0.16280051336826454d], 6, 5), $
        norms: [2.0943951023931957d, 1.0471975511965974d, 0.0049997300166832875d, $
                0.11978141121311216d, 0.68577554794744777d]}
    else: message, "unknown law '" + string(law) + "'; expected 'hybrid2', 'hybrid4' or 'hybrid5'"
  endcase
end

; the tolerances of the method, as python/limbdarkfast.py
function limbdarkfast_const
  compile_opt idl2, hidden
  return, {eps: 2.220446049250313d-16, tiny: 1d-280, $
           y_sw: 1d-3, $          ; spherical pole: series half-width in y
           f_sw: 1d-14, $         ; below this flattening, the spherical forms (error ~ 200 f)
           z_tol: 1d-4, $         ; |z| ~ 1: a candidate limb crossing
           g_tol: 1d-12, $        ; accepted crossing residual, relative
           pair_tol: 1d-6, $      ; crossings closer than this are a tangency
           pair_ratio: 0.1d, $    ; near-double pole roots
           n_series: 18, n_aberth: 8, n_bairstow: 2, $
           qr_sweeps4: 6, qr_sweeps3: 5, qr_balance: 4}
end

; WHERE, but never a single index: one match is returned twice. GDL (the
; weekly arm64 build of November 2026) gets arithmetic between one-element
; arrays wrong -- r = a / b leaves the quotient in a and returns b -- so no
; array in this calculation is allowed to have exactly one element. The
; duplicate index writes the same value twice; IDL itself is unaffected.
function limbdarkfast_where, cond, cnt, complement=comp, ncomplement=ncomp
  compile_opt idl2, hidden
  idx = where(cond, cnt, complement=comp, ncomplement=ncomp)
  if cnt eq 1 then idx = [idx, idx]
  if arg_present(comp) && ncomp eq 1 then comp = [comp, comp]
  return, idx
end

; cond ? a : b, elementwise; a and b are arrays of cond's length or scalars.
; Only the selected elements of a are read, so the others may be anything.
function limbdarkfast_choose, cond, a, b
  compile_opt idl2, hidden
  n = n_elements(cond)
  cplx = size(a, /type) eq 9 || size(b, /type) eq 9
  out = cplx ? dcomplexarr(n) : dblarr(n)   ; a fresh array, not b itself
  out[*] = (n_elements(b) eq 1) ? b[0] : b
  idx = limbdarkfast_where(cond, cnt)
  if cnt gt 0 then out[idx] = (n_elements(a) eq 1) ? a[0] : a[idx]
  return, out
end

; t mod 2 pi in [0, 2 pi), in floating point (FLOOR would return an integer)
function limbdarkfast_wrap2pi, t
  compile_opt idl2, hidden
  r = t mod (2d * !dpi)
  return, r + (2d * !dpi) * (r lt 0)
end

function limbdarkfast_arg, z
  compile_opt idl2, hidden
  return, atan(imaginary(z), double(z))
end

function limbdarkfast_expi, t
  compile_opt idl2, hidden
  return, dcomplex(cos(t), sin(t))
end

; an argument broadcast to n points, always as a fresh array, so that nothing
; done to it here can reach the caller's variable
function limbdarkfast_spread, v, n, name
  compile_opt idl2, hidden
  if n_elements(v) ne 1 && n_elements(v) ne n then $
    message, name + ' must be a scalar or have one value per point'
  out = dblarr(n)
  out[*] = (n_elements(v) eq 1) ? double(v[0]) : double(v)
  return, out
end

; F = 1 - (O_0 - sum_j w_j O_Sj) / (pi - sum_j w_j N_j), from the occulted
; flux occ[*, k] of each generator
function limbdarkfast_flux_from_occ, occ, w, L
  compile_opt idl2, hidden
  blocked = occ[*, 0]
  norm = !dpi
  for j = 0, L.nw - 1 do begin
    for k = 0, 2 + L.K do blocked = blocked - w[j] * L.G[k, j] * occ[*, k]
    norm = norm - w[j] * L.norms[j]
  endfor
  return, 1d - blocked / norm
end

; ---------------------------------------------------------------------------
; spherical planet
; ---------------------------------------------------------------------------

; 16 x (area of the triangle with sides 1, z, r)^2, sorted (Kahan) form;
; negative when no such triangle exists
function limbdarkfast_kite_sqarea, z, r
  compile_opt idl2, hidden
  a = z > r
  c = z < r
  b = a < 1d
  a = a > 1d
  bb = b > c
  cc = b < c
  return, (a + (bb + cc)) * (cc - (a - bb)) * (cc + (a - bb)) * (a + (bb - cc))
end

; occulted flux [n, 3 + K] of [mu^0, mu^2, mu^4, P_eps...] for circles of
; radius r at separation z (arrays of n)
function limbdarkfast_sphere_occ, z_in, r_in, L
  compile_opt idl2, hidden
  C = limbdarkfast_const()
  z = abs(double(z_in))
  n = n_elements(z)
  r = limbdarkfast_spread(r_in, n, 'r')
  occ = dblarr(n, 3 + L.K)

  ; wholly on the disc
  ic = limbdarkfast_where(z le 1d - r and z lt 1d + r, nc)
  if nc gt 0 then begin
    zz = z[ic] & rr = r[ic]
    r2 = rr * rr & z2 = zz * zz
    rm1 = (rr - 1d) * (rr + 1d)    ; r^2 - 1, factored: z^2 is lost against it when r = 1, z small
    alpha = -rm1 - z2
    sqarea = limbdarkfast_kite_sqarea(zz, rr)
    A = ((rr + 1d - zz) * (1d - rr + zz)) > C.tiny          ; 1 - (z - r)^2
    ; 1 - (z + r)^2, the larger of z and r taken from 1 first (exactly, when it
    ; is near 1), so that the smaller is not rounded against 1
    Bp = ((1d - (zz > rr)) - (zz < rr)) * (1d + zz + rr)
    s0 = -!dpi * r2
    s2 = 2d * !dpi * r2 * (r2 + 2d * z2 - 1d)
    M0 = !dpi
    M2 = !dpi * alpha
    M4 = (6d * alpha * M2 + 2d * sqarea * M0) / 4d
    s4 = -(2d * r2 * M4 - (2d / 3d) * (alpha * M4 + sqarea * M2))
    occ[ic, 0] = -s0
    occ[ic, 1] = -(0.5d * s0 + 0.25d * s2)
    occ[ic, 2] = -(s0 / 3d + s2 / 6d + s4 / 6d)
    for k = 0, L.K - 1 do begin
      e = L.eps[k]
      p = 1d + e
      sqQ = sqrt((e + A) * (e + Bp))
      occ[ic, 3 + k] = 2d * !dpi * r2 / (sqQ * (p + r2 - z2 + sqQ))
    endfor
  endif

  ; partial overlap (or, for r > 1, total occultation: kite = 0 there, and
  ; every term lands on its full-disc value)
  ip = limbdarkfast_where(z gt 1d - r and z lt 1d + r, np)
  if np gt 0 then begin
    zz = z[ip] & rr = r[ip]
    r2 = rr * rr & z2 = zz * zz
    rm1 = (rr - 1d) * (rr + 1d)    ; r^2 - 1, factored: z^2 is lost against it when r = 1, z small
    alpha = -rm1 - z2
    sqarea = limbdarkfast_kite_sqarea(zz, rr)
    A = ((rr + 1d - zz) * (1d - rr + zz)) > C.tiny
    Bp = ((1d - (zz > rr)) - (zz < rr)) * (1d + zz + rr)
    tot = zz le rr - 1d
    kite = sqrt(sqarea > C.tiny)           ; sqarea > 0 in the lens, to a few ulps
    it = limbdarkfast_where(tot, nt)
    if nt gt 0 then kite[it] = 0d
    kap0 = atan(kite, rm1 + z2)            ; half-angle at the planet
    kap1 = atan(kite, z2 - rm1)            ; half-angle at the star
    s0 = -(kap1 + r2 * kap0 - 0.5d * kite)
    eta2 = r2 * (r2 + 2d * z2)
    s2 = 2d * s0 + 2d * (kap1 + eta2 * kap0 - 0.25d * kite * (1d + 5d * r2 + z2))
    M0 = kap0
    M2 = kap0 * alpha + kite
    M4 = (6d * alpha * M2 + 2d * sqarea * M0) / 4d
    s4 = -(2d * r2 * M4 - (2d / 3d) * (alpha * M4 + sqarea * M2))
    occ[ip, 0] = -s0
    occ[ip, 1] = -(0.5d * s0 + 0.25d * s2)
    occ[ip, 2] = -(s0 / 3d + s2 / 6d + s4 / 6d)

    ; the smaller of A and -Bp is taken from kite^2 = A (-Bp), exactly
    nBp = kite * kite / A
    ib = limbdarkfast_where(~(~tot and A ge -Bp) and -Bp gt 0, nb)
    if nb gt 0 then nBp[ib] = -Bp[ib]
    for k = 0, L.K - 1 do begin
      e = L.eps[k]
      p = 1d + e
      Kp = p + r2 - z2
      Q = (e + A) * (e + Bp)
      ; J = int_{-kap0}^{kap0} dpsi / (p - z^2 - r^2 + 2 z r cos psi)
      apb = e + A
      U = apb * nBp
      ratio = (e + Bp) / nBp
      y = ratio * kite * kite / U
      J = dblarr(np)
      ia = limbdarkfast_where(y gt C.y_sw, na)
      if na gt 0 then $
        J[ia] = 4d / sqrt(Q[ia]) * atan(kite[ia] * sqrt(ratio[ia]), sqrt(U[ia]))
      il = limbdarkfast_where(y lt -C.y_sw, nl)
      if nl gt 0 then $
        J[il] = 4d / sqrt(-Q[il]) * alog((sqrt(U[il]) + kite[il] * sqrt(-ratio[il])) $
                / (2d * sqrt((zz[il] * rr[il] * e) > C.tiny)))
      is = limbdarkfast_where(abs(y) le C.y_sw, ns)       ; across Q = 0: a series
      if ns gt 0 then begin
        ys = y[is]
        y2 = ys * ys
        G = 1d - ys / 3d + y2 / 5d - ys * y2 / 7d + y2 * y2 / 9d - ys * y2 * y2 / 11d
        J[is] = 4d * (kite[is] / nBp[is]) * G / apb[is]     ; kite/nBp = tan(kap0/2)
      endif
      occ[ip, 3 + k] = kap1 / (p * e) + (Kp * J - 2d * kap0) / (4d * p)
    endfor
  endif
  return, occ
end

; ---------------------------------------------------------------------------
; oblate planet: the ellipse x = x0 + a cos phi, y = y0 + b sin phi. With
; z = e^{i phi}, s(phi) = x^2 + y^2 = S0 + 2 Re(S1 z) + 2 S2 cos 2 phi and
; (x dy - y dx)/dphi = K0 + 2 Re(K1 z); the points where s = p are the roots
; of the self-inversive quartic D(z) = S2 z^4 + S1 z^3 + (S0 - p) z^2 +
; conj(S1) z + S2. Every function takes one value per point.
; ---------------------------------------------------------------------------

pro limbdarkfast_coeffs, x0, y0, a, b, S0, S1, S2, K0, K1
  compile_opt idl2, hidden
  S0 = x0 * x0 + y0 * y0 + 0.5d * (a * a + b * b)
  S1 = dcomplex(a * x0, -b * y0)
  S2 = 0.25d * (a * a - b * b)
  K0 = a * b
  K1 = dcomplex(0.5d * b * x0, -0.5d * a * y0)
end

function limbdarkfast_D, z, S0, S1, S2, p
  compile_opt idl2, hidden
  return, (((S2 * z + S1) * z + (S0 - p)) * z + conj(S1)) * z + S2
end

function limbdarkfast_dD, z, S0, S1, S2, p
  compile_opt idl2, hidden
  return, ((4d * S2 * z + 3d * S1) * z + 2d * (S0 - p)) * z + conj(S1)
end

function limbdarkfast_N, z, K0, K1
  compile_opt idl2, hidden
  return, conj(K1) + K0 * z + K1 * z * z
end

; ---- the quartic's roots: a balanced, shifted 4x4 QR, then Aberth --------

; (c, s), c real: [[c, s], [-conj(s), c]] [x, y] = [rho, 0]
pro limbdarkfast_givens, x, y, c, s
  compile_opt idl2, hidden
  n = n_elements(x)
  xa = abs(x)
  rho = sqrt(xa * xa + abs(y)^2)
  ph = replicate(dcomplex(1d, 0d), n)
  i = limbdarkfast_where(xa gt 0, cnt)
  if cnt gt 0 then ph[i] = x[i] / xa[i]
  c = replicate(1d, n)
  s = dcomplexarr(n)
  i = limbdarkfast_where(rho gt 0, cnt)
  if cnt gt 0 then begin
    c[i] = xa[i] / rho[i]
    s[i] = ph[i] * conj(y[i]) / rho[i]
  endif
end

; the eigenvalue of [[a, b], [c, d]] nearer d
function limbdarkfast_wilkinson, a, b, c, d
  compile_opt idl2, hidden
  h = 0.5d * (a - d)
  disc = sqrt(h * h + b * c)
  den = limbdarkfast_choose(abs(h + disc) ge abs(h - disc), h + disc, h - disc)
  out = d
  i = limbdarkfast_where(abs(den) gt 0, cnt)
  if cnt gt 0 then out[i] = d[i] - b[i] * c[i] / den[i]
  return, out
end

; one shifted QR step on the leading m x m block of H [n, 4, 4]
pro limbdarkfast_qr_sweep, H, m
  compile_opt idl2, hidden
  n = (size(H, /dimensions))[0]
  sig = limbdarkfast_wilkinson(H[*, m - 2, m - 2], H[*, m - 2, m - 1], $
                         H[*, m - 1, m - 2], H[*, m - 1, m - 1])
  for i = 0, m - 1 do H[*, i, i] = H[*, i, i] - sig
  cs = dblarr(n, 3)
  sn = dcomplexarr(n, 3)
  for k = 0, m - 2 do begin
    limbdarkfast_givens, H[*, k, k], H[*, k + 1, k], c, s
    cs[*, k] = c
    sn[*, k] = s
    for j = k, m - 1 do begin
      x = H[*, k, j]
      y = H[*, k + 1, j]
      H[*, k, j] = c * x + s * y
      H[*, k + 1, j] = c * y - conj(s) * x
    endfor
  endfor
  for k = 0, m - 2 do begin
    c = cs[*, k]
    s = sn[*, k]
    for i = 0, (k + 2) < (m - 1) do begin
      x = H[*, i, k]
      y = H[*, i, k + 1]
      H[*, i, k] = c * x + conj(s) * y
      H[*, i, k + 1] = c * y - s * x
    endfor
  endfor
  for i = 0, m - 1 do H[*, i, i] = H[*, i, i] + sig
end

; eigenvalues [n, 4] of D's companion matrix
function limbdarkfast_qr_seeds, S0, S1, S2, p
  compile_opt idl2, hidden
  C = limbdarkfast_const()
  n = n_elements(S0)
  H = dcomplexarr(n, 4, 4)
  inv = 1d / S2
  H[*, 0, 0] = -S1 * inv
  H[*, 0, 1] = -(S0 - p) * inv
  H[*, 0, 2] = -conj(S1) * inv
  H[*, 0, 3] = dcomplex(-1d, 0d)
  H[*, 1, 0] = dcomplex(1d, 0d)
  H[*, 2, 1] = dcomplex(1d, 0d)
  H[*, 3, 2] = dcomplex(1d, 0d)
  for it = 0, C.qr_balance - 1 do begin          ; diagonal similarity
    for i = 0, 3 do begin
      col = dblarr(n)
      row = dblarr(n)
      for j = 0, 3 do begin
        if j eq i then continue
        col = col + abs(H[*, j, i])
        row = row + abs(H[*, i, j])
      endfor
      f = replicate(1d, n)
      ok = limbdarkfast_where(col gt 0 and row gt 0, cnt)
      if cnt gt 0 then f[ok] = sqrt(row[ok] / col[ok])
      for j = 0, 3 do begin
        H[*, j, i] = H[*, j, i] * f
        H[*, i, j] = H[*, i, j] / f
      endfor
    endfor
  endfor
  for it = 0, C.qr_sweeps4 - 1 do limbdarkfast_qr_sweep, H, 4
  for it = 0, C.qr_sweeps3 - 1 do limbdarkfast_qr_sweep, H, 3
  aa = H[*, 0, 0] & bb = H[*, 0, 1] & cc = H[*, 1, 0] & dd = H[*, 1, 1]
  hh = 0.5d * (aa - dd)
  disc = sqrt(hh * hh + bb * cc)
  mid = 0.5d * (aa + dd)
  z = dcomplexarr(n, 4)
  z[*, 0] = mid + disc
  z[*, 1] = mid - disc
  z[*, 2] = H[*, 2, 2]
  z[*, 3] = H[*, 3, 3]
  return, z
end

; the four roots [n, 4] of D at level p: eigenvalue seeds, then Aberth
; iterations (at most n_aberth, stopping once every correction is at
; rounding) and two Newton polishes; a correction that is not finite is
; skipped
function limbdarkfast_roots, S0, S1, S2, p
  compile_opt idl2, hidden
  C = limbdarkfast_const()
  n = n_elements(S0)
  z = limbdarkfast_qr_seeds(S0, S1, S2, p)
  for it = 0, C.n_aberth - 1 do begin
    rat = dcomplexarr(n, 4)
    for i = 0, 3 do rat[*, i] = limbdarkfast_D(z[*, i], S0, S1, S2, p) / limbdarkfast_dD(z[*, i], S0, S1, S2, p)
    znew = z
    busy = 0L
    for i = 0, 3 do begin
      s = dcomplexarr(n)
      for j = 0, 3 do begin
        if j eq i then continue
        diff = z[*, i] - z[*, j]
        ok = limbdarkfast_where(abs(diff) gt 0, cnt)
        if cnt gt 0 then s[ok] = s[ok] + 1d / diff[ok]
      endfor
      step = rat[*, i] / (1d - rat[*, i] * s)
      ok = limbdarkfast_where(finite(step), cnt)
      if cnt gt 0 then begin
        zi = z[*, i]
        zn = zi
        zn[ok] = zi[ok] - step[ok]
        znew[*, i] = zn
        busy = busy + long(total(abs(step[ok]) gt 4d * C.eps * abs(zi[ok])))
      endif
    endfor
    z = znew
    if busy eq 0 then break
  endfor
  for it = 0, 1 do begin
    for i = 0, 3 do begin
      zi = z[*, i]
      nw = zi - limbdarkfast_D(zi, S0, S1, S2, p) / limbdarkfast_dD(zi, S0, S1, S2, p)
      ok = limbdarkfast_where(finite(nw), cnt)
      if cnt gt 0 then zi[ok] = nw[ok]
      z[*, i] = zi
    endfor
  endfor
  return, z
end

; ---- limb crossings and the arcs ------------------------------------------

pro limbdarkfast_g_dg, phi, x0, y0, a, b, g, dg
  compile_opt idl2, hidden
  s = sin(phi)
  c = cos(phi)
  x = x0 + a * c
  y = y0 + b * s
  g = x * x + y * y - 1d
  dg = 2d * (-x * a * s + y * b * c)
end

; number of limb crossings, and the first two (phi_a <= phi_b) on the
; ellipse parameter; crossings closer than pair_tol are a tangency and are
; dropped together
function limbdarkfast_crossings, x0, y0, a, b, phi_a, phi_b
  compile_opt idl2, hidden
  C = limbdarkfast_const()
  n = n_elements(x0)
  limbdarkfast_coeffs, x0, y0, a, b, S0, S1, S2, K0, K1
  z = limbdarkfast_roots(S0, S1, S2, 1d)
  sc = 1d + (abs(x0) + a)^2 + (abs(y0) + b)^2
  phi = dblarr(n, 4)
  valid = bytarr(n, 4)
  for i = 0, 3 do begin
    zi = z[*, i]
    cand = abs(abs(zi) - 1d) lt C.z_tol
    ph = limbdarkfast_choose(cand, limbdarkfast_arg(zi), 0d)
    limbdarkfast_g_dg, ph, x0, y0, a, b, g, dg
    for it = 0, 2 do begin
      trial = ph
      ok = limbdarkfast_where(abs(dg) gt 0, cnt)
      if cnt gt 0 then trial[ok] = ph[ok] - g[ok] / dg[ok]
      limbdarkfast_g_dg, trial, x0, y0, a, b, gtr, dgtr
      ok = limbdarkfast_where(abs(gtr) lt abs(g), cnt)
      if cnt gt 0 then begin
        ph[ok] = trial[ok]
        g[ok] = gtr[ok]
        dg[ok] = dgtr[ok]
      endif
    endfor
    valid[*, i] = cand and (abs(g) le C.g_tol * sc)
    phi[*, i] = limbdarkfast_wrap2pi(ph)
  endfor
  drop = bytarr(n, 4)
  for i = 0, 3 do begin
    for j = i + 1, 3 do begin
      d = abs(limbdarkfast_wrap2pi(phi[*, i] - phi[*, j] + !dpi) - !dpi)
      close = valid[*, i] and valid[*, j] and (d lt C.pair_tol)
      drop[*, i] = drop[*, i] or close
      drop[*, j] = drop[*, j] or close
    endfor
  endfor
  valid = valid and ~drop
  ; the two smallest valid angles per point. The sentinel is finite: GDL's
  ; MIN(..., DIMENSION=2) returns subscript 0 for a row that is all Infinity
  big = 1d300
  key = phi
  bad = limbdarkfast_where(~valid, nbad)
  if nbad gt 0 then key[bad] = big
  phi_a = min(key, imin, dimension=2)
  key[imin] = big
  phi_b = min(key, dimension=2)
  return, long(total(valid, 2))
end

; the ellipse arc [lo, hi] inside the star and the length dth of the limb arc
; inside the ellipse, for two crossings. The ellipse arc is decided at the
; midpoint of the longer of its two arcs; the limb arc inside the ellipse is
; always the shorter one (the planet lies within radius A < 1).
pro limbdarkfast_arcs, x0, y0, a, b, phi_a, phi_b, lo, hi, dth
  compile_opt idl2, hidden
  long_first = (phi_b - phi_a) ge !dpi
  mid = 0.5d * (phi_a + phi_b) + !dpi * (~long_first)
  xm = x0 + a * cos(mid)
  ym = y0 + b * sin(mid)
  mid_in = xm * xm + ym * ym lt 1d
  first = limbdarkfast_choose(long_first, mid_in, ~mid_in)   ; is [phi_a, phi_b] inside?
  lo = limbdarkfast_choose(first, phi_a, phi_b)
  hi = limbdarkfast_choose(first, phi_b, phi_a + 2d * !dpi)
  ta = atan(y0 + b * sin(phi_a), x0 + a * cos(phi_a))
  tb = atan(y0 + b * sin(phi_b), x0 + a * cos(phi_b))
  dt = abs(limbdarkfast_wrap2pi(ta) - limbdarkfast_wrap2pi(tb))
  dth = dt < (2d * !dpi - dt)
end

; ---- even terms: Fourier coefficients of h_n(s) K ---------------------------

; q[*, n, m], n = 0..2, m = 0..5: h_n(s) K = Re sum_m q_m e^{i m phi}, with
; h_n(s) = sum_p alpha[n, p] s^p, alpha[n, p] = (-1)^p C(n, p) / (2 (p + 1))
function limbdarkfast_q_coeffs, S0, S1, S2, K0, K1
  compile_opt idl2, hidden
  alpha = [[0.5d, 0.5d, 0.5d], [0d, -0.25d, -0.5d], [0d, 0d, 1d / 6d]]   ; alpha[n, p]
  n = n_elements(S0)
  T = dcomplexarr(n, 11)                 ; harmonics -5..5 at index m + 5
  T[*, 4] = conj(K1)
  T[*, 5] = dcomplex(K0, 0d * K0)
  T[*, 6] = K1
  Sk = dcomplexarr(n, 5)                 ; j = -2..2 at index j + 2
  Sk[*, 0] = dcomplex(S2, 0d * S2)
  Sk[*, 1] = conj(S1)
  Sk[*, 2] = dcomplex(S0, 0d * S0)
  Sk[*, 3] = S1
  Sk[*, 4] = dcomplex(S2, 0d * S2)
  q = dcomplexarr(n, 3, 6)
  for p = 0, 2 do begin
    if p gt 0 then begin                 ; T <- S * T
      Tn = dcomplexarr(n, 11)
      for i = 0, 10 do begin
        for j = -2, 2 do begin
          if i - j ge 0 && i - j le 10 then Tn[*, i] = Tn[*, i] + Sk[*, j + 2] * T[*, i - j]
        endfor
      endfor
      T = Tn
    endif
    for nn = p, 2 do begin
      for m = 0, 5 do q[*, nn, m] = q[*, nn, m] + alpha[nn, p] * T[*, m + 5]
    endfor
  endfor
  return, q
end

; Psi_n(hi) - Psi_n(lo), Psi_n the exact antiderivative, in difference form
function limbdarkfast_even_arc, q, nn, lo, hi
  compile_opt idl2, hidden
  mid = 0.5d * (lo + hi)
  half = 0.5d * (hi - lo)
  acc = double(q[*, nn, 0]) * (hi - lo)
  for m = 1, 5 do $
    acc = acc + (4d / m) * double(q[*, nn, m] * limbdarkfast_expi(m * mid)) * sin(m * half)
  return, acc
end

; ---- poles: partial fractions over the roots of D ---------------------------

; Newton steps on (sigma, pi) for the quadratic factor z^2 - sigma z + pi
pro limbdarkfast_bairstow, S0, S1, S2, p, sig, pq
  compile_opt idl2, hidden
  C = limbdarkfast_const()
  n = n_elements(S0)
  dc = dcomplexarr(n, 5)                 ; low -> high
  dc[*, 0] = dcomplex(S2, 0d * S2)
  dc[*, 1] = conj(S1)
  dc[*, 2] = dcomplex(S0 - p, 0d * S0)
  dc[*, 3] = S1
  dc[*, 4] = dcomplex(S2, 0d * S2)
  for it = 0, C.n_bairstow - 1 do begin
    zero = dcomplexarr(n)
    al = zero                            ; z^k mod q = al z + be
    be = zero + 1d
    das = zero & dbs = zero & dap = zero & dbp = zero
    r1 = dc[*, 0] * al
    r0 = dc[*, 0] * be
    J11 = zero & J12 = zero & J21 = zero & J22 = zero
    for k = 1, 4 do begin
      nal = al * sig + be
      nbe = -al * pq
      ndas = das * sig + al + dbs
      ndbs = -das * pq
      ndap = dap * sig + dbp
      ndbp = -dap * pq - al
      al = nal & be = nbe & das = ndas & dbs = ndbs & dap = ndap & dbp = ndbp
      r1 = r1 + dc[*, k] * al
      r0 = r0 + dc[*, k] * be
      J11 = J11 + dc[*, k] * das
      J12 = J12 + dc[*, k] * dap
      J21 = J21 + dc[*, k] * dbs
      J22 = J22 + dc[*, k] * dbp
    endfor
    det = J11 * J22 - J12 * J21
    ok = limbdarkfast_where(abs(det) gt 0, cnt)
    if cnt gt 0 then begin
      sig[ok] = sig[ok] - (J22[ok] * r1[ok] - J12[ok] * r0[ok]) / det[ok]
      pq[ok] = pq[ok] - (-J21[ok] * r1[ok] + J11[ok] * r0[ok]) / det[ok]
    endif
  endfor
end

; for the quadratic factor q of D, with N = aN z + bN and Q = D/q = aQ z + bQ
; (mod q): the sum of N/D' over its two roots (A_sum) and the C of
; int (A z + B)/q dz (C_num, before the factor i); symmetric in the roots
pro limbdarkfast_pair_terms, S0, S1, S2, K0, K1, p, sig, pq, A_sum, C_num
  compile_opt idl2, hidden
  Q2 = dcomplex(S2, 0d * S2)
  Q1 = S1 + sig * Q2
  Q0 = sig * Q1 + (S0 - p) - pq * Q2
  aN = K1 * sig + K0
  bN = conj(K1) - K1 * pq
  aQ = Q2 * sig + Q1
  bQ = Q0 - Q2 * pq
  den = aQ * aQ * pq + aQ * bQ * sig + bQ * bQ
  bad = limbdarkfast_where(~(abs(den) gt 0), nbad)
  if nbad gt 0 then den[bad] = dcomplex(1d, 0d)
  A_sum = (aN * bQ - bN * aQ) / den
  C_num = 0.5d * (2d * aN * aQ * pq + (aN * bQ + bN * aQ) * sig + 2d * bN * bQ) / den
end

; a whole ellipse on the disc: -(pi/p) Re sum_{|z_r| < 1} N/D'(z_r), the two
; inner roots summed as one quadratic factor
function limbdarkfast_pole_inside, S0, S1, S2, K0, K1, p, z
  compile_opt idl2, hidden
  n = n_elements(S0)
  sig = dcomplexarr(n)
  pq = replicate(dcomplex(1d, 0d), n)
  for i = 0, 3 do begin
    inner = abs(z[*, i]) lt 1d
    sig = sig + limbdarkfast_choose(inner, z[*, i], dcomplex(0d, 0d))
    pq = pq * limbdarkfast_choose(inner, z[*, i], dcomplex(1d, 0d))
  endfor
  limbdarkfast_bairstow, S0, S1, S2, p, sig, pq
  limbdarkfast_pair_terms, S0, S1, S2, K0, K1, p, sig, pq, A_sum, C_num
  return, -!dpi * double(A_sum) / p
end

; change of log(e^{i phi} - z) for phi from lo to hi, continuous on the arc
function limbdarkfast_dlog, z, lo, hi, inner, w1, w2
  compile_opt idl2, hidden
  dm = alog(abs(w2 - z)) - alog(abs(w1 - z))
  zi = limbdarkfast_choose(inner, z, dcomplex(0.5d, 0d))
  zo = limbdarkfast_choose(inner, dcomplex(2d, 0d), z)
  d_in = (hi - lo) + limbdarkfast_arg(1d - zi / w2) - limbdarkfast_arg(1d - zi / w1)
  d_out = limbdarkfast_arg(1d - w2 / zo) - limbdarkfast_arg(1d - w1 / zo)
  return, dcomplex(dm, limbdarkfast_choose(inner, d_in, d_out))
end

; int du / (u^2 - Delta) = -(1/u) sum_k (Delta/u^2)^k / (2k + 1)
function limbdarkfast_jser, u, Delta
  compile_opt idl2, hidden
  C = limbdarkfast_const()
  x = Delta / (u * u)
  acc = dcomplexarr(n_elements(u))
  for k = C.n_series - 1, 0, -1 do acc = acc * x + 1d / (2d * k + 1d)
  return, -acc / u
end

; (1/2p) int_lo^hi K / (p - s) dphi along one ellipse arc inside the star:
; (1/2p) Re sum_r R_r [log(e^{i phi} - z_r)], R_r = -N(z_r) / (i D'(z_r))
function limbdarkfast_pole_arc, S0, S1, S2, K0, K1, p, z, lo, hi
  compile_opt idl2, hidden
  C = limbdarkfast_const()
  n = n_elements(S0)
  ii = dcomplex(0d, 1d)
  w1 = limbdarkfast_expi(lo)
  w2 = limbdarkfast_expi(hi)
  PJ = [0, 0, 0, 1, 1, 2]
  PL = [1, 2, 3, 2, 3, 3]
  ; near-double pairs: up to two disjoint, scored by |Delta| / d^2, d the
  ; distance from the pair's centre to the arc
  used = bytarr(n, 4)
  pj_r = lonarr(n, 2)
  pl_r = lonarr(n, 2)
  good_r = bytarr(n, 2)
  for rnd = 0, 1 do begin
    score = dblarr(n, 6)
    for m = 0, 5 do begin
      zj = z[*, PJ[m]]
      zl = z[*, PL[m]]
      Dl = 0.25d * (zj - zl)^2
      cen = 0.5d * (zj + zl)
      ang = limbdarkfast_wrap2pi(limbdarkfast_arg(cen))
      on = (ang ge lo and ang le hi) or (ang + 2d * !dpi ge lo and ang + 2d * !dpi le hi)
      d = limbdarkfast_choose(on, abs(abs(cen) - 1d), abs(w1 - cen) < abs(w2 - cen))
      sc = abs(Dl) / ((d * d) > 1d-300)
      free = ~(used[*, PJ[m]] or used[*, PL[m]])
      score[*, m] = limbdarkfast_choose(free, sc, 1d300)    ; finite: see limbdarkfast_crossings
    endfor
    best = min(score, ibest, dimension=2)
    mbest = ibest / n                    ; the pair index (column of score)
    good = best le C.pair_ratio
    pj_r[*, rnd] = PJ[mbest]
    pl_r[*, rnd] = PL[mbest]
    good_r[*, rnd] = good
    g = limbdarkfast_where(good, ng)
    if ng gt 0 then begin
      used[g + n * PJ[mbest[g]]] = 1b
      used[g + n * PL[mbest[g]]] = 1b
    endif
  endfor

  ell = dcomplexarr(n)
  ; isolated roots; one beyond |z| = 2 is taken through its partner
  ; zeta = 1/conj(z) inside the disc, where R(z) = conj(R(zeta)) and the
  ; log's change is that of log(1 - e^{i phi} conj(zeta))
  for k = 0, 3 do begin
    zk = z[*, k]
    a2 = abs(zk)^2
    far = a2 gt 4d
    zeta = limbdarkfast_choose(far, zk / a2, zk)
    nw = zeta - limbdarkfast_D(zeta, S0, S1, S2, p) / limbdarkfast_dD(zeta, S0, S1, S2, p)
    better = finite(nw) and (abs(limbdarkfast_D(nw, S0, S1, S2, p)) lt abs(limbdarkfast_D(zeta, S0, S1, S2, p)))
    zeta = limbdarkfast_choose(better, nw, zeta)
    dDv = limbdarkfast_dD(zeta, S0, S1, S2, p)
    live = ~used[*, k] and (abs(dDv) gt 0)
    R = -limbdarkfast_N(zeta, K0, K1) / (ii * limbdarkfast_choose(live, dDv, dcomplex(1d, 0d)))
    b2 = 1d - w2 * conj(zeta)
    b1 = 1d - w1 * conj(zeta)
    dlw = dcomplex(alog(abs(b2)) - alog(abs(b1)), limbdarkfast_arg(b2) - limbdarkfast_arg(b1))
    dln = limbdarkfast_dlog(zeta, lo, hi, abs(zeta) lt 1d, w1, w2)
    ell = ell + limbdarkfast_choose(live and ~far, R * dln, dcomplex(0d, 0d)) $
              + limbdarkfast_choose(live and far, conj(R) * dlw, dcomplex(0d, 0d))
  endfor
  ; near-double pairs, each as one quadratic factor q = z^2 - sigma z + pi:
  ; int (A z + B)/q dz = (A/2) log q + C J(z - sigma/2), the winding of log q
  ; taken from the two roots' own continued logs
  for rnd = 0, 1 do begin
    g = limbdarkfast_where(good_r[*, rnd], ng)
    if ng eq 0 then continue
    z1 = z[g + n * pj_r[g, rnd]]
    z2 = z[g + n * pl_r[g, rnd]]
    gS0 = S0[g] & gS1 = S1[g] & gS2 = S2[g] & gK0 = K0[g] & gK1 = K1[g]
    sig = z1 + z2
    pq = z1 * z2
    limbdarkfast_bairstow, gS0, gS1, gS2, p, sig, pq
    limbdarkfast_pair_terms, gS0, gS1, gS2, gK0, gK1, p, sig, pq, A_sum, C_num
    lo_ = lo[g] & hi_ = hi[g] & a1 = w1[g] & a2_ = w2[g]
    qh = a2_ * a2_ - sig * a2_ + pq
    ql = a1 * a1 - sig * a1 + pq
    darg = limbdarkfast_arg(qh / ql)
    rough = imaginary(limbdarkfast_dlog(z1, lo_, hi_, abs(z1) lt 1d, a1, a2_)) $
          + imaginary(limbdarkfast_dlog(z2, lo_, hi_, abs(z2) lt 1d, a1, a2_))
    wind = double(round((rough - darg) / (2d * !dpi)))
    LQ = dcomplex(alog(abs(qh)) - alog(abs(ql)), darg + 2d * !dpi * wind)
    Delta = 0.25d * sig * sig - pq
    cen = 0.5d * sig
    Jd = limbdarkfast_jser(a2_ - cen, Delta) - limbdarkfast_jser(a1 - cen, Delta)
    ell[g] = ell[g] + 0.5d * (ii * A_sum) * LQ + (ii * C_num) * Jd
  endfor
  return, double(ell) / (2d * p)
end

; occulted flux [n, 3 + K] of [mu^0, mu^2, mu^4, P_eps...] for ellipses with
; centre (x0, y0) in their principal frame and semi-axes a > b
function limbdarkfast_ellipse_occ, x0, y0, a, b, L
  compile_opt idl2, hidden
  n = n_elements(x0)
  occ = dblarr(n, 3 + L.K)
  d = sqrt(x0 * x0 + y0 * y0)
  fast_in = d + a lt 1d
  general = ~fast_in and (d - a lt 1d)
  cnt = lonarr(n)
  phi_a = dblarr(n)
  phi_b = dblarr(n)
  gi = limbdarkfast_where(general, ng)
  if ng gt 0 then begin
    cnt[gi] = limbdarkfast_crossings(x0[gi], y0[gi], a[gi], b[gi], pa, pb)
    phi_a[gi] = pa
    phi_b[gi] = pb
  endif
  partial = general and (cnt eq 2)
  inside = fast_in or (general and (cnt ne 2) and (d lt 1d))

  ii = limbdarkfast_where(inside, ni)
  if ni gt 0 then begin
    limbdarkfast_coeffs, x0[ii], y0[ii], a[ii], b[ii], S0, S1, S2, K0, K1
    q = limbdarkfast_q_coeffs(S0, S1, S2, K0, K1)
    for nn = 0, 2 do occ[ii, nn] = 2d * !dpi * double(q[*, nn, 0])
    for k = 0, L.K - 1 do begin
      p = 1d + L.eps[k]
      z = limbdarkfast_roots(S0, S1, S2, p)
      occ[ii, 3 + k] = limbdarkfast_pole_inside(S0, S1, S2, K0, K1, p, z)
    endfor
  endif

  ip = limbdarkfast_where(partial, np)
  if np gt 0 then begin
    px0 = x0[ip] & py0 = y0[ip] & sa = a[ip] & sb = b[ip]
    limbdarkfast_arcs, px0, py0, sa, sb, phi_a[ip], phi_b[ip], lo, hi, dth
    limbdarkfast_coeffs, px0, py0, sa, sb, S0, S1, S2, K0, K1
    q = limbdarkfast_q_coeffs(S0, S1, S2, K0, K1)
    for nn = 0, 2 do occ[ip, nn] = limbdarkfast_even_arc(q, nn, lo, hi) + dth / (2d * (nn + 1))
    for k = 0, L.K - 1 do begin
      e = L.eps[k]
      p = 1d + e
      z = limbdarkfast_roots(S0, S1, S2, p)
      occ[ip, 3 + k] = limbdarkfast_pole_arc(S0, S1, S2, K0, K1, p, z, lo, hi) + dth / (2d * p * e)
    endfor
  endif
  return, occ
end

; occulted generator fluxes [n, 3 + K] for oblate planets (arrays of n)
function limbdarkfast_oblate_occ, x0, y0, r, f, L
  compile_opt idl2, hidden
  C = limbdarkfast_const()
  n = n_elements(x0)
  occ = dblarr(n, 3 + L.K)
  is = limbdarkfast_where(f lt C.f_sw, ns, complement=ie, ncomplement=nell)
  if ns gt 0 then occ[is, *] = limbdarkfast_sphere_occ(sqrt(x0[is]^2 + y0[is]^2), r[is], L)
  if nell gt 0 then begin
    a = r[ie] / sqrt(1d - f[ie])
    occ[ie, *] = limbdarkfast_ellipse_occ(x0[ie], y0[ie], a, a * (1d - f[ie]), L)
  endif
  return, occ
end

; called by the public routines so that IDL compiles this file
pro limbdarkfast_core
  compile_opt idl2, hidden
end
