;+
; TEST_LIMBDARKFAST
;
; Checks the IDL version against the shared test vectors. From this directory,
; in IDL or GDL:
;
;   IDL> test_limbdarkfast                      ; vectors in ../tests/vectors
;   IDL> test_limbdarkfast, dir='path/to/tests/vectors'
;
; Prints the largest difference per file, and stops with an error if any
; exceeds 1e-12.
;-

; read a vector file: rows [ncol, nrow] after the '#' header lines
function test_limbdarkfast_read, path, ncol
  compile_opt idl2
  nlines = file_lines(path)
  openr, lun, path, /get_lun
  line = ''
  nhead = 0L
  while ~eof(lun) do begin
    point = (fstat(lun)).cur_ptr
    readf, lun, line
    if strmid(line, 0, 1) ne '#' then begin
      point_lun, lun, point
      break
    endif
    nhead = nhead + 1
  endwhile
  data = dblarr(ncol, nlines - nhead)
  readf, lun, data
  free_lun, lun
  return, data
end

; ---------------------------------------------------------------------------
; brute-force integration, as in Agol, Luger & Foreman-Mackey (2020)
; (test/test_transit_poly.jl of Limbdark.jl): sum the intensity over NS
; annuli of the star, each weighted by the angle Phi(s) it spends inside
; the planet, with the midpoint rule (NS = 5000; 48 samples in angle before
; bisecting for crossings)
; ---------------------------------------------------------------------------

function test_limbdarkfast_total, w, law
  compile_opt idl2
  s = (dindgen(5000) + 0.5d) / 5000d
  return, total(2d * !dpi * s / 5000d * limbdarkfast_intensity(sqrt(1d - s * s), w, law))
end

; spherical planet, by brute force
function test_limbdarkfast_numeric, z, r, w, law
  compile_opt idl2
  if z ge 1d + r then return, 1d
  s1 = 0d > (z - r) & s2 = 1d < (z + r)
  ds = (s2 - s1) / 5000d
  s = s1 + (dindgen(5000) + 0.5d) * ds
  phi = 2d * acos(((s * s + z * z - r * r) / (2d * z * s + 1d-300) > (-1d)) < 1d)
  full = where(s lt r - z, nfull)
  if nfull gt 0 then phi[full] = 2d * !dpi
  occ = total(s * ds * phi * limbdarkfast_intensity(sqrt(1d - s * s), w, law))
  return, 1d - occ / test_limbdarkfast_total(w, law)
end

; The angle of each circle of radius s (array) inside the ellipse
; ((x - x0)/A)^2 + ((y - y0)/B)^2 <= 1: sample 49 angles across the window
; the ellipse can occupy, bisect each sign change, and add up the parts of
; the intervals inside.
function test_limbdarkfast_arc, s, x0, y0, A, B
  compile_opt idl2
  m_arc = 48
  d = sqrt(x0 * x0 + y0 * y0)
  lo = 0d & hi = 2d * !dpi
  if d gt A then begin
    half = asin(A / d)
    lo = atan(y0, x0) - half & hi = atan(y0, x0) + half
  endif
  ns = n_elements(s)
  p = lo + (hi - lo) * dindgen(m_arc + 1) / m_arc
  inside = ((s # cos(p) - x0) / A)^2 + ((s # sin(p) - y0) / B)^2 lt 1d
  ina = inside[*, 0:m_arc - 1]
  part = (ina and inside[*, 1:m_arc]) * ((hi - lo) / m_arc)
  idx = where(ina ne inside[*, 1:m_arc], nc)
  if nc gt 0 then begin
    if nc eq 1 then idx = [idx, idx]          ; never one element (GDL)
    is = idx mod ns & k = idx / ns
    si = s[is] & pa = p[k] & pb = p[k + 1] & sa = ina[idx]
    a_ = pa & b_ = pb
    for it = 0, 49 do begin
      m = 0.5d * (a_ + b_)
      same = (((si * cos(m) - x0) / A)^2 + ((si * sin(m) - y0) / B)^2 lt 1d) eq sa
      a_ = same * m + (1 - same) * a_
      b_ = same * b_ + (1 - same) * m
    endfor
    cross = 0.5d * (a_ + b_)
    part[idx] = sa * (cross - pa) + (1 - sa) * (pb - cross)
  endif
  return, total(part, 2)
end

; oblate planet, by brute force
function test_limbdarkfast_numeric_oblate, X, Y, r, f, theta, w, law
  compile_opt idl2
  A = r / sqrt(1d - f) & B = A * (1d - f)
  x0 = X * cos(theta) + Y * sin(theta)
  y0 = -X * sin(theta) + Y * cos(theta)
  d = sqrt(x0 * x0 + y0 * y0)
  s1 = 0d > (d - A) & s2 = 1d < (d + A)
  if s1 ge s2 then return, 1d
  ds = (s2 - s1) / 5000d
  s = s1 + (dindgen(5000) + 0.5d) * ds
  occ = total(s * ds * test_limbdarkfast_arc(s, x0, y0, A, B) $
              * limbdarkfast_intensity(sqrt(1d - s * s), w, law))
  return, 1d - occ / test_limbdarkfast_total(w, law)
end

; Largest difference from brute force: their table of (r, b) with every
; contact (r <= 2), its oblate twin (f = 0.2, the planet displaced along
; 0.6 rad and tilted by 0.4 rad), then a light curve for r = 0.1 (1000
; points spherical, 201 oblate).
function test_limbdarkfast_numeric_checks
  compile_opt idl2
  laws = ['hybrid2', 'hybrid4', 'hybrid5']
  wlist = list([0.12123702682335066d, 0.22967167533748584d], $
           [0.3709047608728925d, 0.00711438176464726d, 0.3162293789763985d, 0.01888876666845794d], $
           [0.136725258105043d, 0.25114827117024674d, 0.1088217352833381d, $
            0.030866335871229982d, 0.32786772243688317d])
  rs = [0.01d, 0.01d, 0.01d, 0.01d, 0.01d, 0.1d, 0.1d, 0.1d, 0.1d, 0.1d, $
        1d, 1d, 1d, 2d, 2d, 2d, 0.5d, 0.5d, 0.5d, 0.5d]
  bs = [0d, 0.01d, 0.99d, 1d, 1.01d, 0d, 0.1d, 0.9d, 1d, 1.1d, $
        0d, 1d, 2d, 1d, 2d, 3d, 0d, 0.5d, 1d, 1.5d]
  ff = 0.2d & beta = 0.6d & th = 0.4d
  worst = 0d
  for il = 0, 2 do begin
    w = wlist[il]
    for i = 0, 15 do $                        ; spherical: the first 16
      worst = worst > abs(limbdarkfast_hybrid_flux(bs[i], rs[i], w, laws[il]) $
                          - test_limbdarkfast_numeric(bs[i], rs[i], w, laws[il]))
    for i = 0, 19 do begin                    ; oblate: r < 1
      if rs[i] ge 1d then continue
      X = bs[i] * cos(beta) & Y = bs[i] * sin(beta)
      worst = worst > abs(limbdarkfast_hybrid_flux_oblate(X, Y, rs[i], ff, th, w, laws[il]) $
                          - test_limbdarkfast_numeric_oblate(X, Y, rs[i], ff, th, w, laws[il]))
    endfor
  endfor
  w = wlist[1]
  b = abs(3d * ((dindgen(1000) + 1d) / 1000d - 0.5d))
  Fa = limbdarkfast_hybrid_flux(b, 0.1d, w, 'hybrid4')
  for i = 0, 999 do worst = worst > abs(Fa[i] - test_limbdarkfast_numeric(b[i], 0.1d, w, 'hybrid4'))
  X = -1.5d + 3d * dindgen(201) / 200d
  Fa = limbdarkfast_hybrid_flux_oblate(X, 0.3d, 0.1d, ff, th, w, 'hybrid4')
  for i = 0, 200 do worst = worst > $
    abs(Fa[i] - test_limbdarkfast_numeric_oblate(X[i], 0.3d, 0.1d, ff, th, w, 'hybrid4'))
  return, worst
end

pro test_limbdarkfast, dir=dir
  compile_opt idl2
  if n_elements(dir) eq 0 then dir = '../tests/vectors'
  tol = 1d-12
  laws = ['hybrid2', 'hybrid4', 'hybrid5']
  nws = [2, 4, 5]
  kinds = ['sphere', 'oblate']
  failures = 0
  for ik = 0, 1 do begin
    for il = 0, 2 do begin
      nw = nws[il]
      ngeo = ik eq 0 ? 2 : 5
      d = test_limbdarkfast_read(dir + '/' + kinds[ik] + '_' + laws[il] + '.txt', ngeo + nw + 1)
      npts = (size(d, /dimensions))[1]
      Wall = d[ngeo:ngeo + nw - 1, *]
      F = dblarr(npts)
      ; rows come in blocks sharing one weight vector
      start = 0L
      while start lt npts do begin
        stop_ = start
        while stop_ + 1 lt npts do begin
          if total(abs(Wall[*, stop_ + 1] - Wall[*, start])) ne 0 then break
          stop_ = stop_ + 1
        endwhile
        w = reform(Wall[*, start])
        rows = d[*, start:stop_]
        if ik eq 0 then $
          F[start:stop_] = limbdarkfast_hybrid_flux(reform(rows[0, *]), reform(rows[1, *]), w, laws[il]) $
        else $
          F[start:stop_] = limbdarkfast_hybrid_flux_oblate(reform(rows[0, *]), reform(rows[1, *]), $
              reform(rows[2, *]), reform(rows[3, *]), reform(rows[4, *]), w, laws[il])
        start = stop_ + 1
      endwhile
      ; every 25th point again, one call per point (scalar arguments)
      for k = 0L, npts - 1, 25 do begin
        w = reform(Wall[*, k])
        if ik eq 0 then Fk = limbdarkfast_hybrid_flux(d[0, k], d[1, k], w, laws[il]) $
        else Fk = limbdarkfast_hybrid_flux_oblate(d[0, k], d[1, k], d[2, k], d[3, k], d[4, k], w, laws[il])
        if ~(abs(Fk - F[k]) le tol) then F[k] = !values.d_nan      ; counts as a failure
      endfor
      err = abs(F - reform(d[ngeo + nw, *]))
      worst = max(err, /nan)
      if total(~finite(err)) gt 0 then worst = !values.d_infinity
      ok = worst le tol
      if ~ok then failures = failures + 1
      print, kinds[ik], laws[il], npts, worst, ok ? 'ok' : 'FAIL', $
          format='(A-6, 1X, A, ": ", I5, " points, max |dF| = ", E8.2, "  ", A)'
    endfor
  endfor
  ; every shape is 1 at the limb and 0 at the centre
  for il = 0, 2 do begin
    for j = 0, nws[il] - 1 do begin
      w = dblarr(nws[il])
      w[j] = 1d
      if abs(limbdarkfast_intensity(0d, w, laws[il])) gt 1d-12 || $
         abs(limbdarkfast_intensity(1d, w, laws[il]) - 1d) gt 1d-12 then begin
        print, 'shape ', j + 1, ' of ', laws[il], ' is not 1 at the limb and 0 at the centre'
        failures = failures + 1
      endif
    endfor
  endfor
  ; a call must leave the caller's arrays as they were
  xs = [-0.95d, -0.5d, 0.2d, 0.97d] & ys = [0.3d, 0.3d, 0.3d, 0.3d]
  rs = [0.1d, 0.1d, 0.1d, 0.1d] & fs = [0.2d, 0.2d, 0.2d, 0.2d] & ts = [0.5d, 0.5d, 0.5d, 0.5d]
  keep = [xs, ys, rs, fs, ts]
  void = limbdarkfast_hybrid_flux_oblate(xs, ys, rs, fs, ts, [0.3d, 0.2d, 0.1d, 0.05d], 'hybrid4')
  void = limbdarkfast_hybrid_flux(xs, rs, [0.3d, 0.2d, 0.1d, 0.05d], 'hybrid4')
  if max(abs([xs, ys, rs, fs, ts] - keep)) ne 0 then begin
    print, 'a call changed its input arrays'
    failures = failures + 1
  endif
  if failures eq 0 then print, 'law shapes, single-point calls, inputs unchanged: ok'
  worst = test_limbdarkfast_numeric_checks()
  print, worst, worst le 1d-6 ? 'ok' : 'FAIL', $
      format='("brute-force integration: max |dF| = ", E8.2, "  ", A)'
  if ~(worst le 1d-6) then failures = failures + 1
  if failures gt 0 then message, strtrim(failures, 2) + ' check(s) failed'
end
