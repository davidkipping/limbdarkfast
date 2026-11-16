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
  if failures gt 0 then message, strtrim(failures, 2) + ' check(s) failed'
  print, 'law shapes, single-point calls, inputs unchanged: ok'
end
