;+
; EXAMPLE
;
; A spherical and an oblate planet on a circular orbit. limbdarkfast computes the
; flux at given sky positions; the orbit is up to the caller. On a circular
; orbit with orbital phase phi = 2 pi (t - t0) / P, semi-major axis a and
; impact parameter b (both in stellar radii), the planet's sky position,
; X along its direction of motion, is
;
;   X = a sin(phi),   Y = -b cos(phi),   in front of the star when cos(phi) > 0.
;
;   IDL> example
;-
pro example, t=t, f_sphere=f_sphere, f_oblate=f_oblate
  compile_opt idl2
  P = 3d & t0 = 0d & a = 10d & b = 0.3d            ; days, days, R*, R*
  r = 0.1d                                         ; area-equivalent if oblate
  f = 0.2d & theta = 30d * !dtor                   ; flattening; long axis from X
  law = 'hybrid4'
  w = [0.30d, 0.20d, 0.10d, 0.05d]

  t = -0.1d + 0.2d * dindgen(2001) / 2000d
  phi = 2d * !dpi * (t - t0) / P
  X = a * sin(phi)
  Y = -b * cos(phi)
  front = cos(phi) gt 0

  f_sphere = replicate(1d, n_elements(t))
  f_oblate = replicate(1d, n_elements(t))
  i = where(front)
  f_sphere[i] = limbdarkfast_hybrid_flux(sqrt(X[i]^2 + Y[i]^2), r, w, law)
  f_oblate[i] = limbdarkfast_hybrid_flux_oblate(X[i], Y[i], r, f, theta, w, law)

  print, law, ', w = ', w, format='(A, A, 4F6.2)'
  print, 'depth at mid-transit: sphere ', 1d - f_sphere[1000], ', oblate ', 1d - f_oblate[1000], $
      format='(A, E13.6, A, E13.6)'
  print, 'largest oblate - sphere difference: ', max(abs(f_oblate - f_sphere)) * 1d6, ' ppm', $
      format='(A, F6.1, A)'
end
