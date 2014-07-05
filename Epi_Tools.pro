;Epi_Tools.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function R_fit_fn, mdt, params
;
;    Called by mpfitfun.pro. Used to fit R_forced, R_free, w0, and Omega_ps
;    to streamline over time.
;

;Declare time to be a global variable, for streamline fitting
common global_vars, time_mod, m_LR

R_forced = params[0]
R_free = params[1]
w0 = params[2]
Omega_ps = params[3]
w = w0 + Omega_ps*time_mod
phase = mod2pi(mdt - m_LR*w)
dr = -R_forced*cos(mdt) - R_free*cos(phase)

return, dr
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function R_mw_fn, mdt, params
;
;    Called by mpfitfun.pro. Used to fit r[t] - a = -R_epi*cos(mt - mw) to 
;    streamline streamline at one instant of time, with
;    params[0]=R_epi and params[1]=mw.
;

R_epi = params[0]
mw = params[1]
dr = -R_epi*cos(mdt - mw)

return, dr
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function res_fn_mod, r
;
;  Calculate Kappa - epsilon_LR*w as a function of radius r, where Kappa=epicyclic
;  frequence, w=particle's doppler-shifted frequency mm_LR*(Omega - Omega_ps), 
;  epsilon_LR=+1(-1) for an inner(outer) Lindblad resonance. GM=1 in the following, 
;  so the calculated resonance position is for a massless satellite.

;Oblateness and resonance parameters are passed via global variables.
common shared_vars, Rp, J2, mm_LR, epsilon_LR, Omega_ps

;Calculate disk particles' angular velocity Omega(r) and epicyclic frequency Kappa(r).
GM = 1d
orbit_frequencies, GM, Rp, J2, r, Omega, Kappa, Nu

;Doppler shifted frequency.
w = mm_LR*(Omega - Omega_ps)

;Frequency distance from resonance.
res_fn = Kappa - epsilon_LR*w

return, res_fn
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro forcing_function, epsilon_LR, m_LR, GM, Rp, J2, mass_s, a_s, r, Psi
;
;  Calculate perturber's forcing function, from Eqns (26-27) in Hahn et al (2010)
;
;

;Calculate function f that appears in the forcing function, Eqn. (27).
beta = r/a_s
s = 0.5d
lc = lap_coeff(m_LR, s, beta)
dlc = deriv_lap_coeff(m_LR, s, beta)
f = epsilon_LR*beta*beta*dlc + 2*m_LR*beta*lc
if (m_LR eq 1) then f = f - (2*m_LR + epsilon_LR)*beta*beta

;Particle's angular velocity Omega.
orbit_frequencies, GM, Rp, J2, r, Omega, Kappa, Nu

;Forcing function Psi, Eqn. (26).
Psi = epsilon_LR*f*mass_s*r*Omega*Omega

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro forced_eccentricity, epsilon_LR, m_LR, GM, Rp, J2, mass_s, a_s, r, e_forced, $
    Psi=Psi, D=D
;
;  Calculate perturber's forcing function, from Eqn (23a) in Hahn et al (2010)
;
;

;Calculate forcing function
forcing_function, epsilon_LR, m_LR, GM, Rp, J2, mass_s, a_s, r, Psi

;Particle and satellite orbital frequencies.
orbit_frequencies, GM, Rp, J2, r, Omega, Kappa, Nu
orbit_frequencies, GM, Rp, J2, a_s, Omega_s, Kappa_s, Nu_s

;Particle's doppler-shifted frequency w_m
w_m = m_LR*(Omega - Omega_s)

;Frequency distance from resonance squared, Eqn (23b)
D = Kappa^2 - w_m^2

;Forced eccentricity, Eqn (23a)/r
e_forced = -Psi/D/r

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro Omega_Kappa_dwdt, GM, Rp, J2, J4, J6, J8, r, Omega, Kappa, dwdt
;
;  Calculate Omega, Kappa, and dwdt=Omega-Kappa. If GM=1 then frequencies are in
;  units of radians/(2*Pi orbit periods). If inputs are mks, then output units are
;  radians/seconds.
;

;Angular frequency.
Rpa = Rp/r
GM_r3 = GM/(r^3)
Omega2 = GM_r3*( 1d +    (3d/2d)*J2*(Rpa^2)  -     (15d/8d)*J4*(Rpa^4) $ 
                    +  (35d/16d)*J6*(Rpa^6)  -  (315d/128d)*J8*(Rpa^8)  )
Omega = sqrt(Omega2)

;Epicyclic frequency.
Kappa2 = GM_r3*( 1d -    (3d/2d)*J2*(Rpa^2)  +     (45d/8d)*J4*(Rpa^4) $ 
                    - (175d/16d)*J6*(Rpa^6)  + (2205d/128d)*J8*(Rpa^8)  )
Kappa = sqrt(Kappa2)

;Precession rate
dwdt = Omega - Kappa

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro precession_rate, time, w, dwdt
;
;  Calculate dw/dt.
;

Nt = n_elements(time)
dwdt = dblarr(Nt)
for tm = 1, Nt - 1 do begin 
  dw = w[tm] - w[tm - 1]
  if (dw lt -!dpi) then dw = dw + 2*!dpi
  if (dw gt  !dpi) then dw = dw - 2*!dpi
  dt = time[tm] - time[tm - 1]
  dwdt[tm] = dw/dt
endfor
dwdt[0] = dwdt[1]

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


