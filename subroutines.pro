;subroutines.pro
;
;    These functions and procedures are called by v7.0 of epi_int.pro.
;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mod2pi, angle
;
;  Force angle to lie between -Pi < angle < Pi radians
;
twopi = 2*!dpi
angle_mod = angle mod twopi
j = where(angle_mod gt !dpi, Nj)
if (Nj gt 0) then angle_mod[j] = angle_mod[j] - twopi
j = where(angle_mod lt -!dpi, Nj)
if (Nj gt 0) then angle_mod[j] = angle_mod[j] + twopi

return, angle_mod
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro initialize, inputs, quantities, time, ring, satellite, coordinates
;
;  Initialize
;

;Declare the global variables that are used to solve for the resonance radius.
common shared_vars, Rp, J2, m_LR, epsilon_LR, r_sat

;Start time and storage for run time, in minutes.
start_time_min = systime(/seconds)/60d
run_time_min = 0d

;Extract input parameters.
restart = inputs.restart
Nr = inputs.Nr
Nt = inputs.Nt
a_s = inputs.a_s
e_s = inputs.e_s
I_s = inputs.I_s
w_s = inputs.w_s
M_s = inputs.M_s
N_s = inputs.N_s
mass_s_final = inputs.mass_s_final
t_grow_sat = inputs.t_grow_sat
sat_index = inputs.sat_index
m_LR = inputs.m_LR
fourier_satellite = inputs.fourier_satellite
ringlet_disk = inputs.ringlet_disk
a_r_min = inputs.a_r_min
a_r_max = inputs.a_r_max
radial_spacing = inputs.radial_spacing
da1 = inputs.da1
e_r_init = inputs.e_r_init
de_r_init = inputs.de_r_init
I_r_init = inputs.I_r_init
w_r_init = inputs.w_r_init
N_r_init = inputs.N_r_init
dt = inputs.dt
Noutput = inputs.Noutput
Nsteps_per_output = inputs.Nsteps_per_output
output_file = inputs.output_file
M_planet_gm = inputs.M_planet_gm
r1_cm = inputs.r1_cm
Rp_cm = inputs.Rp_cm
J2 = inputs.J2
surface_density = inputs.surface_density
shear_viscosity = inputs.shear_viscosity
bulk_viscosity = inputs.bulk_viscosity
optical_depth = inputs.optical_depth
c0 = inputs.c0
display = inputs.display
precision = inputs.precision

;Streamline semimajor axes, where streamline spacings can be uniform, or the
;spacings can grow exponentially radially outwards or radially inwards.
a_vector = [a_r_min]
if (Nr gt 1) then begin 
    if (radial_spacing eq 'uniform') then begin  
        a_vector = (a_r_max - a_r_min)*dindgen(Nr)/(Nr - 1) + a_r_min
    endif
    if (radial_spacing eq 'expo-outwards') then begin  
        Delta = a_r_max - a_r_min
        factor = alog(Delta/da1)/(Nr - 2)
        da = da1*exp(factor*(dindgen(Nr) - 1))
        da[0] = 0d
        a_vector = a_r_min + da
    endif
    if (radial_spacing eq 'expo-inwards') then begin  
        Delta = a_r_max - a_r_min
        factor = alog(Delta/da1)/(Nr - 2)
        da = da1*exp(factor*(dindgen(Nr) - 1))
        da[0] = 0d
        da = rotate(da, 2)
        a_vector = a_r_max - da
    endif
endif
a_r = dblarr(Nt, Nr) 
for j = 0, Nr - 1 do a_r[*, j] = a_vector[j]

;Initial ring eccentricities, inclinations, longitude of periapse, & ascending nodes,
;for a narrow ringlet or a disk that has a Rayleigh distribution of e's and I's.
if (ringlet_disk eq 'ringlet') then begin
    e_r = dblarr(Nt, Nr)
    for j = 0, Nr - 1 do e_r[*, j] = e_r_init + j*de_r_init/(Nr - 1) 
    I_r = dblarr(Nt, Nr) + I_r_init
    w_r = dblarr(Nt, Nr)
    N_r = dblarr(Nt, Nr)
endif
if (ringlet_disk eq 'disk') then begin
    rn = randomu(seed, Nt, Nr)
    e_r = e_r_init*sqrt(-alog(1d - rn))
    rn = randomu(seed, Nt, Nr)
    I_r = I_r_init*sqrt(-alog(1d - rn))
    w_r = dblarr(Nt, Nr)
    N_r = dblarr(Nt, Nr)
endif

;Ring mean anomaly.
M_vector = [0d]
if (Nt gt 1) then M_vector = !dpi*(2*dindgen(Nt)/Nt - 1)
M_r = dblarr(Nt, Nr)
for j = 0, Nt - 1 do M_r[j, *] = M_vector[j]

;Planet radius in units of r1.
Rp = Rp_cm/r1_cm

;Total number of ring particles.
Np = long(Nt)*long(Nr)

;r1=1 in units of km.
r1_km = r1_cm/1d5

;Mass of each ring particle, assuming a circular ringlet
delta = (shift(a_r, 0, -1) - shift(a_r, 0, 1))/2
delta[0, 0] = a_r[*, 1] - a_r[*, 0]
delta[0, Nr - 1] = a_r[*, Nr - 1] - a_r[*, Nr - 2]
mass_r = 2*!dpi*surface_density*a_r*delta/Nt

;Set the satellites' initial mass, which starts at zero if growing,
;or if one is a fourier satellite.
mass_s = mass_s_final
Ns = n_elements(a_s)
if (t_grow_sat gt 0) then mass_s = dblarr(Ns)
if (fourier_satellite eq 'y') then mass_s[sat_index] = 0

;Give each ring particle in a streamline a unique ID number so they can be tracked over time.
ID_r = intarr(Nt, Nr)
for j = 0, Nt - 1 do ID_r[j, *] = j

;Convert ring and satellite orbit elements to planetocentric polar coordinates.
GM_r = 1d + mass_r
elements2coordinates, a_r, e_r, I_r, w_r, M_r, N_r, GM_r, Rp, J2, precision, $
    r_r, t_r, z_r, vr_r, vt_r, vz_r
GM_s = 1d + mass_s
elements2coordinates, a_s, e_s, I_s, w_s, M_s, N_s, GM_s, Rp, J2, precision, $
    r_s, t_s, z_s, vr_s, vt_s, vz_s

;Orbit period in seconds, days, and years at r=r1. G_cgs=gravitation constant in cgs units.
twopi = 2*!dpi
G_cgs = 6.67259d-8
T1_sec = twopi*sqrt((r1_cm^3)/G_cgs/M_planet_gm)
T1_day = T1_sec/60/60/24
T1_yr = T1_day/365.25d0
print, 'orbital period at r=1 (sec)                            = ', T1_sec

;ring surface density in gm/cm^2
surface_density_cgs = surface_density*M_planet_gm/(r1_cm^2)
print, 'ring surface density (gm/cm^2)                         = ', surface_density_cgs

;Viscosity in cm^2/sec.
shear_viscosity_cgs =  twopi*shear_viscosity*(r1_cm^2)/T1_sec
bulk_viscosity_cgs  =  twopi*bulk_viscosity*(r1_cm^2)/T1_sec

;Toomre wavelength at r=1 assuming G=1
r1 = 1d
GM = 1d
orbit_frequencies, GM, Rp, J2, r1, Omega_1, Kappa_1, Eta_1, Beta_1, Nu_1
lambda_T = 4*!dpi*!dpi*surface_density/Kappa_1^2
print, 'Toomre wavelength                                      = ', lambda_T

;Toomre stability parameter assuming G=1
Q_T = -1d
if (surface_density gt 0) then Q_T = c0*Kappa_1/!dpi/surface_density
print, 'Toomre Q                                               = ', Q_T

;Density wave group speed and crossing time assuming G=1.
v_wave = !dpi*surface_density/Kappa_1
t_wave = -1d
if (v_wave gt 0) then t_wave = (a_r_max - a_r_min)/v_wave
print, 'density wave crossing time                             = ', t_wave

;Ring opacity.
opacity = -1d
if (surface_density gt 0) then opacity = optical_depth/surface_density

;Find the resonance radius r_res. The satellite's approximate resonance radius is
;r_res_approx, and the array r_guess brackets the resonance radius.
sat_index = inputs.sat_index
epsilon_LR = inputs.epsilon_LR
m_LR = inputs.m_LR
r_sat = a_s[sat_index]
r_res_approx = ( ( 1d - double(epsilon_LR)/m_LR )^(2d/3d) )*r_sat
if ((epsilon_LR eq 1) and (m_LR eq 1)) then $ 
    ;The epsilon_LR=1=m_LR is where the particle's precession rate=satellite's angular speed.
    r_res_approx = r_sat*((3*J2/2)*(Rp/r_sat)^2)^(2d/7d)
r_guess = [0.8d, 1 + 1d-8, 1.2d]*r_res_approx
r_res = fx_root( r_guess, 'resonance_fn', /double, tol=1.0e-10)
print, 'resonance radius                                       = ', r_res

;Calculate the laplace coefficients that are used by the fourier satellite. Note that
;this assumes that the ring particles' semimajor axes stay constant or nearly so.
beta = a_r/a_s[sat_index]
s = 0.5d
lc = lap_coeff(m_LR, s, beta)
dlc = deriv_lap_coeff(m_LR, s, beta)

;Calculate the drag coefficient across the ring, even if not needed.
drag_coeff_inner = inputs.drag_coeff_inner
drag_coeff_outer = inputs.drag_coeff_outer
drag_transition_radius = inputs.drag_transition_radius
drag_transition_width = inputs.drag_transition_width
drag_timescale = inputs.drag_timescale
r_ref = inputs.r_ref
da = a_r - r_ref
s = 1d
drag_coeff_min = drag_coeff_inner
if (drag_coeff_outer lt drag_coeff_inner) then begin
    s = -1d
    drag_coeff_min = drag_coeff_outer
endif
drag_coeff = atan(s*(da - drag_transition_radius)/drag_transition_width)/!dpi + 0.5d
drag_coeff = drag_coeff - min(drag_coeff)
drag_coeff = drag_coeff/max(drag_coeff)
drag_coeff = drag_coeff*s*(drag_coeff_outer - drag_coeff_inner) + drag_coeff_min

;Increase number of outputs +1, for storing initial state.
Noutput = Noutput + 1

;Set initial time and create storage for time and system's total angular momentum
t_now = 0d
time = dblarr(Noutput)
tm = 0
time[tm] = t_now
t_start = t_now
ang_mom = dblarr(Noutput)

;Store various quantities in structure.
qt = { Np:Np, G_cgs:G_cgs, T1_sec:T1_sec, T1_day:T1_day, T1_yr:T1_yr, $
    r1_km:r1_km, surface_density_cgs:surface_density_cgs, lambda_T:lambda_T, Q_T:Q_T, $
    shear_viscosity_cgs:shear_viscosity_cgs, bulk_viscosity_cgs:bulk_viscosity_cgs, $
    start_time_min:start_time_min, run_time_min:run_time_min, Rp:Rp, Ns:Ns, $
    opacity:opacity, r_res:r_res, beta:beta, lc:lc, dlc:dlc, drag_coeff:drag_coeff, $
    t_start:t_start, v_wave:v_wave, t_wave:t_wave, ang_mom:ang_mom}
quantities = create_struct(inputs, qt)
quantities.Noutput = Noutput

;Nudge the ring's tangential velocity so that it is in centrifugual
;equilibrium when the additional forces on the ring are accounted for.
accelerations, quantities, t_now, mass_r, r_r, t_r, z_r, vr_r, vt_r, vz_r, $
    mass_s, r_s, t_s, z_s, Ar_r, At_r, Az_r, Ar_s, At_s, Az_s, $
    shear_viscosity=shear_viscosity, bulk_viscosity=bulk_viscosity
for j = 0, Nr - 1 do Ar_r[*, j] = mean(Ar_r[*, j])
vt_r = sqrt(vt_r^2 - r_r*Ar_r)

;Update ring viscosity in cm^2/sec, for when viscosity_law is not 'constant'.
shear_viscosity_cgs =  twopi*median([shear_viscosity])*(r1_cm^2)/T1_sec
bulk_viscosity_cgs  =  twopi*median([bulk_viscosity])*(r1_cm^2)/T1_sec
print, 'ring shear viscosity (cm^2/sec)                        = ', shear_viscosity_cgs
print, 'ring bulk  viscosity (cm^2/sec)                        = ', bulk_viscosity_cgs
quantities.shear_viscosity_cgs = shear_viscosity_cgs
quantities.bulk_viscosity_cgs = bulk_viscosity_cgs

;Convert planetocentric velocities to barycentric velocities.
planeto_rv2bary_v, mass_r, r_r, t_r, z_r, vr_r, vt_r, vz_r, mass_s, r_s, t_s, z_s, $
    vr_s, vt_s, vz_s, vr_r_BC, vt_r_BC, vz_r_BC, vr_s_BC, vt_s_BC, vz_s_BC

;Henceforth all velocities are barycentric.
vr_r = vr_r_BC
vt_r = vt_r_BC
vz_r = vz_r_BC
vr_s = vr_s_BC
vt_s = vt_s_BC
vz_s = vz_s_BC

;Structure contains the ring's current planetocentric coordinates, barycentric
;velocities, mass, and ID.
ring = create_struct('r', r_r, 't', t_r, 'z', z_r, 'vr', vr_r, 'vt', vt_r, $
    'vz', vz_r, 'mass', mass_r, 'ID', ID_r)

;Structure contains the satellites' current planetocentric coordinates, barycentric
;velocities, and mass.
satellite = create_struct('r', r_s, 't', t_s, 'z', z_s, 'vr', vr_s, 'vt', vt_s, $
    'vz', vz_s, 'mass', mass_s)

;If restarting, restore final state from earlier run to use as initial conditions.
;Also use satellite masses from inputs.pro rather than preserving their earlier masses.
if (restart eq 'y') then begin
    restore, filename='restart.dat' 
    time[tm] = t_now
    quantities.t_start = t_now
    satellite.mass = mass_s
endif

;Create an array of structures for storing the ring and satellites' planetocentric
;coordinates and velocities, masses, and ring ID over time.
arr_r = dblarr(Nt, Nr)
arr_s = dblarr(Ns)
structr = create_struct('r_r', r_r, 't_r', t_r, 'z_r', z_r, 'vr_r', vr_r, 'vt_r', vt_r, $
    'vz_r', vz_r, 'mass_r', mass_r, 'ID_r', ID_r, 'r_s', r_s, 't_s', t_s, 'z_s', z_s, $
    'vr_s', vr_s, 'vt_s', vt_s, 'vz_s', vz_s, 'mass_s', mass_s)
coordinates = replicate(structr, Noutput)

;Store system's initial state.
store_results, tm, t_now, time, quantities, ring, satellite, coordinates

;Display initial state.
spawn, 'pwd', ttl
if (display eq 'monitor') then window, xs=900, ys=900, retain=2, title=ttl[0]
if (display eq 'postscript') then begin
    if (!d.window gt -1) then wdelete, !d.window	;delete current plotting window
    ft = file_test('postscript')			;create postscript directory
    if (ft ne 1) then file_mkdir, 'postscript'		;as needed.
endif
display_results, inputs, quantities, tm, time, coordinates

;set IDL's !ERROR_STATE system variable to its default no-errors state
message, /reset

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro orbit_frequencies, GM, Rp, J2, a, Omega, Kappa, Eta, Beta, Nu, $
    Omega2=Omega2, Kappa2=Kappa2, Eta2=Eta2, Beta2=Beta2, Nu2=Nu2
;
;  Calculate mean angular velocity Omega, epicyclic frequency Kappa, and
;  vertical oscillation frequency Nu expanded out to the J2 oblateness term,
;  as a function of distance a. These frequencies are what B-R&L(1994) call
;  Omega_0, Kappa_0, Eta_0, and Nu_0, Eqns. (7,16,17,30). Optional parameters
;  return these frequencies squared.

a2 = a*a
GM_a3 = GM/a2/a					;GM/distance^3
factor = J2*Rp*Rp/a2				;J2*(Rp/a)^2
Omega2 = GM_a3*(1 + factor*1.5d)		;angular frequency^2
Omega = sqrt(Omega2)				;angular frequency
Kappa2 = GM_a3*(1 - factor*1.5d)		;epicyclic frequency^2
Kappa = sqrt(Kappa2)				;epicyclic frequency
Eta2 = GM_a3*(1 - factor*2d)			;eta frequency^2
Eta = sqrt(Eta2)				;eta frequency
Beta2 = GM_a3*(1 + factor*7.5d)			;beta frequency^2
Beta = sqrt(Beta2)				;beta frequency
Nu2 = GM_a3*(1 + factor*4.5d)			;vertical frequency^2
Nu = sqrt(Nu2)					;vertical frequency

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro planet_potential, GM, Rp, J2, r, z, Phi
;
;  The oblate planet's gravitational potential Phi(r,z).
;

rho = sqrt(r*r + z*z)				;distance to planet
Rp_rho = Rp/rho					;planet radius/distance
sin_lat = z/rho					;sin(latitude)
P2 = (3*sin_lat*sin_lat - 1)/2			;second Lagrange polynomial of sin_lat.
Phi = (-GM/rho)*(1 - J2*P2*Rp_rho*Rp_rho)	;potential

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro elements2coordinates, a, e, I, w, M, N, GM, Rp, J2, precision, $
    r, t, z, vr, vt, vz, Omega=Omega, Kappa=Kappa, Eta=Eta, Nu=Nu
;
;  Convert orbit elements to cylindrical coordinates and velocities to first or second
;  order in e and first order in I. Also returns the orbit frequencies as optional parameters. 
;

orbit_frequencies, GM, Rp, J2, a, Omega, Kappa, Eta, Beta, Nu
OK = Omega/Kappa
x = e*cos(M)
y = e*sin(M)
if (precision eq 1) then begin 
    r = a*(1 - x)
    t = mod2pi( w + M + OK*y*2 )
    vr = a*Kappa*y 			
    vt = a*Omega*(1 + x)
endif
if (precision eq 2) then begin 
    EK = Eta/Kappa
    EK2 = EK*EK
    x2 = x*x
    e2 = e*e
    r = a*(1 - x + EK2*(2*e2 - x2))
    t = mod2pi( w + M + OK*y*( 2 + (1.5d + EK2)*x ) )
    vr = a*Kappa*y*(1 + 2*EK2*x) 			
    vt = a*Omega*(1 + x - 2*EK2*e2 + (1 + EK2)*x2)
endif
aI = a*I
phi = w - N + M
z = aI*sin(phi)
vz = aI*Nu*cos(phi)

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro coordinates2elements, r, t, z, vr, vt, vz, GM, Rp, J2, precision, $
    a, e, I, w, M, N, Omega=Omega, Kappa=Kappa, Eta=Eta, Beta=Beta, Nu=Nu
;
;  Convert polar coordinates to orbit elements and return optional orbit frequencies
;  to accuracy e^2 and I^1.
;

;Get semimajor axis a, Energy, potential Phi(a), and orbit frequencies assuming z=0.
zero = 0d
semimajor_axis, r, zero, vr, vt, zero, GM, Rp, J2, a, Energy=Energy
planet_potential, GM, Rp, J2, a, zero, Phi_a
orbit_frequencies, GM, Rp, J2, a, Omega, Kappa, Eta, Beta, Nu

;Various factors
OK = Omega/Kappa
OK2 = OK*OK
aOmega = a*Omega
aK = a*Kappa
vr_factor = vr/aK
ra = r/a
Delta = r - a

;Get mean anomaly M, longitude of periapse w, ascending node N, inclination I
;where x=e*cos(M), y=e*sin(M), xi=I*cos(phi), yi=I*sin(phi) where phi=w-N+M.
;To get e, use Eqn (79) in B-R&L1995.
if (precision eq 1) then begin
    twoI3 = (vr*vr + (Kappa*Delta)^2 )
    e = sqrt(twoI3>0)/aK
    x = 1 - ra
    y = vr_factor
    M = atan(y, x)
    w = mod2pi( t - M - y*OK*2 )
endif
if (precision eq 2) then begin
    twoI3 = (vr*vr + (Kappa*Delta)^2 ) - 2*Eta*Eta*(Delta^3)/a
    e = sqrt(twoI3>0)/aK
    EK = Eta/Kappa
    EK2 = EK*EK
    vt_factor = vt/aOmega
    D = 1 - ra
    x = EK2*( (1 - vt_factor) + D + 2*e*e ) + D
    y = vr_factor/(1 + 2*EK2*x)
    M = atan(y, x)
    w = mod2pi( t - M - y*OK*( 2 + (1.5d + EK2)*x ) )
endif
xi = vz/a/Nu
yi = z/a
phi = atan(yi, xi)
N = mod2pi(w + M - phi)
I = sqrt(xi*xi + yi*yi)

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro semimajor_axis, r, z, vr, vt, vz, GM, Rp, J2, a, Energy=Energy
;
;  Solve for the body's semimajor axis a.
;

;Gravitational potential Phi and two-body energy.
planet_potential, GM, Rp, J2, r, z, Phi
Energy = (vr*vr + vt*vt + vz*vz)/2 + Phi

;Semimajor axis for orbit about oblate planet assuming GM=1.
h = r*vt
g = h*h/2/Rp
g2 = g*g
a = g*Rp*(   1 + sqrt( 1 - (1.5d)*J2/g2 )   )

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function resonance_fn, r
;
;  Calculate Kappa - epsilon_LR*w as a function of radius r, where Kappa=epicyclic
;  frequence, w=particle's doppler-shifted frequency, epsilon_LR=+1(-1) for an
;  inner(outer) Lindblad resonance. GM=1 in the following, so the calculated
;  resonance position is for a massless satellite.

;Oblateness and resonance parameters are passed via global variables.
common shared_vars, Rp, J2, m_LR, epsilon_LR, r_sat

;Calculate disk particles' angular velocity Omega(r) and epicyclic frequency Kappa(r).
GM = 1d
orbit_frequencies, GM, Rp, J2, r, Omega, Kappa, Eta, Beta, Nu

;Calculate satellite's angular velocity Omega_sat.
orbit_frequencies, GM, Rp, J2, r_sat, Omega_sat, Kappa_sat, Eta_sat, Beta_sat, Nu_sat

;Doppler shifted frequency.
w = m_LR*(Omega - Omega_sat)

;Frequency distance from resonance.
res_fn = Kappa - epsilon_LR*w

return, res_fn
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro planeto_rv2bary_v, mass_r, r_r, t_r, z_r, vr_r, vt_r, vz_r, mass_s, r_s, t_s, z_s, $
    vr_s, vt_s, vz_s, vr_r_BC, vt_r_BC, vz_r_BC, vr_s_BC, vt_s_BC, vz_s_BC
;
;  Calculate system's barycentric velocities from planetocentric coordinates
;  and velocities. Inputs are ring mass and polar coordinates and velocities,
;  and satellite mass and polar coordinates and velocities,
;  while outputs are barycentric polar velocities for ring and satellites.

;Convert ring's planetocentric cylindrical coordinates to planetocentric cartesian.
rt2xy, r_r, t_r, vr_r, vt_r, x_r_pc, y_r_pc, vx_r_pc, vy_r_pc
z_r_pc = z_r
vz_r_pc = vz_r

;Convert satellites' planetocentric cylindrical coordinates to planetocentric cartesian.
rt2xy, r_s, t_s, vr_s, vt_s, x_s_pc, y_s_pc, vx_s_pc, vy_s_pc
z_s_pc = z_s
vz_s_pc = vz_s

;Calculate barycentric cartesian coordinates.
planeto2bary_xyz, mass_r, x_r_pc, y_r_pc, z_r_pc, mass_s, x_s_pc, y_s_pc, z_s_pc, $
    x_r_BC, y_r_BC, z_r_BC, x_s_BC, y_s_BC, z_s_BC, x_0, y_0, z_0 

;Calculate barycentric cartesian velocities.
planeto2bary_xyz, mass_r, vx_r_pc, vy_r_pc, vz_r_pc, mass_s, vx_s_pc, vy_s_pc, vz_s_pc, $
    vx_r_BC, vy_r_BC, vz_r_BC, vx_s_BC, vy_s_BC, vz_s_BC, vx_0, vy_0, vz_0 

;Convert ring and satellite's barycentric cartesian coordinates and velocities to 
;barycentric polar coordinates and velocities.
xy2rt, x_r_BC, y_r_BC, vx_r_BC, vy_r_BC, r_r_BC, t_r_BC, vr_r_BC, vt_r_BC
xy2rt, x_s_BC, y_s_BC, vx_s_BC, vy_s_BC, r_s_BC, t_s_BC, vr_s_BC, vt_s_BC

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mixed_rv2planeto_v, mass_r, r_r, t_r, z_r, vr_r_BC, vt_r_BC, vz_r_BC, $
    mass_s, r_s, t_s, z_s, vr_s_BC, vt_s_BC, vz_s_BC, vr_r_pc, vt_r_pc, vz_r_pc, $
    vr_s_pc, vt_s_pc, vz_s_pc
;
;  Calculate system's planetocentric velocities from planetocentric coordinates
;  and barycentric velocities. Inputs are ring mass and polar coordinates and velocities,
;  and satellite mass and polar coordinates and velocities,
;  while outputs are planetocentric polar velocities for ring and satellites.

;Convert ring and satellite's planetocentric cylindrical coordinates to
;planetocentric cartesian.
x_r_pc = r_r*cos(t_r)
y_r_pc = r_r*sin(t_r)
z_r_pc = z_r
x_s_pc = r_s*cos(t_s)
y_s_pc = r_s*sin(t_s)
z_s_pc = z_s

;Calculate system's barycentric cartesian coordinates and longitudes.
planeto2bary_xyz, mass_r, x_r_pc, y_r_pc, z_r_pc, mass_s, x_s_pc, y_s_pc, z_s_pc, $
    x_r_BC, y_r_BC, z_r_BC, x_s_BC, y_s_BC, z_s_BC, x_0, y_0, z_0 
t_r_BC = atan(y_r_BC, x_r_BC)
t_s_BC = atan(y_s_BC, x_s_BC)

;Convert ring and satellite's barycentric polar velocities to barycentric cartesian velocities.
vx_r_BC = vr_r_BC*cos(t_r_BC) - vt_r_BC*sin(t_r_BC)
vy_r_BC = vr_r_BC*sin(t_r_BC) + vt_r_BC*cos(t_r_BC)
vx_s_BC = vr_s_BC*cos(t_s_BC) - vt_s_BC*sin(t_s_BC)
vy_s_BC = vr_s_BC*sin(t_s_BC) + vt_s_BC*cos(t_s_BC)

;Calculate planetocentric cartesian velocities.
bary2planeto_xyz, mass_r, vx_r_BC, vy_r_BC, vz_r_BC, mass_s, vx_s_BC, vy_s_BC, vz_s_BC, $
    vx_r_pc, vy_r_pc, vz_r_pc, vx_s_pc, vy_s_pc, vz_s_pc, vx_0, vy_0, vz_0 

;Convert ring and satellite's planetocentric cartesian coordinates and velocities to 
;planetocentric polar coordinates and velocities.
xy2rt, x_r_pc, y_r_pc, vx_r_pc, vy_r_pc, r_r_pc, t_r_pc, vr_r_pc, vt_r_pc
xy2rt, x_s_pc, y_s_pc, vx_s_pc, vy_s_pc, r_s_pc, t_s_pc, vr_s_pc, vt_s_pc

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro planeto2bary_xyz, mass_r, x_r_pc, y_r_pc, z_r_pc, mass_s, x_s_pc, y_s_pc, z_s_pc, $
    x_r_BC, y_r_BC, z_r_BC, x_s_BC, y_s_BC, z_s_BC, x_0, y_0, z_0 
;
;  Convert the planetocentric cartesian coordinates of ring (x_r_pc,y_r_pc,z_r_pc)
;  and satellites (x_s_pc,y_s_pc,z_s_pc) into barycentric coordinates for ring
;  and satellite (x_r_BC, etc). Other inputs are the ring and satellite
;  masses mass_r,mass_s. Other outputs are the central planet's 
;  coordinates x_0,y_0,z_0 relative to the barycenter. Replace coordinates with
;  velocities to convert convert planetocentric velocities into barycentric
;  velocities.

;Calculate the central planet's position relative to the barycenter.
mass_0 = 1						;central mass is 1
mass_t = mass_0 + total(mass_r) + total(mass_s)		;total mass
x_0 = (-1/mass_t)*( total(mass_r*x_r_pc) + total(mass_s*x_s_pc) )
y_0 = (-1/mass_t)*( total(mass_r*y_r_pc) + total(mass_s*y_s_pc) )
z_0 = (-1/mass_t)*( total(mass_r*z_r_pc) + total(mass_s*z_s_pc) )

;Calculate barycentric coordinates of the disk and satellites.
x_r_BC = x_r_pc + x_0
y_r_BC = y_r_pc + y_0
z_r_BC = z_r_pc + z_0
x_s_BC = x_s_pc + x_0
y_s_BC = y_s_pc + y_0
z_s_BC = z_s_pc + z_0

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro bary2planeto_xyz, mass_r, x_r_bc, y_r_bc, z_r_bc, mass_s, x_s_bc, y_s_bc, z_s_bc, $
    x_r_PC, y_r_PC, z_r_PC, x_s_PC, y_s_PC, z_s_PC, x_0, y_0, z_0
;
;  Convert barycentric cartesian coordinates for ring (x_r_bc,y_r_bc,z_r_bc) and
;  satellites (x_s_bc,y_s_bc,z_s_bc) into planetocentric coordinates for disk and
;  satellites (x_r_PC, etc). Other inputs are the ring and satellite masses
;  mass_r,mass_s, and other outputs are the central planet's coordinates x_0,y_0,z_0
;  relative to the barycenter. Replace coordinates with velocities to convert
;  barycentric velocities into planetocentric velocities.

;Calculate the central planet's position relative to the barycenter.
mass_0 = 1						;central mass is 1
x_0 = (-1/mass_0)*( total(mass_r*x_r_bc) + total(mass_s*x_s_bc) )
y_0 = (-1/mass_0)*( total(mass_r*y_r_bc) + total(mass_s*y_s_bc) )
z_0 = (-1/mass_0)*( total(mass_r*z_r_bc) + total(mass_s*z_s_bc) )

;Calculate planetocentric velocities of the disk and satellites.
x_r_PC = x_r_bc - x_0
y_r_PC = y_r_bc - y_0
z_r_PC = z_r_bc - z_0
x_s_PC = x_s_bc - x_0
y_s_PC = y_s_bc - y_0
z_s_PC = z_s_bc - z_0

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rt2xy, r, t, vr, vt, x, y, vx, vy
;
;  Convert cylindrical coordinates and velocities r,t,vr,vt to cartesian x,y,vx,vy.
;  Optional arguments provide the sin and cos of longitude t.
;

sint = sin(t)
cost = cos(t)
x = r*cost
y = r*sint
vx = vr*cost - vt*sint
vy = vr*sint + vt*cost

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro xy2rt, x, y, vx, vy, r, t, vr, vt
;
;  Convert cartesian coordinates and velocities x,y,vx,vy to cylindrical r,t,vr,vt.

r = sqrt(x*x + y*y)
t = atan(y, x)
sint = sin(t)
cost = cos(t)
vr =  vx*cost + vy*sint
vt = -vx*sint + vy*cost

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sort_ring_longitudes, r, t, z, vr, vt, vz, mass, ID
;
;  Sort ring particles within each streamline j so that all particle longitudes t_d[idx, j]
;  increase with index idx.
;

;Extract array dimensions, noting that r is a 1D vector if Nr=1 streamline.
sz = size(r)
Ndim = sz[0]
Nt = sz[1]
Nr = sz[2]
if (Ndim eq 1) then Nr = 1

;Sort particles in each streamline according to their longitude.
if (Nt gt 1) then for j = 0, Nr - 1 do begin 
    s = sort(t[*, j])
    r[0, j] = r[s, j]
    t[0, j] = t[s, j]
    z[0, j] = z[s, j]
    vr[0, j] = vr[s, j]
    vt[0, j] = vt[s, j]
    vz[0, j] = vz[s, j]
    mass[0, j] = mass[s, j]
    ID[0, j] = ID[s, j]
endfor

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro longitude_wrap, t, f, tw, fw, in_array=in_array, wrapped_array=wrapped_array
;
;  Wrap ring particle arrays around in longitude. Inputs t,f are 2D arrays
;  where t=ring longitudes and f=some function of ring longitude (such as the
;  ring radius array). Array fw = copy of the input f array but
;  `wrapped' around in longitude such that f's rightmost column is also copied
;  and appended to the left side of fw, and f's leftmost column is copied and
;  appended to fw's right side. Ditto for the t array except the column appended
;  to the left has 2*Pi subtracted to it, while the appended column on the right
;  has twopi added to it. The optional in_array also gets wrapped as well.
;

sz = size(t)
Nt = sz[1]
Nr = sz[2]
twopi = 2*!dpi
tw = [ t[Nt - 1, *] - twopi, t, t[0, *] + twopi ]
fw = [ f[Nt - 1, *], f, f[0, *] ]
if (keyword_set(in_array) eq 1) then $
    wrapped_array = [ in_array[Nt - 1, *], in_array, in_array[0, *] ]

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro unwrap_angle, a, a_unw
;
;  Unwrap angle a
;
twopi = 2*!dpi
N = n_elements(a)
a_unw = a
for j = 1, N - 1 do begin 
    delta = a[j] - a[j - 1]
    if (delta gt !dpi) then delta = delta - twopi
    if (delta lt -!dpi) then delta = delta + twopi
    a_unw[j] = a_unw[j - 1] + delta
endfor

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro interpolate_f, t, f, n, interpolation_method, f_n
;
;  Inputs are t=ring longitude and f=ring property evaluate at each ring particle's
;  longitude t. Let t[j,k]=longitude theta for particle j in streamline k. The following
;  uses various interpolation methods to evaluate f at longitude t[j,k] along 
;  streamline k+n, with results returned in array f_n. 
;  When interpolation_method='approximate' then f_n[j,k]=f[j,k+n] (ie, array f is
;  simply shifted -n rows). When interpolation_method='polynomial' then a
;  second degree lagrange polynomial is fitted to points f[j-1:j+1,k+n], with that
;  fit then used to evaluate f at longitude t[j,k] along streamline k+n.

;extract array dimensions
sz = size(t)
Nt = sz[1]
Nr = sz[2]

case interpolation_method of

    ;Least accurate but fastest method.
    'approximate': f_n = shift(f, 0, -n)

    ;Use a 2nd degree Lagrange polynomial to interpolate. More accurate but slower. 
    'polynomial': begin
        longitude_wrap, t, f, tw, fw
        t0 = shift(tw,  1, -n)
        t1 = shift(tw,  0, -n)
        t2 = shift(tw, -1, -n)
        f0 = shift(fw,  1, -n)
        f1 = shift(fw,  0, -n)
        f2 = shift(fw, -1, -n)
        lagrange_poly_fit, t0, t1, t2, f0, f1, f2, tw, fw_fit
        f_n = fw_fit[1:Nt, *]
    end

endcase

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lagrange_poly_fit, x0, x1, x2, y0, y1, y2, x, y
; 
;  Calculate the 2nd order lagrange polynomial fit to data (x0,y0),(x1,y1),(x2,y2)
;  and use that polynomial to calculate y(x).
;
l0 = ((x - x1)/(x0 - x1))*((x - x2)/(x0 - x2))
l1 = ((x - x0)/(x1 - x0))*((x - x2)/(x1 - x2))
l2 = ((x - x0)/(x2 - x0))*((x - x1)/(x2 - x1))
y = y0*l0 + y1*l1 + y2*l2

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro derivative_dfdr, r, t, f, difference_method, dfdr
;
;  Calculate df/dr via finite differences and interpolation. 
;

;df/dr in ring ring interior.
interpolate_f, t, r,  1, difference_method, r_plus
interpolate_f, t, r, -1, difference_method, r_minus
interpolate_f, t, f,  1, difference_method, f_plus
interpolate_f, t, f, -1, difference_method, f_minus
dfdr = (f_plus - f_minus)/(r_plus - r_minus)

;At inner edge.
dfdr[0, 0] = (f_plus[*, 0] - f[*, 0])/(r_plus[*, 0] - r[*, 0])

;At outer edge.
sz = size(r)
Nt = sz[1]
Nr = sz[2]
j = Nr - 1
dfdr[0, j] = (f[*, j] - f_minus[*, j])/(r[*, j] - r_minus[*, j])

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro radiative_bc, rbc_indices, quantities, r, t, z, vr, vt, vz
;
;  Used to allow a spiral density wave to propagate out of the simulated ring region,
;  rather than reflecting at the ring's edge. Need to explain how this works...
;

;Extract fitting annulus
idx_fit_in = rbc_indices[0]
idx_fit_out = rbc_indices[1]
r_fit = r[*, idx_fit_in:idx_fit_out]
t_fit = t[*, idx_fit_in:idx_fit_out]
z_fit = z[*, idx_fit_in:idx_fit_out]
vr_fit = vr[*, idx_fit_in:idx_fit_out]
vt_fit = vt[*, idx_fit_in:idx_fit_out]
vz_fit = vz[*, idx_fit_in:idx_fit_out]

;Get orbit elements of fitting annulus
GM = 1d
Rp = quantities.Rp
J2 = quantities.J2
precision = quantities.precision
coordinates2elements, r_fit, t_fit, z_fit, vr_fit, vt_fit, vz_fit, GM, Rp, J2, $
    precision, a_fit, e_fit, I_fit, w_fit, M_fit, N_fit

;Extract longitudinally-averaged orbit elements a,e, from the fitting annulus, and
;assume longitude of periapse along streamline is linear in longitude t: w[t] = w0 + dwdt*t.
sz = size(r_fit)
Nt_fit = sz[1]
Nr_fit = sz[2]
a_avg = dblarr(Nr_fit)
e_avg = dblarr(Nr_fit)
w0_avg = dblarr(Nr_fit)
dwdt_avg = dblarr(Nr_fit)
for j = 0, Nr_fit - 1 do begin
    a_avg[j] = mean(a_fit[*, j])
    e_avg[j] = mean(e_fit[*, j])
    unwrap_angle, w_fit[*, j], w_fit_unw
    coeff = poly_fit(t_fit[*, j], w_fit_unw, 1)
    w0_avg[j] = mod2pi(coeff[0])
    dwdt_avg[j] = coeff[1]
endfor

;Semimajor axis spacing.
da = (a_avg[Nr_fit - 1] - a_avg[0])/(Nr_fit - 1)

;Semimajor axis of extrapolated annulus.
idx_extrap_in = rbc_indices[2]
idx_extrap_out = rbc_indices[3]
Nr_extrap = idx_extrap_out - idx_extrap_in + 1
Nt_extrap = Nt_fit
a_extrap = dblarr(Nt_extrap, Nr_extrap)
if (quantities.inner_rbc eq 'y') then $
    for j = 0, Nr_extrap - 1 do a_extrap[*, j] = a_avg[0] + da*(j - Nr_extrap)
if (quantities.outer_rbc eq 'y') then $
    for j = 0, Nr_extrap - 1 do a_extrap[*, j] = a_avg[Nr_fit - 1] + da*(j + 1)

;Eccentricity of extrapolated annulus.
coeff = poly_fit(a_avg, e_avg, 1, status=status)
e_extrap = dblarr(Nt_extrap, Nr_extrap)
for j = 0, Nr_extrap - 1 do e_extrap[0, j] = coeff[0] + coeff[1]*a_extrap[*, j]

;w0 in extrapolated annulus.
unwrap_angle, w0_avg, w0_avg_unwrapped
coeff = poly_fit(a_avg, w0_avg_unwrapped, 1, status=status)
w0_extrap = dblarr(Nt_extrap, Nr_extrap)
for j = 0, Nr_extrap - 1 do w0_extrap[0, j] = coeff[0] + coeff[1]*a_extrap[*, j]

;dwdt of extrapolated annulus
dwdt_extrap = dblarr(Nt_extrap, Nr_extrap)
for j = 0, Nr_extrap - 1 do dwdt_extrap[*, j] = mean(dwdt_avg)

;Longitude along each streamline in the extrapolated annulus
t_extrap = dblarr(Nt_extrap, Nr_extrap)
t_vector = !dpi*(2*dindgen(Nt_extrap)/Nt_extrap - 1)
for j = 0, Nr_extrap - 1 do t_extrap[0, j] = t_vector

;Longitude of periapse in the extrapolated annulus
w_extrap = mod2pi(w0_extrap + dwdt_extrap*t_extrap)

;Mean anomaly.
M_extrap = mod2pi(t_extrap - w_extrap)

;Vertical orbit elements are zero.
I_extrap = dblarr(Nt_extrap, Nr_extrap)
N_extrap = dblarr(Nt_extrap, Nr_extrap)

;Mass and ID of extrapolated region
mass_extrap = dblarr(Nt_extrap, Nr_extrap)
ID_extrap = intarr(Nt_extrap, Nr_extrap)

;Convert elements to coordinates and sort according to longitude.
elements2coordinates, a_extrap, e_extrap, I_extrap, w_extrap, M_extrap, N_extrap, GM, Rp, J2, $
    precision, r_extrap, t_extrap, z_extrap, vr_extrap, vt_extrap, vz_extrap
sort_ring_longitudes, r_extrap, t_extrap, z_extrap, vr_extrap, vt_extrap, vz_extrap, $
    mass_extrap, ID_extrap

;Insert extrapolated coordinates back into input arrays.
r[0, idx_extrap_in] = r_extrap
t[0, idx_extrap_in] = t_extrap
z[0, idx_extrap_in] = z_extrap
vr[0, idx_extrap_in] = vr_extrap
vt[0, idx_extrap_in] = vt_extrap
vz[0, idx_extrap_in] = vz_extrap

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro surface_density, r, t, lambda, difference_method, sigma
;
;  Calculate surface density at each ring particle's positions, ie across non-uniform grid.
;

;Surface density in the ring's interior.
interpolate_f, t, r,  1, difference_method, r_plus
interpolate_f, t, r, -1, difference_method, r_minus
sigma = 2*lambda/(r_plus - r_minus)

;Surface density along the ring's inner edge.
sigma[0, 0] = lambda[*, 0]/(r_plus[*, 0] - r[*, 0])

;Surface density along the ring's outer edge.
sz = size(r)
Nt = sz[1]
Nr = sz[2]
j = Nr - 1
sigma[0, j] = lambda[*, j]/(r[*, j] - r_minus[*, j])

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kick_coordinates, dt, r, s
;
;  Calculate the position kick on the disk particles and satellites that are
;  associated with the central planet's motion about the barycenter. Inputs
;  are the integrator's *full* timestep dt (noting that the kick that occurs here spans
;  timestep dt/2), and r=ring structure, s=satellite

;Convert mixed polar coordinates into mixed cartesian coordinates, for ring and satellite.
rt2xy, r.r, r.t, r.vr, r.vt, x_r, y_r, vx_r, vy_r
rt2xy, s.r, s.t, s.vr, s.vt, x_s, y_s, vx_s, vy_s

;Extract the z coordinates, velocities, and masses.
z_r = r.z
vz_r = r.vz
mass_r = r.mass
z_s = s.z
vz_s = s.vz
mass_s = s.mass

;Cartesian position kicks due to central planet's motion about barycenter.
m0 = 1d							;central planet's mass=1
factor = dt/2/m0
dx = factor*( total(mass_r*vx_r) +  total(mass_s*vx_s) )
dy = factor*( total(mass_r*vy_r) +  total(mass_s*vy_s) )
dz = factor*( total(mass_r*vz_r) +  total(mass_s*vz_s) )

;Kick the cartesian coordinates.
x_r = x_r + dx
y_r = y_r + dy
z_r = z_r + dz
x_s = x_s + dx
y_s = y_s + dy
z_s = z_s + dz

;Convert cartesian to polar coordinates, for ring and satellite.
xy2rt, x_r, y_r, vx_r, vy_r, r_r, t_r, vr_r, vt_r
xy2rt, x_s, y_s, vx_s, vy_s, r_s, t_s, vr_s, vt_s

;Store updated polar coordinates.
r.r = r_r
r.t = t_r
r.z = z_r
s.r = r_s
s.t = t_s
s.z = z_s

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kick_velocities, quantities, time, dt, ring, sat
;
;  Calculate the accelerations A on all ring particles and satellites, and kick
;  their velocities over timestep dt by amount dv=A*dt, and increment time by dt.
;

;Extract ring and satellite mixed coordinates and velocities, masses, and ID.
r_r = ring.r
t_r = ring.t
z_r = ring.z
vr_r = ring.vr
vt_r = ring.vt
vz_r = ring.vz
mass_r = ring.mass
ID_r = ring.ID
r_s = sat.r
t_s = sat.t
z_s = sat.z
vr_s = sat.vr
vt_s = sat.vt
vz_s = sat.vz
mass_s = sat.mass

;Sort ring particles according to their longitudes.
sort_ring_longitudes, r_r, t_r, z_r, vr_r, vt_r, vz_r, mass_r, ID_r

;Apply radiative boundary conditions, as needed.
if (quantities.inner_rbc eq 'y') then begin
    N_sl_rbc = quantities.N_sl_rbc
    rbc_indices = [N_sl_rbc, 2*N_sl_rbc - 1, 0, N_sl_rbc - 1]
    radiative_bc, rbc_indices, quantities, r_r, t_r, z_r, vr_r, vt_r, vz_r 
endif
if (quantities.outer_rbc eq 'y') then begin
    N_sl_rbc = quantities.N_sl_rbc
    Nrm1 = quantities.Nr - 1
    rbc_indices = [Nrm1 - 2*N_sl_rbc + 1, Nrm1 - N_sl_rbc, Nrm1 - N_sl_rbc + 1, Nrm1]
    radiative_bc, rbc_indices, quantities, r_r, t_r, z_r, vr_r, vt_r, vz_r 
endif

;Calculate accelerations on ring and satellites.
accelerations, quantities, time, mass_r, r_r, t_r, z_r, vr_r, vt_r, vz_r, $
    mass_s, r_s, t_s, z_s, Ar_r, At_r, Az_r, Ar_s, At_s, Az_s

;Update ring coordinates and mass and ID, and kick their velocities.
if (quantities.Np gt 0) then begin
    ring.r = r_r
    ring.t = t_r
    ring.z = z_r
    ring.mass = mass_r
    ring.ID = ID_r
    ring.vr = vr_r + Ar_r*dt
    ring.vt = vt_r + At_r*dt
    ring.vz = vz_r + Az_r*dt
endif

;Kick satellite velocities and update their mass in case they are growing.
if (quantities.Ns gt 0) then begin 
    sat.vr = sat.vr + Ar_s*dt
    sat.vt = sat.vt + At_s*dt
    sat.vz = sat.vz + Az_s*dt
    sat.mass = mass_s
endif

;advance time
time = time + dt

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro accelerations, quantities, time, mass_r, r_r, t_r, z_r, vr_r, vt_r, vz_r, $
    mass_s, r_s, t_s, z_s, Ar_r, At_r, Az_r, Ar_s, At_s, Az_s, $
    Ar_rs=Ar_rs, At_rs=At_rs, Az_rs=Az_rs, Ar_ss=Ar_ss, At_ss=At_ss, Az_ss=Az_ss, $
    Ar_rr=Ar_rr, At_rr=At_rr, Az_rr=Az_rr, Ar_rv=Ar_rv, At_rv=At_rv, Az_rv=Az_rv, $
    Ar_rp=Ar_rp, At_rp=At_rp, Az_rp=Az_rp, Ar_rd=Ar_rd, At_rd=At_rd, Az_rd=Az_rd, $
    sigma_r=sigma_r, shear_viscosity=shear_viscosity, bulk_viscosity=bulk_viscosity, c=c
;
;  Calculate accelerations on ring and satellites, where Ar_r=radial accelerations of ring,
;  At_s=tangential acceleration of satellites, etc. Optional parameters provide accelerations
;  due to specific phenomen, such as At_rs=tangential acceleration of ring due to satellites,
;

;Storage for accelerations on ring and satellites.
Nt = quantities.Nt
Nr = quantities.Nr
Ns = quantities.Ns
Ar_r = dblarr(Nt, Nr)
At_r = dblarr(Nt, Nr)
Az_r = dblarr(Nt, Nr)
Ar_s = dblarr(Ns)
At_s = dblarr(Ns)
Az_s = dblarr(Ns)

;The ring's linear mass density to zero order in eccentricity.
Rp = quantities.Rp
J2 = quantities.J2
GM = 1d
semimajor_axis, r_r, z_r, vr_r, vt_r, vz_r, GM, Rp, J2, a_r
lambda_r = (Nt/!dpi/2d)*mass_r/a_r

;Ring surface density, as needed.
viscosity = quantities.viscosity
pressure = quantities.pressure
if ((viscosity eq 'y') or (pressure eq 'y')) then $
    surface_density, r_r, t_r, lambda_r, quantities.difference_method, sigma_r

;Grow satellite masses, as needed.
t_grow_sat = quantities.t_grow_sat
mass_s_final = quantities.mass_s_final
if (t_grow_sat gt 0) then mass_s = mass_s_final*( 1 - exp(-time/t_grow_sat) )

;Calculate acceleration of ring due to satellites' gravity, plus acceleration of
;each satellite due to ring gravity, plus satellite-satellite gravities, as needed,
;where Ar_rs=radial acceleration of ring due to satellites.
ring_gravity = quantities.ring_gravity
if (quantities.satellite_gravity eq 'y') then begin
    sat_grav, quantities, time, mass_r, r_r, t_r, z_r, mass_s, r_s, t_s, z_s, $
        Ar_rs, At_rs, Az_rs, Ar_sr, At_sr, Az_sr, Ar_ss, At_ss, Az_ss
    Ar_r = Ar_r + Ar_rs
    At_r = At_r + At_rs
    Az_r = Az_r + Az_rs
    Ar_s = Ar_s + Ar_ss
    At_s = At_s + At_ss
    Az_s = Az_s + Az_ss
    if (ring_gravity eq 'y') then begin
        Ar_s = Ar_s + Ar_sr
        At_s = At_s + At_sr
        Az_s = Az_s + Az_sr
    endif
endif

;Calculate the acceleration of ring particles due to ring gravity, where 
;At_rr = tangential acceleration of ring due to ring self gravity, etc.
if (ring_gravity eq 'y') then begin
    ring_gravity, quantities, lambda_r, r_r, t_r, vr_r, vt_r, Ar_rr, At_rr, Az_rr
    Ar_r = Ar_r + Ar_rr
    At_r = At_r + At_rr
    Az_r = Az_r + Az_rr
endif

;Acceleration of ring particles due to drag, as needed.
if (quantities.drag_force eq 'y') then begin
    drag, quantities, time, mass_r, r_r, t_r, z_r, vr_r, vt_r, vz_r, Ar_rd, At_rd, Az_rd
    Ar_r = Ar_r + Ar_rd
    At_r = At_r + At_rd
    Az_r = Az_r + Az_rd
endif

;Calculate accelerations of ring particles due to viscosity.
if (viscosity eq 'y') then begin
    viscosity, quantities, lambda_r, sigma_r, r_r, t_r, z_r, vr_r, vt_r, vz_r, $
        Ar_rv, At_rv, Az_rv, c=c, shear_viscosity=shear_viscosity, $
        bulk_viscosity=bulk_viscosity
    Ar_r = Ar_r + Ar_rv
    At_r = At_r + At_rv
    Az_r = Az_r + Az_rv
endif

;Calculate accelerations of ring particles due to pressure.
if (pressure eq 'y') then begin
    pressure, quantities, lambda_r, sigma_r, r_r, t_r, z_r, vr_r, vt_r, vz_r, c, $
        Ar_rp, At_rp, Az_rp
    Ar_r = Ar_r + Ar_rp
    At_r = At_r + At_rp
    Az_r = Az_r + Az_rp
endif

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sat_grav, quantities, time, mass_r, r_r, t_r, z_r, mass_s, r_s, t_s, z_s, $
    Ar_rs, At_rs, Az_rs, Ar_sr, At_sr, Az_sr, Ar_ss, At_ss, Az_ss
;
;  Calculate accelerations of ring Ar_rs,At_rs,Az_rs that are due to the satellites,
;  and the acceleration of each satellite Ar_sr,At_sr,Az_sr due to ring gravity,
;  and accelerations Ar_ss,At_ss,Az_ss due to satellite interactions.

;Sin and cos of ring & satellite longitudes, and their cartesian coordinates.
sin_t_r = sin(t_r)
cos_t_r = cos(t_r)
x_r = r_r*cos_t_r
y_r = r_r*sin_t_r
sin_t_s = sin(t_s)
cos_t_s = cos(t_s)
x_s = r_s*cos_t_s
y_s = r_s*sin_t_s

;loop over each non-fourier satellite and calculate the cartesian acceleration that satellite
;exerts on the ring Ax_rs,Ay_rs,Az_rs, the acceleration the ring exerts on the
;satellite Ax_sr,Ay_sr,Az_sr, and the satellite-satellite accelerations Ax_ss,Ay_ss,Az_ss,
Ax_rs = 0d
Ay_rs = 0d
Az_rs = 0d
Ns = quantities.Ns
Ax_sr = dblarr(Ns)
Ay_sr = dblarr(Ns)
Az_sr = dblarr(Ns)
Ax_ss = dblarr(Ns)
Ay_ss = dblarr(Ns)
Az_ss = dblarr(Ns)
sat_index = quantities.sat_index
fourier_satellite = quantities.fourier_satellite
t_grow_sat = quantities.t_grow_sat
mass_s_final = quantities.mass_s_final
for j = 0, Ns - 1 do begin

    ;Proceed if this satellite is not a fourier satellite.
    if ((fourier_satellite ne 'y') or (j ne sat_index)) then begin

        ;ring-satellite separation squared and cubed.
        dz = z_s[j] - z_r
        delta_2 = r_r*r_r + r_s[j]*r_s[j] - 2*r_r*r_s[j]*cos(t_s[j] - t_r) + dz^2
        delta_3 = delta_2*sqrt(delta_2)

        ;acceleration of ring due to satellite j, in cartesian coordinates, assuming G=1
        factor = mass_s[j]/delta_3
        Ax_rs_j = factor*(x_s[j] - x_r)
        Ay_rs_j = factor*(y_s[j] - y_r)
        Az_rs_j = factor*dz
        Ax_rs = Ax_rs + Ax_rs_j
        Ay_rs = Ay_rs + Ay_rs_j
        Az_rs = Az_rs + Az_rs_j

        ;acceleration of satellite due to ring, in cartesian coordinates
        if (mass_s[j] gt 0) then begin 
            factor = -mass_r/mass_s[j]
            Ax_sr[j] = total(factor*Ax_rs_j)
            Ay_sr[j] = total(factor*Ay_rs_j)
            Az_sr[j] = total(factor*Az_rs_j)
        endif

        ;gravity on satellite j due to other satellites k
        for k = 0, Ns - 1 do begin
            if (k ne j) then begin
                dz = z_s[k] - z_s[j]
                delta_2 = r_s[j]*r_s[j] + r_s[k]*r_s[k] - $
                    2*r_s[j]*r_s[k]*cos(t_s[k] - t_s[j]) + dz*dz
                delta_3 = delta_2*sqrt(delta_2)
                factor = mass_s[k]/delta_3
                Ax_ss[j] = Ax_ss[j] + factor*(x_s[k] - x_s[j])
                Ay_ss[j] = Ay_ss[j] + factor*(y_s[k] - y_s[j])
                Az_ss[j] = Az_ss[j] + factor*dz
            endif
        endfor

    endif

endfor

;Acceleration of ring and satellites in polar coordinates.
Ar_rs =  Ax_rs*cos_t_r + Ay_rs*sin_t_r
At_rs = -Ax_rs*sin_t_r + Ay_rs*cos_t_r
Ar_sr =  Ax_sr*cos_t_s + Ay_sr*sin_t_s
At_sr = -Ax_sr*sin_t_s + Ay_sr*cos_t_s
Ar_ss =  Ax_ss*cos_t_s + Ay_ss*sin_t_s
At_ss = -Ax_ss*sin_t_s + Ay_ss*cos_t_s

;Alternatively, only the m^th fourier component of the satellite's gravity perturbs the ring.
;This assumes that the ring semimajor axes have drifted negligibly since the run started.
if (fourier_satellite eq 'y') then begin 
    m_LR = quantities.m_LR
    beta = quantities.beta
    lc = quantities.lc
    dlc = quantities.dlc
    fn = dlc
    if (m_LR eq 1) then fn = fn - 1
    angle = m_LR*(t_r - t_s[sat_index])
    if (t_grow_sat lt 0) then mass_s[sat_index] = mass_s_final[sat_index]
    factor = mass_s[sat_index]/(r_s[sat_index]*r_s[sat_index])
    Ar_rs = Ar_rs + factor*fn*cos(angle)
    fn = lc
    if (m_LR eq 1) then fn = fn - 1
    At_rs = At_rs - factor*fn*sin(angle)*m_LR/beta
    mass_s[sat_index] = 0
endif

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ring_gravity, quantities, lambda, r, t, vr, vt, Ar, At, Az
;
;  calculate acceleration Ar,At,Az due to ring gravity. Note that this code
;  assumes z=0 so Az=0.
;

;Storage for ring's accelerations.
sz = size(r)
Nt = sz[1]
Nr = sz[2]
Ar = dblarr(Nt, Nr)
At = dblarr(Nt, Nr)
Az = dblarr(Nt, Nr)

;Loop over all streamlines and calculate ring gravity assuming G=1.
;Calculate streamline separations, the gravity that one
;streamline exerts on another, then the r,t components of that gravity. Note that
;the current code ignores the ring's vertical displacements!
difference_method = quantities.difference_method
for j = 1, Nr - 1 do begin
    interpolate_f, t, r, j, difference_method, r_j
    g = 2*shift(lambda, 0, -j)/(r_j - r)
    Ar = Ar + g
    At = At - g*shift(vr/vt, 0, -j)				;approximate
    ;interpolate_f, t, vr/vt, j, difference_method, vr_vt_j	;exact
    ;At = At - g*vr_vt_j
endfor

;Acceleration due to adjacent particles in streamline.
longitude_wrap, t, r, tw, rw, in_array=lambda, wrapped_array=lambda_w
rtw = rw*tw
delta_plus = shift(rtw, -1, 0) - rtw
delta_minus = rtw - shift(rtw, 1, 0)
lambda_plus = shift(lambda_w, -1, 0)
lambda_minus = shift(lambda_w, 1, 0)
At_w_sl = 2*(lambda_plus/delta_plus - lambda_minus/delta_minus)
At_sl = At_w_sl[1:Nt, *]
At = At + At_sl
Ar = Ar +  At_sl*vr/vt

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro viscosity, quantities, lambda, sigma, r, t, z, vr, vt, vz, Ar, At, Az, $
    c=c, shear_viscosity=shear_viscosity, bulk_viscosity=bulk_viscosity
;
;  Acceleration due to viscosity. The choices of viscosity_law='constant',
;  'collisions', and 'fake c' are experimental and not thoroughly tests.
;

;Number of particles in each streamline Nt and number of streamlines Nr.
sz = size(r)
Nt = sz[1]
Nr = sz[2]

;Ring angular velocity vt/r.
Omega = vt/r

;Quantities used below.
difference_method = quantities.difference_method
derivative_dfdr, r, t, Omega, difference_method, dOmega_dr
rSigma = r*sigma
rdOmega_dr = r*dOmega_dr
c = quantities.c0

;Get ring viscosity from viscosity law.
viscosity_law = quantities.viscosity_law
if (viscosity_law eq 'constant') then begin
    shear_viscosity = quantities.shear_viscosity
    bulk_viscosity = quantities.bulk_viscosity
endif
if (viscosity_law eq 'collisions') then begin
    tau = quantities.opacity*sigma
    tau0 = quantities.opacity*quantities.surface_density
    c0 = quantities.c0
    factor = 1.5d*Omega/rdOmega_dr
    c = (tau0*c0)*factor*factor/tau
    shear_viscosity = (0.46d)*tau*c*c/Omega/(1 + tau*tau)
    bulk_viscosity = shear_viscosity
endif
if (viscosity_law eq 'fake c') then begin
    x0 = 1d-5
    cfactor = 100d
    c0 = quantities.c0
    tau = quantities.opacity*sigma
    x = r - quantities.r_ref
    c = c0*cfactor*x0/sqrt(x^2 + x0^2)
    shear_viscosity = (0.46d)*tau*c*c/Omega/(1 + tau*tau)
    bulk_viscosity = shear_viscosity
endif
if (viscosity_law eq 'wakes') then begin
    sigma_ratio_sq = (sigma/quantities.surface_density)^2
    shear_viscosity = quantities.shear_viscosity*sigma_ratio_sq
    bulk_viscosity = quantities.bulk_viscosity*sigma_ratio_sq
endif

;Viscous angular momentum flux F, and radial gradient dF_dr.
F = (-shear_viscosity)*rSigma*rdOmega_dr
derivative_dfdr, r, t, F, difference_method, dF_dr

;Linear momentum flux G due to viscosity, and radial gradient dG_dr.
derivative_dfdr, r, t, vr, difference_method, dvr_dr
f1 = -(4*shear_viscosity/3 + bulk_viscosity)
f2 = -(bulk_viscosity - 2*shear_viscosity/3)
G = (f1*dvr_dr + f2*vr/r)*sigma
derivative_dfdr, r, t, G, difference_method, dG_dr

;Radial and tangential accelerations in the ring's interior.
Ar = -dG_dr/sigma
At = -dF_dr/rSigma

;Acceleration at the disk's sharp inner edge, which is zero if the edge is confined.
j = 0
if (quantities.confine_inner_edge eq 'y') then begin 
    Ar[*, j] = 0d
    At[*, j] = 0d
endif else begin 
    Ar[0, j] = -G[*, j]/lambda[*, j]
    At[0, j] = -F[*, j]/r[*, j]/lambda[*, j]
endelse

;Acceleration at the disk's sharp outer edge, which is zero if the edge is confined.
j = Nr - 1
if (quantities.confine_outer_edge eq 'y') then begin
    Ar[*, j] = 0d
    At[*, j] = 0d
endif else begin 
    interpolate_f, t, G, -1, difference_method, G_minus
    interpolate_f, t, F, -1, difference_method, F_minus
    Ar[0, j] = G_minus[*, j]/lambda[*, j]
    At[0, j] = F_minus[*, j]/r[*, j]/lambda[*, j]
endelse

;Vertical accelerations are zero for now.
Az = 0d

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro pressure, quantities, lambda, sigma, r, t, z, vr, vt, vz, c, Ar, At, Az
;
;  acceleration due to pressure.
;

;Number of particles in each streamline Nt and number of streamlines Nr.
sz = size(r)
Nt = sz[1]
Nr = sz[2]

;Pressure.
if (keyword_set(c) eq 0) then c = quantities.c0
P = c*c*sigma

;Radial acceleration in the ring's interior.
difference_method = quantities.difference_method
derivative_dfdr, r, t, P, difference_method, dP_dr
Ar = -dP_dr/sigma

;Acceleration at the ring's sharp inner edge.
j = 0
Ar[0, j] = -P[*, j]/lambda[*, j]

;Acceleration at the ring's sharp outer edge.
interpolate_f, t, P, -1, difference_method, P_minus
j = Nr - 1
Ar[0, j] = P_minus[*, j]/lambda[*, j]

;Tangential acceleration.
At = -Ar*vr/vt

;Additional acceleration due to pressure gradient along streamline.
longitude_wrap, t, r, tw, rw
longitude_wrap, t, P, tw, Pw
longitude_wrap, t, sigma, tw, sigma_w
delta_Pw = shift(Pw, -1, 0) - shift(Pw, 1, 0)
delta_tw = shift(tw, -1, 0) - shift(tw, 1, 0)
Ap_sl_w = -delta_Pw/(delta_tw*rw*sigma_w)
Ap_sl = Ap_sl_w[1:Nt, *]
At = At + Ap_sl
Ar = Ar + Ap_sl*vr/vt

;Vertical accelerations are zero for now.
Az = 0d

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro drag, quantities, time, mass, r, t, z, vr, vt, vz, Ar, At, Az
;
;  The radial Ar and tangential At acceleration due to drag. Inputs are the
;  current time and polar positions and velocities. 
;

;Calculate semimajor axis a and angular velocity Omega(a).
GM = 1 + mass
Rp = quantities.Rp
J2 = quantities.J2
semimajor_axis, r, z, vr, vt, vz, GM, Rp, J2, a
orbit_frequencies, GM, Rp, J2, a, Omega, Kappa, Eta, Beta, Nu

;Adjust drag coefficient so it decays exponentially with time when drag_timescale>0.
drag_factor = quantities.drag_coeff
drag_timescale = quantities.drag_timescale
t_start = quantities.t_start
if (drag_timescale gt 0) then drag_factor = drag_factor*exp(-(time - t_start)/drag_timescale)

;Radial, tangential, and vertical accelerations due to drag.
factor = -drag_factor*Omega
Ar = factor*vr
At = factor*(vt - a*Omega)
Az = factor*vz

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro drift, quantities, dt, p
;
;  Particles experience unperturberd epicyclic orbital motion during timestep dt.
;  Structure p contains the particles cylindrical coordinates and velocities.
;

;Get oblateness parameters
Rp = quantities.Rp
J2 = quantities.J2
precision = quantities.precision

;Convert coordinates to elements assuming GM=1, and get orbit frequencies.
GM = 1d
coordinates2elements, p.r, p.t, p.z, p.vr, p.vt, p.vz, GM, Rp, J2, precision, $
    a, e, I, w, M, N, Omega=Omega_0, Kappa=Kappa_0, Eta=Eta_0, Beta=Beta_0, Nu=Nu_0

;Include the e^2 modifications to these frequences, from Eqns (50-52) of B-R&L(1994).
;The I^2 modifications are not done, nor is Nu modified.
if (precision eq 1) then begin 
    Omega = Omega_0
    Kappa = Kappa_0
    Nu = Nu_0
endif
if (precision eq 2) then begin 
    EK = Eta_0/Kappa_0
    EK2 = EK*EK
    EK4 = EK2*EK2
    OK = Omega_0/Kappa_0
    OK2 = OK*OK
    BK = Beta_0/Kappa_0
    BK2 = BK*BK
    e2 = e*e
    Omega = Omega_0*( 1 + (1.5d - EK2*3d)*e2 )
    Kappa = Kappa_0*( 1 + ( (OK2-EK4)*3.75d - BK2*1.5d)*e2 )
    Nu = Nu_0
endif

;Advance longitude of periapse, longitude of ascending node, and mean anomaly.
w = mod2pi( w + (Omega - Kappa)*dt )
N = mod2pi( N + (Omega - Nu)*dt )
M = mod2pi( M + Kappa*dt )

;Convert orbit elements back to coordinates.
elements2coordinates, a, e, I, w, M, N, GM, Rp, J2, precision, r, t, z, vr, vt, vz

;Store updated coordinates and velocities in structure p.
p.r = r
p.t = t
p.z = z
p.vr = vr
p.vt = vt
p.vz = vz

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro store_results, tm, t_now, time, quantities, r, s, c
;
;  r=ring structure, s=satellite structure, c=coordinate array of structures.
;

;Update time.
time[tm] = t_now

;Extract ring and satellites' planetocentric coordinates, barycentric velocities, mass, and ID.
r_r = r.r
t_r = r.t
z_r = r.z
vr_r_BC = r.vr
vt_r_BC = r.vt
vz_r_BC = r.vz
mass_r = r.mass
ID_r = r.ID
r_s = s.r
t_s = s.t
z_s = s.z
vr_s_BC = s.vr
vt_s_BC = s.vt
vz_s_BC = s.vz
mass_s = s.mass

;Convert barycentric velocities to planetocentric velocities.
mixed_rv2planeto_v, mass_r, r_r, t_r, z_r, vr_r_BC, vt_r_BC, vz_r_BC, mass_s, r_s, t_s, z_s, $
    vr_s_BC, vt_s_BC, vz_s_BC, vr_r, vt_r, vz_r, vr_s, vt_s, vz_s

;Sort ring particles according to their longitudes.
sort_ring_longitudes, r_r, t_r, z_r, vr_r, vt_r, vz_r, mass_r, ID_r

;Create a structure containing all current planetocentric coordinates, masses, and ring ID.
c[tm] = create_struct('r_r', r_r, 't_r', t_r, 'z_r', z_r, 'vr_r', vr_r, 'vt_r', vt_r, $
    'vz_r', vz_r, 'mass_r', mass_r, 'ID_r', ID_r, 'r_s', r_s, 't_s', t_s, 'z_s', z_s, $
    'vr_s', vr_s, 'vt_s', vt_s, 'vz_s', vz_s, 'mass_s', mass_s)

;convert planetocentric polar coordinates to planetocentric cartesian coordinates
rt2xy, r_r, t_r, vr_r, vt_r, x_r, y_r, vx_r, vy_r
rt2xy, r_s, t_s, vr_s, vt_s, x_s, y_s, vx_s, vy_s

;convert planetocentric cartesian coordinates and velocities to barycentric
;cartesian coordinates and velocities.
planeto2bary_xyz, mass_r, x_r, y_r, z_r, mass_s, x_s, y_s, z_s, $
    x_r_BC, y_r_BC, z_r_BC, x_s_BC, y_s_BC, z_s_BC, x_0, y_0, z_0 
planeto2bary_xyz, mass_r, vx_r, vy_r, vz_r, mass_s, vx_s, vy_s, vz_s, $
    vx_r_BC, vy_r_BC, vz_r_BC, vx_s_BC, vy_s_BC, vz_s_BC, vx_0, vy_0, vz_0 

;system's total angular z-component of momentum assuming central mass M_0=1
M_0 = 1d
Lz = total( mass_r*(x_r_BC*vy_r_BC - y_r_BC*vx_r_BC) ) + $
     total( mass_s*(x_s_BC*vy_s_BC - y_s_BC*vx_s_BC) ) +  M_0*(x_0*vy_0 - y_0*vx_0)
quantities.ang_mom[tm] = Lz

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro extract_results, quantities, coordinates, ID_r, mass_r, r_r, t_r, z_r, $
    vr_r, vt_r, vz_r, a_r, e_r, I_r, w_r, M_r, N_r, mass_s, r_s, t_s, z_s, $
    vr_s, vt_s, vz_s, a_s, e_s, I_s, w_s, M_s, N_s, sort_ID=sort_ID
;
;
;

;Get array sizes
Nt = quantities.Nr
Nr = quantities.Nr
Noutput = quantities.Noutput

;Extract ring's mass, ID, and planetocentric coordinates and velocities.
mass_r = coordinates.mass_r
ID_r = coordinates.ID_r
r_r = coordinates.r_r
t_r = coordinates.t_r
z_r = coordinates.z_r
vr_r = coordinates.vr_r
vt_r = coordinates.vt_r
vz_r = coordinates.vz_r

;Sort ring particles along each streamline according to their ID number,
;rather than by their longitude, if desired.
if (keyword_set(sort_id) eq 1) then begin
    for tm = 0, Noutput - 1 do begin
        for j = 0, Nr - 1 do begin
            s = sort(ID_r[*, j, tm])
            ID_r[0, j, tm] = ID_r[s, j, tm]
            mass_r[0, j, tm] = mass_r[s, j, tm]
            r_r[0, j, tm] = r_r[s, j, tm]
            t_r[0, j, tm] = t_r[s, j, tm]
            z_r[0, j, tm] = z_r[s, j, tm]
            vr_r[0, j, tm] = vr_r[s, j, tm]
            vt_r[0, j, tm] = vt_r[s, j, tm]
            vz_r[0, j, tm] = vz_r[s, j, tm]
        endfor 
    endfor
endif

;Extract satellites mass and planetocentric coordinates and velocities.
mass_s = coordinates.mass_s
r_s = coordinates.r_s
t_s = coordinates.t_s
z_s = coordinates.z_s
vr_s = coordinates.vr_s
vt_s = coordinates.vt_s
vz_s = coordinates.vz_s

;Convert coordinates to orbit elements.
J2 = quantities.J2
Rp = quantities.Rp
precision = quantities.precision
GM_r = 1 + mass_r
coordinates2elements, r_r, t_r, z_r, vr_r, vt_r, vz_r, GM_r, Rp, J2, precision, $
    a_r, e_r, I_r, w_r, M_r, N_r
GM_s = 1 + mass_s
coordinates2elements, r_s, t_s, z_s, vr_s, vt_s, vz_s, GM_s, Rp, J2, precision, $
    a_s, e_s, I_s, w_s, M_s, N_s

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro surface_density_grid, lambda, quantities, r, theta, Nr_map, Nt_map, r_map, $
    theta_map, sd_map
;
;  
;

;extract parameters
Nr = quantities.Nr
Nt = quantities.Nt
r_range = quantities.r_range

;This procedure requires at least Nr=3 streamlines to work. If not, then zero
;arrays are returned.
if (Nr le 3) then begin
    r_map = dblarr(Nt_map, Nr_map)
    theta_map = dblarr(Nt_map, Nr_map)
    sd_map = dblarr(Nt_map, Nr_map)
    return
endif

;wrap r, theta around in longitude
;Nw = 1
longitude_wrap, theta, r, theta_w, r_w

;resample r(theta) along uniform longitude grid
theta_u = !dpi*(2*dindgen(Nt_map)/Nt_map - 1)
r_u = dblarr(Nt_map, Nr)
t_u = dblarr(Nt_map, Nr)
for j = 0, Nr - 1 do begin $
    r_u[0, j] = spline(theta_w[*, j], r_w[*, j], theta_u) & $
    t_u[0, j] = theta_u & $
endfor

;resampled lambda array
lambda_u = dblarr(Nt_map, Nr)
for j = 0, Nr - 1 do lambda_u[*, j] = lambda[0, j]

;calculate surface density over a uniformly sampled longitudinal grid and an
;irregularly sampled radial grid
difference_method = 'polynomial'
surface_density, r_u, t_u, lambda_u, difference_method, sigma_u

;reference radius
r1 = 1

;map's inner and outer radii
r_min = r1 + min(r_range)
r_max = r1 + max(r_range)

;uniform radial grid where surface density will be sampled
r_v = (r_max - r_min)*dindgen(Nr_map)/(Nr_map - 1) + r_min

;use spline to interpolate surface density across uniform radial grid
r_map = dblarr(Nt_map, Nr_map)
theta_map = dblarr(Nt_map, Nr_map)
sd_map = dblarr(Nt_map, Nr_map)
for j = 0, Nt_map - 1 do begin $
    r_map[j, 0] = transpose(r_v) & $
    theta_map[j, *] = theta_u[j] & $
    r0 = reform(r_u[j,*]) & $
    sd0 = reform(sigma_u[j,*]) & $
    k = where( (r_v ge min(r0)) and (r_v le max(r0)), Nk ) & $
    if (Nk gt 0) then begin $
        sd_spline = spline(r0, sd0, r_v[k]) & $
        sd_map[j, k] = transpose(sd_spline) & $
    endif & $
endfor

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro display_results, inputs, quantities, tm, time, coordinates, $
    sd_byte=sd_byte, t_w=t_w, a_w=a_w, dr=dr
;
;  Display results.
;

;Setup for showing plots on monitor or else in a postscript file.
display = quantities.display
if (display eq 'monitor') then begin
    erase
    chrsize = 1.3
    thck = 1
endif
if (display eq 'postscript') then begin
    set_plot,'ps'
    device, xsize=20.0, ysize=18.2, xoffset=1, yoffset=4, /color, /portrait, $
        bits_per_pixel=8
    chrsize=0.8
    thck = 2
endif

;These coordinates are used to position various plots.
offset = [0.09, 0.05, -0.02, -0.03]
x0 = 0.000
x1 = 0.333
x2 = 0.667
x3 = 1.000
y0 = 0.000
y1 = 0.333
y2 = 0.667
y3 = 1.000
pos00 = [ x0, y0, x1, y1 ] + offset		;lower left figure
pos10 = [ x1, y0, x2, y1 ] + offset		;lower right
pos01 = [ x0, y1, x1, y2 ] + offset		;middle left
pos11 = [ x1, y1, x2, y2 ] + offset		;middle right
pos02 = [ x0, y2, x1, y3 ] + offset		;upper left
pos12 = [ x1, y2, x2, y3 ] + offset		;upper right
pos20 = [ x2, y0, x3, y3 ] + offset		;rightmost column

;Extract the ring and satellites' planetocentric coordinates and velocities, masses, 
;ring ID, and number of bodies.
r_r = coordinates[tm].r_r
t_r = coordinates[tm].t_r
z_r = coordinates[tm].z_r
vr_r = coordinates[tm].vr_r
vt_r = coordinates[tm].vt_r
vz_r = coordinates[tm].vz_r
mass_r = coordinates[tm].mass_r
ID_r = coordinates[tm].ID_r
r_s = coordinates[tm].r_s
t_s = coordinates[tm].t_s
z_s = coordinates[tm].z_s
vr_s = coordinates[tm].vr_s
vt_s = coordinates[tm].vt_s
vz_s = coordinates[tm].vz_s
mass_s = coordinates[tm].mass_s
Nr = quantities.Nr
Nt = quantities.Nt

;Rotate longitudes to reference frame that corotates with satellite's longitude, as needed.
sat_index = inputs.sat_index
corotate_sat = inputs.corotate_sat
if (corotate_sat eq 'y') then begin 
    t_r = mod2pi(t_r - t_s[sat_index])
    t_s = mod2pi(t_s - t_s[sat_index])
endif

;Convert ring's planetocentric polar coordinates into orbit elements assuming GM=1.
GM = 1d
Rp = quantities.Rp
J2 = quantities.J2
precision = quantities.precision
coordinates2elements, r_r, t_r, z_r, vr_r, vt_r, vz_r, GM, Rp, J2, precision, $
    a_r, e_r, I_r, w_r, M_r, N_r

;Rotate all longitudes to reference frame that corotates with ring's middle streamline's
;longitude of periapse, as needed.
corotate_ring = inputs.corotate_ring
if (corotate_ring eq 'y') then begin 
    streamline_index = Nr/2 & $
    t_ref = median(w_r[*, streamline_index]) & $
    t_r = mod2pi(t_r - t_ref)
    w_r = mod2pi(w_r - t_ref)
    t_s = mod2pi(t_s - t_ref)
endif

;Sort ring particles according to their longitudes.
sort_ring_longitudes, r_r, t_r, z_r, vr_r, vt_r, vz_r, mass_r, ID_r

;Calculate ring's linear surface density, with a fix for when ring is massless
lambda_r = (Nt/!dpi/2d)*mass_r/a_r
if (total(lambda_r) eq 0) then begin 
    surface_density = 1
    dr = shift(r_r, 0, -1) - r_r
    dr[0, Nr - 1] = r_r[*, Nr - 1] - r_r[*, Nr - 2]
    lambda_r = surface_density*dr
endif

;Generate surface density map relative to disk's undisturbed surface density,
;and radial distance relative to r_rel.
N_pixels = 256
surface_density_grid, lambda_r, quantities, r_r, t_r, N_pixels, N_pixels, $
    r_map, t_map, sd_map
surface_density = quantities.surface_density
r_ref = inputs.r_ref
if (surface_density gt 0) then sd_map = sd_map/surface_density
dr_map = r_map - r_ref

;Convert surface density into bytes, and display.
sd_range = inputs.sd_range
sd_min = inputs.sd_range[0]
sd_max = inputs.sd_range[1]
sd_byte = (sd_map > sd_min) < sd_max
sd_byte = sd_byte - sd_min
if (max(sd_byte) gt 0) then begin
    sd_byte = sd_byte/(sd_max - sd_min)
    sd_byte = fix(255*sd_byte)
    if (display eq 'postscript') then begin
        tv, sd_byte, 0.06, 0.7, XSIZE=0.26, YSIZE=0.275, /NORMAL
    endif else begin
        tv, sd_byte, 0.044, 0.693, /normal
    endelse
    xyouts, 0.14, 0.684, 'longitude !7h!3 !9 6!3', /NORMAL, charsize=chrsize, charthick=thck 
    xyouts, 0.034, 0.80, 'radius r !9 6!3', /NORMAL, charsize=chrsize, charthick=thck, $
        orientat=90
endif

;Create string containing the current time.
time_now = time[tm]
ttl = 'time t = 0.0'
if (time_now gt 0) then begin 
    str_len = round( alog10(time_now) + 1 )
    ttl = 'time t = ' + strmid(strtrim(string(time_now), 2), 0, str_len)
endif

;Define psym=8 as filled dots.
Nd = 41
phi_d = 2*!pi*findgen(Nd)/(Nd - 1)
xd = cos(phi_d)
yd = sin(phi_d)
usersym, xd, yd, fill=1

;Plot radial profile of surface density.
r_range = inputs.r_range
plot, dr_map[0, *], sd_map[0, *], ystyle=1, xstyle=1, charsize=chrsize, psym=8, $
    xtitle='radial distance r - r!lREF!n', noerase=1, ytitle='!7r!3(r)/!7r!3!lo!n', $
    position=pos12, yrange=sd_range, title=ttl, xticks=2, symsize=0.4, nodata=1, $
    thick=thck, xthick=thck, ythick=thck, charthick=thck, xrange=r_range
sd_longitude = inputs.sd_longitude
for j = 0, n_elements(sd_longitude) - 1 do begin $
    idx = where(t_map[*, 0] ge sd_longitude[j], Nidx) & $
    if (Nidx gt 0) then begin $
        idx = idx[0] & $
        l = where(sd_map[idx, *] gt 0, Nl) & $
        if (Nl gt 0) then oplot, dr_map[idx, l], sd_map[idx, l], thick=thck & $
    endif & $
endfor
plots, [0, 0], sd_range, linestyle=1, color=128, thick=thck

;Plot tangential cut of surface density, following streamline
ttl = '!7r!3 along middle streamline'
plot, t_map[*, 0]/!dpi, sd_map[*, 0], ystyle=1, xstyle=1, charsize=chrsize, $
    xtitle='!7h/p!3    (rad)', ytitle='!7r!3(!7h!3)/!7r!3!lo!n', $
    position=pos01, yrange=sd_range, xrange=[-1, 1], $
    symsize=0.4, nodata=1, noerase=1, title=ttl, $
    thick=thck, xthick=thck, ythick=thck, charthick=thck 
idx = Nr/2
if (Nt ge 3) then begin 
    r_spline = spline(t_r[*, idx], r_r[*, idx], t_map[*, 0])
    sigma_plot = dblarr(N_pixels)
    for j = 0, N_pixels - 1 do begin
        k = where(r_map[j, *] ge r_spline[j])
        k = k[0]
        if (k ge 0) then sigma_plot[j] = sd_map[j, k]
    endfor
    oplot, t_map[*, 0]/!dpi, sigma_plot, thick=thck
endif
oplot, [0, 0], sd_range, color=128, linestyle=1, thick=thck

;Plot e versus a - r_ref.
idx = Nt/2
ecc_range = inputs.ecc_range
da = a_r - r_ref
e_log = 0
if (inputs.ecc_axis eq 'log') then e_log = 1
plot, da[idx, *], e_r[idx, *], ylog=e_log, yrange=ecc_range, ystyle=1, $
    xstyle=1, charsize=chrsize, symsize=0.4, xtitle='semimajor axis a- r!lREF!n', $
    ytitle='eccentricity e', noerase=1, position=pos11, xticks=2, psym=8, $
    thick=thck, xthick=thck, ythick=thck, charthick=thck, xrange=r_range
oplot, da[idx, *], e_r[idx, *], color=128, thick=thck
plots, [0, 0], ecc_range, linestyle=1, color=128, thick=thck

;Plot w versus a - r_ref
j = Nt/2
wpi_range = inputs.wpi_range
plot, da[j, *], w_r[j, *]/!dpi, yrange=wpi_range, ystyle=1, psym=8, $
    xstyle=1, charsize=chrsize, symsize=0.4, xtitle='radial distance a - r!lREF!n', $
    ytitle='long. of peri. !7!Sx!R!u!3__!n/!7p!3', noerase=1, position=pos10, $
    xticks=2, thick=thck, xthick=thck, ythick=thck, charthick=thck, xrange=r_range
oplot, da[j, *], w_r[j, *]/!dpi, color=128, thick=thck
oplot, r_range, t_s[sat_index]/!dpi + 0*r_range, color=128, thick=thck
plots, [0, 0], [-1, 1], color=128, thick=thck, linestyle=1

;indicate type of reference frame
ttl = 'nonrotating reference frame   '
if (corotate_sat eq 'y') then ttl = 'frame corotates w/satellite   '
if (corotate_ring eq 'y') then ttl = 'frame corotates w/ring !7!Sx!R!u!3__!n   '

;plot r - r_ref vs theta for all streamlines
longitude_wrap, t_r, r_r, t_w, r_w
longitude_wrap, t_r, a_r, t_w, a_w
dr = r_w - r_ref
da = a_w - r_ref
plot, t_w[*, 0]/!dpi, dr[*, 0], xstyle=1, charsize=chrsize, xrange=[-1, 1], $
    xtitle='!7h!3/!7p!3    (rad)', nodata=1, yrange=r_range, $
    ytitle='r(!7h!3) - r!lREF!n', noerase=1, ystyle=1, position=pos20, title=ttl, $
    thick=thck, xthick=thck, ythick=thck, charthick=thck
;oplot, t_s[sat_index]/!dpi + [0, 0], r_range, color=128, linestyle=1, thick=thck
for j = 0, Nr - 1 do begin $
    oplot, t_w[*, j]/!dpi, da[*, j], color=128, thick=thck & $
    oplot, t_w[*, j]/!dpi, dr[*, j], thick=thck & $
endfor
oplot, [-1, 1], [0, 0], color=128, linestyle=1, thick=thck
for j = 0, n_elements(sigma_radii) - 1 do $
    oplot, [-1, 1], sigma_radii[j] + [0, 0], color=128, linestyle=2, thick=thck
for j = 0, n_elements(sigma_lngtude) - 1 do $
    oplot, sigma_lngtude[j]*[1,1]/!dpi, r_range, color=128, linestyle=2, thick=thck
if (Nt le 24) then oplot, t_w/!dpi, dr, psym=8, symsize=0.5
oplot, t_s[sat_index]/!dpi + [0,0], r_range, color=128 & $

;Plot fractional change in system's total angular momentum
Lz = quantities.ang_mom
if ((tm gt 0) and (Lz[0] gt 0)) then begin $
    dLz = abs(Lz/Lz[0] - 1d)
    plot, time[0:tm], dLz[0:tm], charsize=chrsize, xtitle='time', $
        ytitle='!7D!3L!lz!n/L!lz!n', noerase=1, position=pos00, thick=thck, $
        xthick=thck, ythick=thck, charthick=thck
endif

;Or output plots to postscript file that is numbered sequentually and padded with zeros.
if (display eq 'postscript') then begin $
    device, /close & $
    frame_number = strtrim(string(tm), 2) & $
    N_zeros = fix(alog10(quantities.Noutput) + 1) - strlen(frame_number) & $
    for j = 1, N_zeros do frame_number = '0' + frame_number & $
    filename = 'postscript/' + frame_number + '.ps'
    spawn, 'mv idl.ps ' + filename & $
    set_plot, 'x' & $
endif

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


