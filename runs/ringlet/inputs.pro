;inputs.pro  
;
;    called by v7.0 of epi_int.pro.
;

;Set restart='y' to resume intgration from last saved output.
restart = 'n'

;Set the number of ring particles, where Nr=number of streamlines along radial
;direction and Nt=number of particles in tangential direction along one streamline.
Nr = 16l
Nt = 32l

;Set time parameters, where dt=timestep, Nsteps_per_output=number of timesteps per
;output, and Noutput=number of outputs +1 to store initial state.
dt = 0.1d
Nsteps_per_output = 2000l
Noutput = 100l

;Set difference_method='approximate' (faster but less accurate) or 'polynomial' (slower
;but more accurate). Set precision=1 for drift step to use equations accurate to e^1,
;or precision=2 (recommended) to use e^2 equations.
difference_method = 'polynomial'
precision = 2

;Set the physical parameters that are used to convert model units to into physical
;units of time, length, and mass, with M_planet_gm=the central planet's mass in grams
;and r1_cm=distance in centimeters that corresponds to where the simulated distance r=1.
M_planet_gm = 5.6846d29
r1_cm = 1.178145d10

;Set the ring's initial conditions, with a_r_min and a_r_max = semimajor axes of the
;ring's inner and outer edges, e_r_init=initial eccentricity,
;I_r_init=initial inclintation, w_r_init=initial longitude of periapse,
;and N_r_init=initial longitude of ascending node. All angles are in radians
;and length are in units of r1. Set corotate_ring='y' to make all plots
;co-precess with the ring's mean longitude of periapse.
da = 1.0d-4
a_r_max = 1 + da
a_r_min = 1 - 10*da
da1 = (a_r_max - a_r_min)/200	;smallest spacing used when 'expo-outwards' or 'expo-inwards'
radial_spacing = 'uniform'	;'uniform' or 'expo-outwards' or 'expo-inwards'
ringlet_disk = 'ringlet'
e_r_init = 0d
de_r_init = 0d
I_r_init = 0d
w_r_init = 0d
N_r_init = 0d
corotate_ring = 'n'

;Radiative boundary conditions
outer_rbc = 'n'
inner_rbc = 'n'
N_sl_rbc = 5

;Set ring_gravity='y' to turn on ring gravity with surface_density=ring surface density
;if it were initially circular.
ring_gravity = 'y'
surface_density = 5.0d-8

;Set arrays of satellites' initial planetocentric orbit elements a_s,e_s,I_s,w_s,N_s,M_s. 
;If timescale t_grow_sat>0 then the satellites' masses grows exponentially over that
;timescale to final mass mass_s_final. Set satellite_gravity='y' to turn on satellite
;gravity, corotate_sat='y' for all plotted longitudes to corotate with satellite[sat_idx]
;longitude. Set epsilon_LR=1 when interested in the inner Lindblad resonance of satellite
;sat_idx, and epsilon_LR=-1 when interested in its outer Lindblad resonance. If
;fourier_satellite='y' then satellite sat_index has zero mass (and thus doesn't alter
;the system's barycenter), yet its m_LR^th fourier component of gravity perturbs the ring.
a_s = [ 1.587401052d ]/1.00574985d
e_s = [ 0d ]
I_s = [ 0d ]
w_s = [ 0d ]
M_s = [ 1d ]
N_s = [ 0d ]
mass_s_final = [ 0d ]
t_grow_sat = -5d3
sat_index = 0
corotate_sat = 'n'
m_LR = 2
epsilon_LR = 1
satellite_gravity = 'y'
fourier_satellite = 'n'

;Parameters for the central planet's oblateness, where Rp_cm=central planet's radius
;in centimeters and J2.
Rp_cm = 60330d5
J2 = 0d

;Set pressure='p' to turn on, and set ring particles' initial dispersion speed c.
pressure = 'n'
c0 = 2*!dpi*surface_density

;Set viscosity='y' to turn on, and set values for shear and bulk viscosity.
;Also set confine_inner(outer)_edge='y' to confine ring's inner(outer) edge by zeroing
;the viscous tangential acceleration there.
viscosity = 'n'
viscosity_law = 'constant'		; 'constant' or 'collisions' or 'wakes' or 'fake c'
shear_viscosity = 1.0d-13
bulk_viscosity = shear_viscosity
optical_depth = 1d
confine_outer_edge = 'y'
confine_inner_edge = 'y'

;Drag force is turned on by setting drag_force='y' and setting frag coefficient drag_coeff.
;The drage force decays away over exponential timescale drag_timescale when positive, and
;is constant when given a negative value.
drag_force = 'n'
drag_coeff_inner = 4d-4
drag_coeff_outer = 4d-4
drag_transition_radius = -0.00035
drag_transition_width = 1d-4
drag_timescale = -4.0d4

;Parameters for the display. Set display='monitor' to display plots on monitor,
;or 'postscript' to send plots to numbered postscript files.
;Set r_ref=reference radial distance such that all radial distances are
;measured from r_ref, units of vertical axes measured in plot_units='r1' or 'km',
;eccentricity axis ecc_axis='log' or 'linear', ecc_range=vertical range in
;eccentricity plot, wpi_range=range of ring longitude-of-periapse/PI plot,
;sd_range=range of surface density plot, sd_lngtude=longitudes where
;radial surface density cuts are plotted, sd_radii=
;r_range=radial range of epicyclic amplitude versus longitude plot, 
;z_range=vertical range of height versus longitude plot.
display = 'monitor'
r_ref = 1d
plot_units = 'r1'
ecc_axis = 'log'
ecc_range = [1.0d-10, 1.0d-4]
wpi_range = [-1, 1]
sd_range = [0.0, 1.5]
sd_radii = [0]
sd_longitude = [-0.25, 0.25]*!dpi
r_range = [-0.0012, 0.0005]
z_range = [-0.00002, 0.00002]

;Name of file where simulation output data and restart data is stored
output_file  = 'results.dat'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;store all inputs in a structure
inputs = { restart:restart, Nr:Nr, Nt:Nt, dt:dt, Nsteps_per_output:Nsteps_per_output, $
           Noutput:Noutput, difference_method:difference_method, $
           M_planet_gm:M_planet_gm, r1_cm:r1_cm, precision:precision, $
           a_s:a_s, e_s:e_s, I_s:I_s, w_s:w_s, M_s:M_s, N_s:N_s, $
           mass_s_final:mass_s_final, t_grow_sat:t_grow_sat, sat_index:sat_index, $
           corotate_sat:corotate_sat, m_LR:m_LR, epsilon_LR:epsilon_LR, $
           satellite_gravity:satellite_gravity, fourier_satellite:fourier_satellite, $
           radial_spacing:radial_spacing, da1:da1, ringlet_disk:ringlet_disk, $
           a_r_min:a_r_min, a_r_max:a_r_max, e_r_init:e_r_init, de_r_init:de_r_init, $
           I_r_init:I_r_init, w_r_init:w_r_init, N_r_init:w_r_init, $
           corotate_ring:corotate_ring, $
           inner_rbc:inner_rbc, outer_rbc:outer_rbc, N_sl_rbc:N_sl_rbc, $
           ring_gravity:ring_gravity, surface_density: surface_density, $
           pressure:pressure, c0:c0, viscosity:viscosity, viscosity_law:viscosity_law, $
           shear_viscosity:shear_viscosity, bulk_viscosity:bulk_viscosity, $
           optical_depth:optical_depth, confine_inner_edge:confine_inner_edge, $
           confine_outer_edge:confine_outer_edge, Rp_cm:Rp_cm, J2:J2, $ 
           drag_force:drag_force, drag_coeff_inner:drag_coeff_inner, $
           drag_coeff_outer:drag_coeff_outer, drag_transition_radius:drag_transition_radius, $
           drag_transition_width:drag_transition_width, drag_timescale:drag_timescale, $
           display:display, r_ref:r_ref, plot_units:plot_units, ecc_axis:ecc_axis, $
           ecc_range:ecc_range, wpi_range:wpi_range, sd_range:sd_range, $
           sd_radii:sd_radii, sd_longitude:sd_longitude, $
           r_range:r_range, z_range:z_range, output_file:output_file }

