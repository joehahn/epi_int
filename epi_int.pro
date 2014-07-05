;epi_int.pro
;
;    version 7.0, by Joe Hahn, jhahn@spacescience.org, October 27, 2011.
;

;Compile routines in subroutines.pro and in the ILD library.
.run lap_coeff, subroutines, fx_root, spline, mean, moment

;Read input parameters and initialize system.
@inputs
initialize, inputs, quantities, time, ring, satellite, coordinates

;Execute main loop.
t_now = time[0]
for tm = 1l, Noutput do begin $

    ;Kick velocities due to pertubations and advance time a half step.
    t_step = dt/2 & $
    kick_velocities, quantities, t_now, t_step, ring, satellite & $

    for j = 1l, Nsteps_per_output do begin $

        ;Kick coordinates to account for central body's motion about center of mass,
        ;with a full timestep input here though the kick spans a half-timestep.
        kick_coordinates, dt, ring, satellite & $
 
        ;Ring and satellites execute unperturbed epicyclic orbital drift during timestep dt.
        drift, quantities, dt, ring & $
        drift, quantities, dt, satellite & $

        ;Kick coordinates to account for central body's motion about center of mass,
        ;with a full timestep input here though the kick spans a half-timestep.
        kick_coordinates, dt, ring, satellite & $

        ;Kick velocities due to perturbations and advance time a full step.
        t_step = dt & $
        kick_velocities, quantities, t_now, t_step, ring, satellite & $

        ;Break out of inner loop if Esc is pressed.
        if (get_kbrd(0) eq string(27B)) then break & $

    endfor & $

    ;Break out of main loop if Esc is pressed.
    if (get_kbrd(0) eq string(27B)) then break & $

    ;Kick velocities due to pertubations and advance time back a half step.
    t_step = -dt/2 & $
    kick_velocities, quantities, t_now, t_step, ring, satellite & $

    ;Convert mixed coordinates into planetocentric coordinates, sort disk paricles
    ;according to their longitudes, and store results in coordinates.
    store_results, tm, t_now, time, quantities, ring, satellite, coordinates & $

    ;Plot ring surface density, orbit elements, and streamlines.
    display_results, inputs, quantities, tm, time, coordinates & $

    ;If there are no errors, store system coordinates and time in file for possible restart.
    if (!error_state.name eq 'IDL_M_SUCCESS') then $
        save, filename='restart.dat', t_now, ring, satellite & $

endfor

;Display execution time
quantities.run_time_min = systime(/seconds)/60d - quantities.start_time_min
print, 'execution time (min) = ', quantities.run_time_min

;Save results in output_file.
save, filename=output_file, quantities, time, coordinates

