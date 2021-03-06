;results.pro
;
; Read and inspect the output generated by epi_int.pro version 7.0
;

;Compile subroutines.pro, read output, and adjust Noutput by 1.
.run subroutines
@inputs
restore, filename=output_file
@inputs
Noutput = Noutput + 1

;disk surface density in gm/cm^2
surface_density_cgs = surface_density*M_planet_gm/(r1_cm^2)
print, 'ring surface density (gm/cm^2)  = ', surface_density_cgs
print, 'execution time (minutes)        = ', quantities.run_time_min

;Extract planetocentric coordinates and orbit elements of ring and satellites.
extract_results, quantities, coordinates, ID_r, mass_r, r_r, t_r, z_r, $
    vr_r, vt_r, vz_r, a_r, e_r, I_r, w_r, M_r, N_r, mass_s, r_s, t_s, z_s, $
    vr_s, vt_s, vz_s, a_s, e_s, I_s, w_s, M_s, N_s, sort_ID=0

;Show time evolution.
if (quantities.display eq 'monitor') then begin $
    window, xs=900, ys=900, retain=2 & $
    for tm = 0, Noutput - 1 do begin $
        display_results, inputs, quantities, tm, time, coordinates & $
        wait, 0.05 & $
    endfor & $
endif

