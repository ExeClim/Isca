
       ! Socrates interface

       !Set tide-locked flux - should be set by namelist!
       soc_stellar_constant = 3500000.0
       fms_stellar_flux = soc_stellar_constant*COS(rlat)*COS(rlon)
       WHERE (fms_stellar_flux < 0.0) fms_stellar_flux = 0.0


       ! GCM mode
       hires_mode = .FALSE.

       ! LW calculation
       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
       soc_mode = .TRUE.
       CALL socrates_interface(Time, rlat, rlon, soc_mode, hires_mode,    &
            tg_tmp, t_surf, p_full, p_half, n_profile, n_layer,     &
            output_heating_rate, output_net_surf_sw_down, output_surf_lw_down, fms_stellar_flux )

       tg_tmp = tg_tmp + output_heating_rate*delta_t
       surf_lw_down(:,:) = output_surf_lw_down(:,:)



       ! SW calculation
       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
       soc_mode = .FALSE.
       CALL socrates_interface(Time, rlat, rlon, soc_mode, hires_mode,    &
            tg_tmp, t_surf, p_full, p_half, n_profile, n_layer,     &
            output_heating_rate, output_net_surf_sw_down, output_surf_lw_down, fms_stellar_flux )

       tg_tmp = tg_tmp + output_heating_rate*delta_t
       net_surf_sw_down(:,:) = output_net_surf_sw_down(:,:)
