module driver

  use sps_vars; use sps_utils
  implicit none
  save

  !f2py intent(hide) pset
  type(PARAMS) :: pset

  !f2py intent(hide) ocompsp
  type(COMPSPOUT), dimension(ntfull) :: ocompsp

  integer :: is_setup=0

  !f2py intent(hide) has_ssp
  integer, dimension(nz) :: has_ssp=0

  !f2py intent(hide) has_ssp_age
  integer, dimension(nz,nt) :: has_ssp_age=0

  
contains

  subroutine setup(compute_vega_mags0, vactoair_flag0)

    ! Load all the data files/templates into memory.

    implicit none

    integer, intent(in) :: compute_vega_mags0, vactoair_flag0
         

    compute_vega_mags = compute_vega_mags0
    vactoair_flag = vactoair_flag0
    call sps_setup(-1)
    is_setup = 1

    ! We will only compute mags when asked for through get_mags.
    pset%mag_compute=0

  end subroutine

  subroutine set_ssp_params(imf_type0,imf_upper_limit0,imf_lower_limit0,&
                            imf1,imf2,imf3,vdmc,mdave,dell,&
                            delt,sbss,fbhb,pagb,add_stellar_remnants0,&
                            tpagb_norm_type0,add_agb_dust_model0,agb_dust,&
                            redgb,agb,masscut,fcstar,evtype,smooth_lsf0)
 
    ! Set the parameters that affect the SSP computation.

    implicit none

    integer, intent(in) :: imf_type0,add_stellar_remnants0,tpagb_norm_type0,&
                           add_agb_dust_model0,smooth_lsf0
    double precision, intent(in) :: imf_upper_limit0, imf_lower_limit0,&
                                    imf1,imf2,imf3,vdmc,mdave,dell,&
                                    delt,sbss,fbhb,pagb,agb_dust,&
                                    redgb,agb,masscut,fcstar,evtype

    imf_type=imf_type0
    imf_upper_limit=imf_upper_limit0
    imf_lower_limit=imf_lower_limit0
    add_stellar_remnants=add_stellar_remnants0
    tpagb_norm_type=tpagb_norm_type0
    add_agb_dust_model=add_agb_dust_model0
    smooth_lsf=smooth_lsf0
    pset%imf1=imf1
    pset%imf2=imf2
    pset%imf3=imf3
    pset%vdmc=vdmc
    pset%mdave=mdave
    pset%dell=dell
    pset%delt=delt
    pset%sbss=sbss
    pset%fbhb=fbhb
    pset%pagb=pagb
    pset%agb_dust=agb_dust 
    pset%redgb=redgb
    pset%agb=agb
    pset%masscut=masscut
    pset%fcstar=fcstar
    pset%evtype=evtype

    has_ssp(:) = 0
    has_ssp_age(:,:) = 0
    
  end subroutine

  subroutine set_csp_params(smooth_velocity0,redshift_colors0,&
                            compute_light_ages0,nebemlineinspec0,&
                            dust_type0,add_dust_emission0,add_neb_emission0,&
                            add_neb_continuum0,cloudy_dust0,add_igm_absorption0,&
                            zmet,sfh,wgp1,wgp2,wgp3,tau,&
                            const,tage,fburst,tburst,dust1,dust2,&
                            logzsol,zred,pmetals,dust_clumps,frac_nodust,&
                            dust_index,dust_tesc,frac_obrun,uvb,mwr,&
                            dust1_index,sf_start,sf_trunc,sf_slope,&
                            duste_gamma,duste_umin,duste_qpah,&
                            sigma_smooth,min_wave_smooth,max_wave_smooth,&
                            gas_logu,gas_logz,igm_factor,fagn,agn_tau)

    ! Set all the parameters that don't affect the SSP computation.

    implicit none
    
    integer, intent(in) :: smooth_velocity0,redshift_colors0,&
                           compute_light_ages0,nebemlineinspec0,&
                           dust_type0,add_dust_emission0,add_neb_emission0,&
                           add_neb_continuum0,cloudy_dust0,add_igm_absorption0,&
                           zmet,sfh,wgp1,wgp2,wgp3
    double precision, intent(in) :: tau,&
                            const,tage,fburst,tburst,dust1,dust2,&
                            logzsol,zred,pmetals,dust_clumps,frac_nodust,&
                            dust_index,dust_tesc,frac_obrun,uvb,mwr,&
                            dust1_index,sf_start,sf_trunc,sf_slope,&
                            duste_gamma,duste_umin,duste_qpah,&
                            sigma_smooth,min_wave_smooth,max_wave_smooth,&
                            gas_logu,gas_logz,igm_factor,fagn,agn_tau

    smooth_velocity=smooth_velocity0
    redshift_colors=redshift_colors0
    compute_light_ages=compute_light_ages0
    nebemlineinspec=nebemlineinspec0
    dust_type=dust_type0
    add_dust_emission=add_dust_emission0
    add_neb_emission=add_neb_emission0
    add_neb_continuum=add_neb_continuum0
    cloudy_dust=cloudy_dust0
    add_igm_absorption=add_igm_absorption0

    pset%zmet=zmet
    pset%sfh=sfh
    pset%wgp1=wgp1
    pset%wgp2=wgp2
    pset%wgp3=wgp3

    pset%tau=tau
    pset%const=const
    pset%tage=tage
    pset%fburst=fburst
    pset%tburst=tburst
    pset%dust1=dust1
    pset%dust2=dust2
    pset%logzsol=logzsol
    pset%zred=zred
    pset%pmetals=pmetals
    pset%dust_clumps=dust_clumps
    pset%frac_nodust=frac_nodust
    pset%dust_index=dust_index
    pset%dust_tesc=dust_tesc
    pset%frac_obrun=frac_obrun
    pset%uvb=uvb
    pset%mwr=mwr
    pset%dust1_index=dust1_index
    pset%sf_start=sf_start
    pset%sf_trunc=sf_trunc
    pset%sf_slope=sf_slope
    pset%duste_gamma=duste_gamma
    pset%duste_umin=duste_umin
    pset%duste_qpah=duste_qpah
    pset%sigma_smooth=sigma_smooth
    pset%min_wave_smooth=min_wave_smooth
    pset%max_wave_smooth=max_wave_smooth
    pset%gas_logu=gas_logu
    pset%gas_logz=gas_logz
    pset%igm_factor=igm_factor
    pset%fagn=fagn
    pset%agn_tau=agn_tau
    
  end subroutine

  subroutine ssps

    ! Loop over the metallicity grid and compute all the SSPs.

    implicit none
    integer :: zi
    do zi=1,nz
      call ssp(zi)
    enddo

  end subroutine

  subroutine ssp(zi)

    ! Compute a SSP at a single metallicity.

    implicit none
    integer, intent(in) :: zi
    pset%zmet = zi
    call ssp_gen(pset, mass_ssp_zz(:,zi),lbol_ssp_zz(:,zi),&
         spec_ssp_zz(:,:,zi))
    if (minval(pset%ssp_gen_age) .eq. 1) then
       has_ssp(zi) = 1
    endif
    has_ssp_age(zi,:) = pset%ssp_gen_age
    
  end subroutine

  subroutine compute_zdep(ns,n_age,ztype)

    ! Compute the full CSP (and the SSPs if they aren't already cached).
    ! After interpolation in metallicity

    implicit none
    integer, intent(in) :: ns,n_age,ztype
    double precision, dimension(ns,n_age) :: spec
    double precision, dimension(n_age) :: mass,lbol
    integer :: zlo,zmet
    double precision :: zpos
    character(100) :: outfile

    if (ztype .eq. 0) then
       ! Build the SSP for one metallicity, then feed to compsp
       zmet = pset%zmet
       if (has_ssp(zmet) .eq. 0) then
          call ssp(zmet)
       endif
       !mass = mass_ssp_zz(:,zmet)
       !lbol = lbol_ssp_zz(:,zmet)
       !spec = spec_ssp_zz(:,:,zmet)
       call compsp(0,1,outfile,mass_ssp_zz(:,zmet),lbol_ssp_zz(:,zmet),&
            spec_ssp_zz(:,:,zmet),pset,ocompsp)
    endif

    if (ztype .eq. 1) then
       zpos = pset%logzsol
       ! Find the bracketing metallicity indices and generate ssps if
       ! necessary, then interpolate, and feed the result to compsp
       zlo = max(min(locate(log10(zlegend/0.0190),zpos),nz-1),1)
       do zmet=zlo,zlo+1
          if (has_ssp(zmet) .eq. 0) then
             call ssp(zmet)
          endif
       enddo
       call ztinterp(zpos,spec,lbol,mass)
       call compsp(0,1,outfile,mass,lbol,spec,pset,ocompsp)
    endif
    
    if (ztype .eq. 2) then
       zpos = pset%logzsol
       ! Build the SSPs for *every* metallicity if necessary, then
       ! comvolve with the MDF, and then feed to compsp
       do zmet=1,nz
          if (has_ssp(zmet) .eq. 0) then
             call ssp(zmet)
          endif
       enddo
       call ztinterp(zpos,spec,lbol,mass,zpow=pset%pmetals)
       call compsp(0,1,outfile,mass,lbol,spec,pset,ocompsp)
    endif

    if (ztype .eq. 3) then
       ! Build the SSPs for *every* metallicity and feed all of them to compsp
       ! for z-dependent tabular
       do zmet=1,nz
          if (has_ssp(zmet) .eq. 0) then
             call ssp(zmet)
          endif
       enddo
       call compsp(0,nz,outfile,mass_ssp_zz,lbol_ssp_zz,&
            spec_ssp_zz,pset,ocompsp)
    endif

    
  end subroutine

  subroutine get_spec(ns,n_age,spec_out)

    ! Get the grid of spectra for the computed CSP at all ages.

    implicit none
    integer :: i
    integer, intent(in) :: ns,n_age
    double precision, dimension(n_age,ns), intent(out) :: spec_out
    do i=1,n_age
      spec_out(i,:) = ocompsp(i)%spec
    enddo

  end subroutine

  subroutine get_mags(ns,n_age,n_bands,z_red,mc,mags)

    ! Get the photometric magnitudes in all the recognized bands.
    
    implicit none
    integer :: i
    integer, intent(in) :: ns, n_age, n_bands
    double precision, intent(in) :: z_red
    integer, dimension(n_bands), intent(in) :: mc
    double precision, dimension(n_age,n_bands), intent(out) :: mags
    double precision, dimension(ns) :: tspec
    do i=1,n_age
      tspec = ocompsp(i)%spec
      call getmags(z_red,tspec,mags(i,:),mc)
    enddo

  end subroutine

  subroutine interp_ssp(ns,zpos,tpos,spec,mass,lbol)

    ! Return the SSPs interpolated to the target metallicity 
    !(zpos) and target age (tpos)
    
    implicit none

    integer, intent(in) :: ns
    double precision, intent(in) :: zpos
    double precision, intent(in) :: tpos

    double precision, dimension(ns,1), intent(inout) :: spec
    double precision, dimension(1), intent(inout) :: mass,lbol

    double precision, dimension(nt) :: time
    
    integer :: zlo,zmet,tlo

    zlo = max(min(locate(log10(zlegend/0.0190),zpos),nz-1),1)
    time = timestep_isoc(zlo,:)
    tlo = max(min(locate(time,tpos),nt-1),1)

    do zmet=zlo,zlo+1
       if ((has_ssp_age(zmet,tlo) .eq. 0) .or. (has_ssp_age(zmet,tlo+1) .eq. 0)) then
          pset%ssp_gen_age = 0
          pset%ssp_gen_age(tlo:tlo+1) = 1
          call ssp(zmet)
          pset%ssp_gen_age = 1
       endif
    enddo

    call ztinterp(zpos,spec,lbol,mass,tpos=tpos)

    end subroutine

  subroutine smooth_spectrum(ns,wave,spec,sigma_broad,minw,maxw)
    
    ! Smooth the spectrum by a gaussian of width sigma_broad
    
    implicit none
    integer, intent(in) :: ns
    double precision, intent(in) :: sigma_broad,minw,maxw
    double precision, dimension(ns), intent(in) :: wave
    double precision, dimension(ns), intent(inout) :: spec
    
    call smoothspec(wave,spec,sigma_broad,minw,maxw)

  end subroutine

  subroutine stellar_spectrum(ns,mact,logt,lbol,logg,phase,ffco,lmdot,wght,spec_out)
    
    ! Get a stellar spectrum for a given set of parameters
    
    implicit none
    integer :: i
    integer, intent(in) :: ns
    double precision, intent(in) :: mact, logt, lbol, logg, phase, ffco, lmdot, wght
    double precision, dimension(ns), intent(inout) :: spec_out

    call getspec(pset,mact,logt,lbol,logg,phase,ffco,lmdot,wght,spec_out)
    
  end subroutine 
  
  subroutine get_ssp_spec(ns,n_age,n_z,ssp_spec_out,ssp_mass_out,ssp_lbol_out)

    ! Return the contents of the ssp spectral array,
    ! regenerating the ssps if necessary

    implicit none
    integer, intent(in) :: ns,n_age,n_z
    integer :: zi
    double precision, dimension(ns,n_age,n_z), intent(inout) :: ssp_spec_out
    double precision, dimension(n_age,n_z), intent(inout) :: ssp_mass_out, ssp_lbol_out
    do zi=1,nz
       if (has_ssp(zi) .eq. 0) then
          call ssp(zi)
       endif
    enddo

    ssp_spec_out = spec_ssp_zz
    ssp_mass_out = mass_ssp_zz
    ssp_lbol_out = lbol_ssp_zz

  end subroutine

  subroutine set_sfh_tab(ntab, age, sfr, met)

    ! Fill the sfh_tab array

    implicit none
    integer, intent(in) :: ntab
    double precision, dimension(ntab), intent(in) :: age, sfr, met
    ntabsfh = ntab
    sfh_tab(1,1:ntabsfh) = age
    sfh_tab(2, 1:ntabsfh) = sfr
    sfh_tab(3, 1:ntabsfh) = met

  end subroutine set_sfh_tab

  subroutine set_ssp_lsf(nsv, sigma, wlo, whi)

    ! Fill the lsfinfo structure

    implicit none
    integer, intent(in) :: nsv
    double precision, dimension(nsv), intent(in) :: sigma
    double precision, intent(in) :: wlo, whi
    lsfinfo%minlam = wlo
    lsfinfo%maxlam = whi
    lsfinfo%lsf = sigma

  end subroutine set_ssp_lsf
  
  subroutine get_setup_vars(cvms, vta_flag)

    implicit none
    integer, intent(out) :: cvms, vta_flag
    cvms = compute_vega_mags
    vta_flag = vactoair_flag

  end subroutine

  subroutine get_nz(n_z)

    ! Get the number of metallicity bins (hard coded in sps_vars).
    implicit none
    integer, intent(out) :: n_z
    n_z = nz

  end subroutine

  subroutine get_zlegend(n_z,z_legend)

    ! Get the available metallicity values.
    implicit none
    integer, intent(in) :: n_z
    double precision, dimension(n_z), intent(out) :: z_legend
    z_legend = zlegend

  end subroutine

  subroutine get_timefull(n_age,timefull)

    ! Get the actual time steps of the SSPs.
    implicit none
    integer, intent(in) :: n_age
    double precision, dimension(n_age), intent(out) :: timefull

    timefull = time_full

  end subroutine

  subroutine get_ntfull(n_age)

    ! Get the total number of time steps (hard coded in sps_vars).
    implicit none
    integer, intent(out) :: n_age
    n_age = ntfull

  end subroutine

  subroutine get_nspec(ns)

    ! Get the number of wavelength bins in the spectra (hard coded in
    ! sps_vars).
    implicit none
    integer, intent(out) :: ns
    ns = nspec

  end subroutine

  subroutine get_nbands(nb)

    ! Get the number of known filters (hard coded in sps_vars).
    implicit none
    integer, intent(out) :: nb
    nb = nbands

  end subroutine

  subroutine get_nemline(nline)

    ! Get the number of emission lines (hard coded in sps_vars).
    implicit none
    integer, intent(out) :: nline
    nline = nemline 

  end subroutine

  subroutine get_emlambda(nline,em_lambda)

    ! Get the emission line wavelengths
    implicit none
    integer, intent(in) :: nline
    double precision, dimension(nline), intent(out) :: em_lambda
    if (vactoair_flag .eq. 1) then
       em_lambda = vactoair(nebem_line_pos)
    else
       em_lambda = nebem_line_pos
    endif

  end subroutine

  subroutine get_lambda(ns,lambda)

    ! Get the grid of wavelength bins.
    implicit none
    integer, intent(in) :: ns
    double precision, dimension(ns), intent(out) :: lambda
    if (vactoair_flag .eq. 1) then
       lambda = vactoair(spec_lambda)
    else
       lambda = spec_lambda
    endif
    
  end subroutine

  subroutine get_libraries(isocname,specname)

    implicit none

    character(4), intent(out) :: isocname
    character(5), intent(out) :: specname
    isocname = isoc_type
    specname = spec_type

  end subroutine  

  subroutine get_isochrone_dimensions(n_age,n_mass)

    implicit none

    ! Get the dimensions of the produced isochrones.
    integer, intent(out) :: n_age,n_mass
    n_age = nt
    n_mass = n_mass

  end subroutine

  subroutine get_nmass_isochrone(zz, tt, nmass)

    implicit none

    ! Get the number of masses included in a specific isochrone.
    integer, intent(in) :: zz,tt
    integer, intent(out) :: nmass
    nmass = nmass_isoc(zz,tt)

  end subroutine

  subroutine get_stats(n_age,nline,age,mass_csp,lbol_csp,sfr,mdust,mformed,emlines)

    implicit none

    ! Get some stats about the computed SP.
    integer :: i
    integer, intent(in) :: n_age,nline
    double precision, dimension(n_age), intent(out) :: age,mass_csp,&
                                                       lbol_csp,sfr,mdust,&
                                                       mformed
    double precision, dimension(n_age,nline), intent(out) :: emlines

    do i=1,n_age
      age(i)      = ocompsp(i)%age
      mass_csp(i) = ocompsp(i)%mass_csp
      lbol_csp(i) = ocompsp(i)%lbol_csp
      sfr(i)      = ocompsp(i)%sfr
      mdust(i)    = ocompsp(i)%mdust
      mformed(i)  = ocompsp(i)%mformed
      emlines(i,:)  = ocompsp(i)%emlines
    enddo

  end subroutine

  subroutine get_filter_data(nb, wave_eff, mag_vega, mag_sun)

    !get info about the filters
    implicit none
    integer, intent(in) :: nb
    double precision, dimension(nb), intent(out) :: wave_eff,mag_vega,mag_sun
    wave_eff = filter_leff
    mag_vega = magvega - magvega(1)
    mag_sun = magsun

  end subroutine

  subroutine write_isoc(outfile)

    implicit none

    character(100), intent(in)  :: outfile
    
    call write_isochrone(outfile, pset)

  end subroutine

end module
