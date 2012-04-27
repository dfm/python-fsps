! This is a module that will be used in Python as an iterface to FSPS.
! Rock on, eh?

module driver

    use sps_vars; use nrtype; use sps_utils
    implicit none

    !f2py intent(hide) pset
    type(PARAMS) :: pset

    !f2py intent(hide) ocompsp
    type(COMPSPOUT), dimension(ntfull) :: ocompsp

    contains

    subroutine setup
        ! Load all the data files/templates into memory.
        call sps_setup(-1)
    end subroutine

    subroutine ssps(imf,imf1,imf2,imf3,vdmc,mdave,dell,delt,sbss,fbhb,pagb)
        ! Calculate all of the SSPs in one go.
        integer :: zi

        integer, intent(in) :: imf
        real, intent(in) :: imf1, imf2, imf3, vdmc, mdave
        real, intent(in) :: dell, delt, sbss, fbhb, pagb

        imf_type   = imf
        pset%imf1  = imf1
        pset%imf2  = imf2
        pset%imf3  = imf3
        pset%vdmc  = vdmc
        pset%mdave = mdave
        pset%dell  = dell
        pset%delt  = delt
        pset%sbss  = sbss
        pset%fbhb  = fbhb
        pset%pagb  = pagb

        ! Loop over the metallicities and generate the SSPs.
        do zi=1,nz
            pset%zmet = zi
            call ssp_gen(pset, mass_ssp_zz(zi,:),lbol_ssp_zz(zi,:),&
                     spec_ssp_zz(zi,:,:))
        enddo
    end subroutine

    subroutine compute(dust,zmet,sfh,tau,const,fburst,tburst,dust_tesc,dust1,&
            dust2,dust_clumps,frac_no_dust,dust_index,mwr,wgp1,wgp2,wgp3,tage)
        ! Compute the stellar population given a set of physical parameters.
        integer, intent(in) :: dust, zmet, sfh
        real, intent(in) :: tau,const,fburst,tburst,dust_tesc,dust1,dust2,&
            dust_clumps,frac_no_dust,dust_index,mwr,tage
        integer, intent(in) :: wgp1,wgp2,wgp3

        dust_type = dust
        pset%zmet = zmet
        pset%sfh = sfh
        pset%tau = tau
        pset%const = const
        pset%tage = tage
        pset%fburst = fburst
        pset%tburst = tburst
        pset%dust_tesc = dust_tesc
        pset%dust1 = dust1
        pset%dust2 = dust2
        pset%dust_clumps = dust_clumps
        pset%frac_nodust = frac_no_dust
        pset%dust_index = dust_index
        pset%mwr = mwr
        pset%wgp1 = wgp1
        pset%wgp2 = wgp2
        pset%wgp3 = wgp3

        call compsp(0,1,'',mass_ssp_zz(zmet,:),lbol_ssp_zz(zmet,:),&
            spec_ssp_zz(zmet,:,:),pset,ocompsp)
    end subroutine

    subroutine get_ntfull(nt)
        ! Get the total number of time steps (hard coded in sps_vars).
        integer, intent(out) :: nt
        nt = ntfull
    end subroutine

    subroutine get_nspec(ns)
        ! Get the number of wavelength bins in the spectra.
        integer, intent(out) :: ns
        ns = nspec
    end subroutine

    subroutine get_lambda(ns,lambda)
        ! Get the grid of wavelength bins.
        integer, intent(in) :: ns
        real, dimension(ns), intent(out) :: lambda
        lambda = spec_lambda
    end subroutine

    subroutine get_stats(nt,age,mass_csp,lbol_csp,sfr,mdust)
        ! Get some stats about the computed SP.
        integer :: i
        integer, intent(in) :: nt
        real, dimension(nt), intent(out) :: age,mass_csp,lbol_csp,sfr,mdust

        do i=1,nt
            age(i)      = ocompsp(i)%age
            mass_csp(i) = ocompsp(i)%mass_csp
            lbol_csp(i) = ocompsp(i)%lbol_csp
            sfr(i)      = ocompsp(i)%sfr
            mdust(i)    = ocompsp(i)%mdust
        enddo
    end subroutine

    subroutine get_spec(ns,nt,spec_out)
        ! Get the set of spectra as a function of time.
        integer :: i
        integer, intent(in) :: ns,nt
        real, dimension(nt,ns), intent(out) :: spec_out
        do i=1,nt
            spec_out(i,:) = ocompsp(i)%spec
        enddo
    end subroutine

end module

