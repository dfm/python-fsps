MODULE fsps

    USE sps_vars
    ! USE sps_utils

    IMPLICIT NONE

    ! MAGIC: Make sure that these numbers stay the same as in `sps_vars.f90`.
    INTEGER, PARAMETER :: nspec2=1963,ntfull2=188,nbands2=83

    REAL, DIMENSION(ntfull2,nspec2) :: spec_ssp=0.
    REAL, DIMENSION(ntfull2) :: mass_ssp=0., lbol_ssp=0.
    REAL, DIMENSION(ntfull2) :: age=0.,mass_csp=0.,lbol_csp=0.,sfr=0.,mdust=0.
    REAL, DIMENSION(nbands2,ntfull2) :: mags=0.
    REAL, DIMENSION(nspec2,ntfull2) :: spec=0.

    CONTAINS

    SUBROUTINE get_dims(ntf, ns, nb)

        !f2py intent(out) ntf
        !f2py intent(out) ns
        !f2py intent(out) nb
        INTEGER :: ntf,ns,nb

        ntf = ntfull
        ns  = nspec
        nb  = nbands

    END SUBROUTINE get_dims

    SUBROUTINE compute(imf_type_in,zmet_in,sfh_in)

        INTEGER :: i
        TYPE(PARAMS) :: pset
        TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp

        !f2py intent(in) imf_type_in
        !f2py intent(in) zmet_in
        !f2py intent(in) sfh_in
        INTEGER :: imf_type_in,zmet_in,sfh_in

        imf_type   = imf_type_in
        pset%zmet  = zmet_in
        pset%sfh   = sfh_in
        pset%zred  = 0.0   !redshift
        pset%dust1 = 0.0   !dust parameter 1
        pset%dust2 = 0.0   !dust parameter 2
        pset%dell  = 0.0   !shift in log(L) for TP-AGB stars
        pset%delt  = 0.0   !shift in log(Teff) for TP-AGB stars
        pset%fbhb  = 0.0   !fraction of blue HB stars
        pset%sbss  = 0.0   !specific frequency of BS stars

        CALL SPS_SETUP(pset%zmet)
        CALL SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)
        CALL COMPSP(0,1,'',mass_ssp,lbol_ssp,spec_ssp,pset,ocompsp)

        DO i = 1, ntfull
            age(i)      = ocompsp(i)%age
            mass_csp(i) = ocompsp(i)%mass_csp
            lbol_csp(i) = ocompsp(i)%lbol_csp
            sfr(i)      = ocompsp(i)%sfr
            mdust(i)    = ocompsp(i)%mdust
            mags(:,i)   = ocompsp(i)%mags
            spec(:,i)   = ocompsp(i)%spec
        ENDDO

    END SUBROUTINE compute

END MODULE

