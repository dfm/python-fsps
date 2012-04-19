MODULE fsps

    USE sps_vars
    ! USE sps_utils

    IMPLICIT NONE

    ! MAGIC!
    INTEGER, PARAMETER :: nspec2=1963, ntfull2=188

    REAL, DIMENSION(ntfull2,nspec2) :: spec_ssp ! shape = ntfull,nspec
    REAL, DIMENSION(ntfull2) :: mass_ssp, lbol_ssp ! shape = ntfull

    CONTAINS

    SUBROUTINE pyfsps(mags)

        INTEGER, PARAMETER :: nspec2=1963, ntfull2=188, nbands2=83
        INTEGER :: i, j
        CHARACTER(100) :: file1=''
        TYPE(PARAMS) :: pset
        TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp

        !f2py intent(out) mags
        REAL, DIMENSION(nbands2,ntfull2) :: mags

        imf_type  = 1

        pset%zmet = 20


        CALL SPS_SETUP(pset%zmet) !read in the isochrones and spectral libraries

        pset%sfh   = 0     !set SFH to "SSP"
        pset%zred  = 0.0   !redshift
        pset%dust1 = 0.0   !dust parameter 1
        pset%dust2 = 0.0   !dust parameter 2
        pset%dell  = 0.0   !shift in log(L) for TP-AGB stars
        pset%delt  = 0.0   !shift in log(Teff) for TP-AGB stars
        pset%fbhb  = 0.0   !fraction of blue HB stars
        pset%sbss  = 0.0   !specific frequency of BS stars

        CALL SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)
        file1 = 'SSP.OUT'
        CALL COMPSP(3,1,file1,mass_ssp,lbol_ssp,spec_ssp,pset,ocompsp)

        DO i = 1, ntfull
            DO j = 1, nbands
                mags(j,i) = ocompsp(i)%mags(j)
            ENDDO
        ENDDO

    END SUBROUTINE pyfsps

END MODULE

