#include <cctk.h>
#include <cctk_Arguments.h>

// check sanity of EOS parameters when EOS is used
void CCTK_FCALL
CCTK_FNAME(EOS_Omni_Check_EOS_Params)(const CCTK_INT *eos_handle);

CCTK_INT EOS_Omni_GetHandle_(CCTK_STRING name)
{
    // invalid name, better would have been a negative number
    CCTK_INT eos_handle = 0;

    if (CCTK_EQUALS(name, "2D_Polytrope"))
        eos_handle = 1;
    else if (CCTK_EQUALS(name, "Ideal_Fluid"))
        eos_handle = 2;
    else if (CCTK_EQUALS(name, "Hybrid"))
        eos_handle = 3;
    else if (CCTK_EQUALS(name, "nuc_eos"))
        eos_handle = 4;
    else if (CCTK_EQUALS(name, "cold_tabulated"))
        eos_handle = 5;
    else if (CCTK_EQUALS(name, "barotropic_tabulated"))
        eos_handle = 6;

    if(eos_handle)
        CCTK_FNAME(EOS_Omni_Check_EOS_Params)(&eos_handle);

    return eos_handle;
}
