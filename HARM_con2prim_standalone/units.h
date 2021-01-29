/* == Definition of constants in CGS units == */
#define C_pi       (M_PI                        ) /* pi */
#define C_2pi      (2*M_PI                      ) /* 2*pi */
#define C_4pi      (4*M_PI                      ) /* 4*pi */
#define C_pi_2     (M_PI_2                      ) /* pi/2 */
#define C_euler    (0.5772156649015328606065    ) /* Euler's constant */
#define C_3_to_1_2 (1.73205080756887729352744634) /* 3^(1/2) */
#define C_sqrthalfpi (1.25331413731550025121    ) /* sqrt(pi/2) */
#define C_4pi__sqrt3_gamma1over3_2to1over3 (2.14952824153447863671) /* 4 pi / (sqrt(3) Gamma(1/3) 2^(1/3)) */
#define C_e        (4.80320680e-10              ) /* electron charge */
#define C_c        (2.99792458e10               ) /* speed of light */
#define C_me       (9.1093826e-28               ) /* electron mass */
#define C_mp       (1.67262171e-24              ) /* proton mass */
#define C_mn       (1.67492728e-24              ) /* neutron mass */
#define C_amu      (1.66053886e-24              ) /* atomic mass unit */
#define C_mAtom    (1.00794 * C_amu             ) /* atomic mass unit */
#define C_h        (6.6260693e-27               ) /* Planck constant */
#define C_k        (1.3806505e-16               ) /* Boltzmann constant */
#define C_G        (6.6742e-8                   ) /* Gravitational constant */
#define C_c2       (8.9875517873681764e20       ) /* c^2 */
#define C_c3       (2.6944002417373989539336e31 ) /* c^3 */
#define C_c5       (2.4216061708512206534320e52 ) /* c^5 */
#define C_e2       (2.3070795563566238161025e-19) /* e^2 */
#define C_e3       (1.1081380213233118296897e-28) /* e^3 */
#define C_me2      (8.2980851353182789400132e-55) /* me^2 */
#define C_me3      (7.5590432344986989028955e-82) /* me^3 */
#define C_me_c2    (8.1871047868450580395355e-7 ) /* me*c^2 */
#define C_e__2pi_me_c (2.799250542245159283e6   ) /* e / (2 Pi m_e c) : this factor times B is nu_c */
#define C_k__me_c2 (1.686372088724713698e-10    ) /* k / (me c^2) : this factor times temperature is Tdl */
#define C_4pi_m_c__3_e (2.3815898455876043578e-7) /* 4 pi m_e c / (3 e) */
#define C_1__me_plus_mp (5.9753838745824098154e23) /* 1 / (m_e + m_p) */
#define C_e3__8_sqrt3_pi2_me (8.895176161867948939e-4) /* e^3 / (8 Sqrt(3) pi^2 m_e) */
#define C_1AU      (1.49597870691e13            ) /* Astronomical Unit */
#define C_nA       (6.0221415e23                ) /* Avogadro number */ 
#define C_eV       (1.60217653e-12              ) /* electron volt in erg */
#define C_sigma    (5.670400e-5                 ) /* Stefan-Boltzmann constant */
#define C_alpha    (7.297352568e-3              ) /* fine structure constant */
#define C_rydberg  (2.17987209e-11              ) /* Rydberg constant times hc in erg */
#define C_sigmaT   (0.665245873e-24             ) /* Thomson cross section in cm^2 */ 
#define C_kappaT   (0.39750994621566976273      ) /* Thomson mass absorption coef. (sigmaT/m_H) in cm^2 g^(-1) */ 
#define C_Jy       (1.e-23                      ) /* Jansky (flux/freq. unit) in cgs */
#define C_bohr_mag (9.2740094980e-27            ) /* Bohr Magneton (electron) */ 
#define C_nucl_mag (5.0507834343e-30            ) /* Nuclear Magneton (proton) */ 
#define C_pc       (3.085678e18                 ) /* parsec */
#define C_yr       (31536000.                   ) /* No. of seconds in year */
#define C_day      (86400.                      ) /* No. of seconds in day  */
#define C_hr       (3600.                       ) /* No. of seconds in hour */
#define C_ly       (9.454254955488e17           ) /* light-year (distance in cm) */
#define C_mSUN     (1.99e33                     ) /* solar mass */
#define C_rSun     (6.96e10                     ) /* Radius of Sun */
#define C_LSun     (3.9e33                      ) /* Luminousity of Sun */
#define C_TSun     (5.78e3                      ) /* Temperature of Sun's photosphere */
#define C_MEarth   (5.976e27                    ) /* Earth's mass */
#define C_rEarth   (6.378e8                     ) /* Earth's radius */
#define C_dSgrA    (8.3e3 * C_pc                ) /* Distance from Earth to Sgr A*  */


#define atomicWeight (1.) /* all Helium & metal are broken down to H^+ due to high temp.'s */
