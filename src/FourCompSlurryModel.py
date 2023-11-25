import numpy as np


# Specific gravity terms
def Sf(Sl, Xf, Cv, Ss):
    """Carrier fluid specific gravity

    Args:
        Sl (-): SG of liquid (without solids)
        Xf (-): Carrier volume fraction
        Cv (-): total delivered volumetric concentration of solids (all fractions)
        Ss (-): SG of solids

    Returns:
        -: SG of liquid + Xf solids
    """
    return Sl + ((Xf * Cv * (Ss - Sl)) / (1 - Cv * (1 - Xf)))


def Sfp(Sl, Xf, Xp, Cv, Ss):
    """Carrier fluid + pseudo-homogeneous specific gravity

    Args:
        Sl (-): SG of liquid (without solids)
        Xf (-): Carrier volume fraction
        Xp (-): Pseudo-homogeneous volume fraction
        Cv (-): total delivered volumetric concentration of solids (all fractions)
        Ss (-): SG of solids

    Returns:
        -: SG of liquid + Xf and Xp solids
    """
    return Sl + (((Xf + Xp) * Cv * (Ss - Sl)) / (1 - Cv * (1 - Xf - Xp)))


def Sfph(Sl, Xf, Xp, Xh, Cv, Ss):
    """Carrier fluid + pseudo-homogeneous + heterogeneous specific gravity

    Args:
        Sl (-): SG of liquid (without solids)
        Xf (-): Carrier volume fraction
        Xp (-): Pseudo-homogeneous volume fraction
        Xh (-): Heterogeneous volume fraction
        Cv (-): total delivered volumetric concentration of solids (all fractions)
        Ss (-): SG of solids

    Returns:
        -: SG of liquid + Xf, Xp and Xh solids
    """
    return Sl + (((Xf + Xp + Xh) * Cv * (Ss - Sl)) / (1 - Cv * (1 - Xf - Xp - Xh)))


def Sm(Sl, Cv, Ss):
    """Mixture specific gravity

    Args:
        Sl (-): SG of liquid (without solids)
        Cv (-): total delivered volumetric concentration of solids (all fractions)
        Ss (-): SG of solids

    Returns:
        -: SG of total mixture
    """
    return Sl + Cv * (Ss - Sl)


# Concentration terms
def Cvf(Cv, Xf):
    """Carrier fluid volume concentration

    Args:
        Cv (-): total delivered volumetric concentration of solids (all fractions)
        Xf (-): Carrier volume fraction


    Returns:
        -: Volume concentration of Xf solids within the mixture of liquid + Xf particles.
    """
    return (Xf * Cv) / (1 - Cv * (1 - Xf))


def Cvp(Cv, Xf, Xp):
    """Carrier fluid volume concentration

    Args:
        Cv (-): total delivered volumetric concentration of solids (all fractions)
        Xf (-): Carrier volume fraction
        Xp (-): Pseudo-homogeneous volume fraction


    Returns:
        -: Volume concentration of Xp solids within the mixture of liquid + Xp and Xf particles.
    """
    return (Xp * Cv) / (1 - Cv * (1 - Xf - Xp))


def Cvh(Cv, Xf, Xp, Xh):
    """Pseudo-homogeneous volume concentration

    Args:
        Cv (-): total delivered volumetric concentration of solids (all fractions)
        Xf (-): Carrier volume fraction
        Xp (-): Pseudo-homogeneous volume fraction
        Xh (-): Heterogeneous volume fraction


    Returns:
        -: Volume concentration of Xp solids within the mixture of liquid + Xh, Xp and Xf particles.
    """
    return (Xh * Cv) / (1 - Cv * (1 - Xf - Xp - Xh))


def Cvs(Cv, Xs):
    """Fully stratified volume concentration

    Args:
        Cv (-): total delivered volumetric concentration of solids (all fractions)
        Xs (-): Fully stratified volume fraction

    Returns:
        -: Volume concentration of Xs solids within the total mixture.
    """
    return Xs * Cv


# Pressure gradient - Carrier fluid fraction
def i_f(Sf, f_darcy, Vm, D, g=9.81):
    """Pressure gradient for the Carrier Fluid - particles < 40 μm

    Args:
        Sf (-): Carrier fluid specific gravity
        ff (-): Darcy-Weisbach friction factor.
        Vm (m/s): average slurry mixture velocity in the pipe
        D (m): pipe inner diameter
        g (m/s/s): acceleration of gravity. Defaults to 9.81.

    Returns:
        m/m: Carrier Fluid pressure gradient
    """
    return Sf * f_darcy * ((Vm**2) / (2 * g * D))


def mu_f(mu_l, Cvf):
    """The “Carrier Fluid” fraction - viscosity correction

    Args:
        mu_l (Pa.s): dynamic viscosity of the liquid (without solids)
        Cvf (-): Carrier fluid volume concentration

    Returns:
        Pa.s: dynamic viscosity of the “Carrier Fluid” (liquid + fines)
    """
    return mu_l * (1 + 2.5 * Cvf + 10 * Cvf**2 + 0.0019 * np.exp(20 * Cvf))


# Misc
def re_number(Sf, mu_f, D, Vm):
    """Reynolds number

    Args:
        Sf (-): Carrier fluid specific gravity
        mu_f (Pa.s): Dynamic viscosity of the “Carrier Fluid” (liquid + fines)
        D (m): Inner diameter of pipe
        Vm (m/s): Velocity of fluid / mixture

    Returns:
        -: Reynolds number
    """
    return (Sf * 1000 * Vm * D) / mu_f


def swamee_jain(e, D, Sf, mu_f, Vm):
    """Swamee-Jain equation. Applicable for 10E5 < Re < 10E7

    Args:
        e (m): Pipe roughness
        D (m): Inner pipe diameter
        Sf (-): Carrier fluid specific gravity
        mu_f (Pa.s): Dynamic viscosity of the “Carrier Fluid” (liquid + fines)
        D (m): Inner diameter of pipe
        Vm (m/s): Velocity of fluid / mixture

    Returns:
        -: Darcy-Weisbach friction factor
    """

    re = re_number(Sf, mu_f, D, Vm)

    return 0.25 / (np.log10((e / (3.7 * D)) + (5.74 / (re**0.9)))) ** 2


# Pressure gradient - Pseudo-homogeneous fraction
def A_prime(Xf, Xp, Af=1.0, Ap=0.5):
    """Empirical parameter, accounts for near wall, hydrodynamic lift effects which may be experienced by the “pseudo-homogeneous” particles.

    Args:
        Xf (-): Carrier volume fraction
        Xp (-): Pseudo-homogeneous volume fraction
        Af (float, optional): Based on laboratory testing. Defaults to 1.0.
        Ap (float, optional): Based on laboratory testing. Defaults to 0.5.
    Returns:
        A_prime (-): Empirical parameter.
    """
    return 1 - (Af * Xf + Ap * Xp)


def i_p(Sfp, Sf, i_f, A_prime):
    """Pressure gradient for the Pseudo-homogeneous fraction - 40 μm < particles < 200 μm

    Args:
        Sfp (-): Carrier fluid + pseudo-homogeneous specific gravity
        Sf (-): Carrier fluid specific gravity
        A_prime (-): Empirical parameter, accounts for near wall, hydrodynamic lift effects which may be experienced by the “pseudo-homogeneous” particles.
        i_f (m/m): Pressure gradient for the Carrier Fluid - particles < 40 μm

    Returns:
        m/m: Pseudo-homogeneous pressure gradient
    """
    return A_prime * (Sfp - Sf) * i_f / Sf


# Pressure gradient - Heterogeneous fraction
# why is hetrogeneous C" instead of B"? :(
def C_prime(Xf, Xp, V100_s, Vm, Vsm_h, Cf=1.0, Cp=0.5, nC=0.5):
    """The empirical parameter for particle-particle interactions between components is based on the volumes of the finer Xf and Xp fractions. Applicable to the Heterogeneous fraction.

    Args:
        Xf (-): Carrier volume fraction
        Xp (-): Pseudo-homogeneous volume fraction
        V100_s (m/s): Pseudo-homogeneous velocity for 0.015D sized particles
        Vm (m/s): Velocity of fluid / mixture
        Vsm_h (m/s): Maximum velocity at limit of stationary deposition for d50h sized particles
        Cf (float, optional): Based on laboratory testing. Defaults to 1.0.
        Cp (float, optional): Based on laboratory testing. Defaults to 0.5.
        nC (float, optional): Based on laboratory testing. Defaults to 0.5.

    Returns:
        C_prime (-): Empirical parameter.
    """

    return 1 - (Cf * Xf + Cp * Xp) * ((V100_s - Vm) / (V100_s - Vsm_h)) ** nC


def i_h(C_prime, Sfph, Sfp, Sf, Ss, Vm, d50_h, mu_f, mu_w=0.001, M=1.0, mu_s=0.5):
    """The excess pressure gradient contribution for the Heterogeneous component is based on the Wilson V50 model. 200 μm < particles < 0.015D

    Args:
        C_prime (-): Empirical parameter.
        Sfph (-): Carrier fluid + pseudo-homogeneous + heterogeneous specific gravity
        Sfp (-): Carrier fluid + pseudo-homogeneous specific gravity
        Sf (-): Carrier fluid specific gravity
        Ss (-): SG of solids
        Vm (m/s): Velocity of fluid / mixture
        d50_h (m): d50 of the Xh fraction solids
        mu_f (Pa.s): Dynamic viscosity of the “Carrier Fluid” (liquid + fines)
        mu_w (float, optional): viscosity of water @ 20 degC. Defaults to 0.001.
        M (float, optional): 1.0 when using the four component model. Defaults to 1.0.
        mu_s (float, optional): sliding friction coef. between particles and pipe wall. Defaults to 0.5 for harder, angular particles. 0.4 for softer, rounded particles.

    Returns:
        m/m: Heterogeneous Pressure gradient
    """
    vr = mu_f / (mu_w * Sf)

    V50_h = 44.1 * (d50_h**0.35 / vr**0.25) * ((Ss - Sfp) / 1.65)

    return C_prime * (mu_s / 2) * (Sfph - Sfp) * (V50_h / Vm) ** M


# Pressure gradient - Fully Stratified fraction
def B_prime(Xf, Xp, Xh, V100_s, Vm, Vsm_s, Bf=1.0, Bp=1.0, Bh=0.5, nB=0.5):
    """The empirical parameter for particle-particle interactions between components is based on the volumes of the finer Xf, Xp and Xh fractions.

    Args:
        Xf (-): Carrier volume fraction
        Xp (-): Pseudo-homogeneous volume fraction
        Xh (_type_): _description_
        V100_s (m/s): Pseudo-homogeneous velocity for 0.015D sized particles
        Vm (m/s): Velocity of fluid / mixture
        Vsm_s (m/s): Maximum velocity at limit of stationary deposition for 0.015D sized particles

        Bf (float, optional): Based on laboratory testing. Defaults to 1.0.
        Bp (float, optional): Based on laboratory testing. Defaults to 1.0.
        Bh (float, optional): Based on laboratory testing. Defaults to 0.5.
        nB (float, optional): Based on laboratory testing. Defaults to 0.5.

    Returns:
        B_prime (-): Empirical parameter.

    """
    return 1 - (Bf * Xf + Bp * Xp + Bh * Xh) * ((V100_s - Vm) / (V100_s - Vsm_s)) ** nB


def i_s(B_prime, Sfph, Ss, Cvs, Vsm_s, Vm, mu_s=0.5):
    """The excess pressure gradient contribution for the “Fully Stratified” component is based on the model proposed by Wilson and Addie, (1995) - 0.015D < d

    Args:
        B_prime (-): Empirical parameter.
        Sfph (-): Carrier fluid + pseudo-homogeneous + heterogeneous specific gravity
        Ss (-): SG of solids
        Cvs (-): Volume concentration of Xs solids within the total mixture.
        Vsm_s (m/s): Maximum velocity at limit of stationary deposition for 0.015D sized particles
        Vm (m/s): Velocity of fluid / mixture
        mu_s (float, optional): sliding friction coef. between particles and pipe wall. Defaults to 0.5 for harder, angular particles. 0.4 for softer, rounded particles.

    Returns:
        m/m: Fully Stratified Pressure gradient
    """
    return B_prime * 2 * mu_s * Cvs * (Ss - Sfph) * (Vsm_s / Vm) ** 0.25


def Vsm(d, D, Ss, Sf, mu_f, e, mu_s=0.5, g=9.81):
    """Determination of Vsm (the deposition velocity)

    Args:
        d (m): the particle diameter. d50_h for Vsm_h , or ds for Vsm_s
        D (m): Inner pipe diameter
        Ss (-): SG of solids
        Sf (-): Carrier fluid specific gravity
        mu_f (Pa.s): Dynamic viscosity of the “Carrier Fluid” (liquid + fines)
        mu_s (float, optional): sliding friction coef. between particles and pipe wall. Defaults to 0.5 for harder, angular particles. 0.4 for softer, rounded particles.
        g (m/s/s, optional): acceleration of gravity. Defaults to 9.81.
        e (m): Pipe roughness

    Returns:
        m/s: Deposition velocity
    """

    d_mm = d * 1000

    Vsm_norm = (
        8.8
        * ((mu_s * ((Ss - Sf) / (0.66 * Sf))) ** 0.55)
        * ((D**0.7 * d_mm**1.75) / (d_mm**2 + 0.11 * D**0.7))
    )

    V_initalise = 0.01  # m/s

    for i in range(10):
        # Iterates through 10 times, probably a better way to do this but it works.

        ff = swamee_jain(e, D, Sf, mu_f, V_initalise)
        V_initalise = (0.018 / ff) ** 0.13 * (2 * g * D * ((Ss / Sf) - 1)) ** 0.5

    Vsm_max = V_initalise

    return min(Vsm_norm, Vsm_max)


# Determination of V100_s (the pseudo-homogeneous velocity)
def V100_s(D, Ss, Sf, g=9.81):
    """Pseudo-homogeneous velocity of the 0.15D sized particle.

    Args:
        D (m): Inner pipe diameter
        Ss (-): SG of solids
        Sf (-): Carrier fluid specific gravity
        g (m/s/s, optional): acceleration of gravity. Defaults to 9.81.

    Returns:
        m/s: Pseudo-homogeneous velocity
    """
    ds = 0.015 * D

    xi = 0.4 * ds**-0.04

    vt_s = 1.73 * xi * (g * ds * (Ss - Sf)) ** 0.5

    return (1800 * g * D * vt_s) ** (1.0 / 3)
