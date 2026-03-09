# importing packages
import pyomo.environ as pyo
import numpy as np
from matplotlib import pyplot as plt
import math
import scipy.stats as stats
from scipy.optimize import fsolve, bisect
import pandas as pd


# Shedlovsky model for predicting the equivalent conductivity of salt (or ionic) solutions
def _shedlovsky(conc, temp, epsilon, eta, lambda_0, a, z_1, z_2, lambda_0_cation, lambda_0_anion):
    """
    Calculates equivalent conductivity of ionic solutions

    Parameters
    ----------
    conc: list or Pandas series
        A list or Pandas series containing the concentrations of the salt in M
    temp: int or float
        Temperature of the solution in K
    epsilon: int or float
        Dielectric constant of the solvent
    eta: int or float
        Viscosity of the solvent in poise
    lambda_0: int or float
        Limiting equivalent conductivity of the salt in cm^2.S/equiv
    a: int or float
        Distance of closest approach of ions in cm
    z_1: int
        Valency of the cation
    z_2: int
        Valency of the anion
    lambda_0_cation: int or float
        Limiting equivalent conductivity of the cation in cm^2.S/equiv
    lambda_0_anion: int or float
        Limiting equivalent conductivity of the anion in cm^2.S/equiv

    Returns
    -------
    equiv_cond: list
        A list containing the equivalent conductivity of the solution at various salt concentrations in cm^2.S/equiv
    """
    # calculate the B and q constants
    B = 50.29 * (10 ** 8) * (epsilon * temp) ** (-0.5)
    q = np.abs(z_1 * z_2) * (lambda_0_cation + lambda_0_anion) / (
                (np.abs(z_1) + np.abs(z_2)) * (np.abs(z_2) * lambda_0_cation + np.abs(z_1) * lambda_0_anion))

    # calculate the relaxation correction, B1, and the hydrodynamic correction, B2
    B_1 = 2.801 * 10 ** 6 * np.abs(z_1 * z_2) * q / ((epsilon * temp) ** (3 / 2) * (1 + math.sqrt(q)))
    B_2 = 41.25 * (np.abs(z_1) + np.abs(z_2)) / (eta * (epsilon * temp) ** (1 / 2))

    # ion concentration
    cation_conc = conc
    Cl_conc = [np.abs(z_1) * conc[i] for i in range(len(conc))]

    # calculate the ionic strength
    I = np.zeros(len(conc))
    for i in range(len(conc)):
        conc_all = [cation_conc[i], Cl_conc[i]]
        z = [z_1, z_2]
        I[i] = 1 / 2 * (sum(conc_all[j] * z[j] ** 2 for j in range(len(conc_all))))

    # calculate the salt equivalent conductivity
    equiv_cond = []
    for i in range(len(conc)):
        equiv_cond.append(lambda_0 - ((B_1 * lambda_0 + B_2) * math.sqrt(I[i]) / (1 + a * B * math.sqrt(I[i]))))

    return equiv_cond


# variant Shedlovsky model for predicting the specific conductivity of salt (or ionic) solutions
def variant_shedlovsky(conc, temp, epsilon, eta, lambda_0, a, z_1, z_2, lambda_0_cation, lambda_0_anion):
    """
    Calculates specific conductivity from equivalent conductivity

    Parameters
    ----------
    conc: list or Pandas series
        A list or Pandas series containing the concentrations of the salt in M
    temp: int or float
        Temperature of the solution in K
    epsilon: int or float
        Dielectric constant of the solvent
    eta: int or float
        Viscosity of the solvent in poise
    lambda_0: int or float
        Limiting equivalent conductivity of the salt in cm^2.S/equiv
    a: int or float
        Distance of closest approach of ions in cm
    z_1: int
        Valency of the cation
    z_2: int
        Valency of the anion
    lambda_0_cation: int or float
        Limiting equivalent conductivity of the cation in cm^2.S/equiv
    lambda_0_anion: int or float
        Limiting equivalent conductivity of the anion in cm^2.S/equiv

    Returns
    -------
    specific_cond_calc: list
        A list containing the specific conductivity of the solution at various salt concentrations in milli.S/cm
    """

    # calculate the salt equivalent conductivity
    equiv_cond = _shedlovsky(conc, temp, epsilon, eta, lambda_0, a, z_1, z_2, lambda_0_cation, lambda_0_anion)

    # cation concentration
    cation_conc = conc

    # convert the cation concentrations from M to equivalent per litre
    equiv_conc = [cation_conc[i] * np.abs(z_1) for i in range(len(conc))]

    # calculate the salt specific conductivity in Siemen per cm
    specific_cond = []
    for i in range(len(conc)):
        specific_cond.append((equiv_conc[i] / 1000) * equiv_cond[i])

    # specific conductivity in milli.S/cm
    specific_cond_calc = [specific_cond[i] * 10 ** 3 for i in range(len(conc))]

    return specific_cond_calc


# Mean spherical approximation (MSA) model for multi-salt specific conductivity predictions
def msa_transport(valency, diameters, diff_coeff, temp, eta, epsilon, lambda_0,
                  salt_1_conc, salt_2_conc=None, salt_3_conc=None):
    """Calculates the specific conductivity of single, binary, and ternary salt solutions

    Parameters
    ----------
    valency: list
        A list of the valency of ions in decreasing order
    diameters: list
        A list of the corresponding hard sphere diameter of the ions in m
    diff_coeff: list
        A list of the corresponding diffusion coefficient of the ions at infinite dilution in m^2/s
    temp: int or float
        Temperature of the solutions in K
    eta: int or float
        Viscosity of the solvent in Pa.s
    epsilon: int or float
        Relative permittivity of the solvent
    lambda_0: list
        A list of the limiting molar conductivity of the ions in S.m^2/mol
    salt_1_conc: Pandas series or list
        A list or Pandas series containing the concentrations of salt 1 in mM
    salt_2_conc: Pandas series or list
        A list or Pandas series containing the concentrations of salt 2 in mM
    salt_3_conc: Pandas series or list
        A list or Pandas series containing the concentrations of salt 3 in mM

    Returns
    -------
    cond_calc_con: list
        A list containing the specific conductivity of the solution at various salt concentrations in milli.S/cm"""

    # constants
    Avogadros_num = 6.022 * 10 ** 23  # Avogadro's constant
    k_B = 1.381 * 10 ** (-23)  # Boltzmann constant in J/K
    charge = 1.602 * 10 ** (-19)  # elementary charge in C
    epsilon_0 = 8.854 * 10 ** (-12)  # permittivity of free space in F/m
    Faraday = 96500  # Faraday's constant in C/mol

    # number of data points
    ndata = len(salt_1_conc)

    # number of species
    n_species = len(valency)

    # species charge
    ion_charge = [valency[i] * charge for i in range(len(valency))]

    # converting molar concentrations to number density
    if salt_2_conc is None and salt_3_conc is None:
        number_density_salt_1 = [salt_1_conc[i] * Avogadros_num for i in range(ndata)]

        # valency of cation
        valency_cation = valency[0]

        # evaluating the number density of individual species
        number_density_cat_1 = number_density_salt_1
        number_density_an = [valency_cation * number_density_salt_1[i] for i in range(ndata)]
    elif salt_2_conc is not None and salt_3_conc is None:
        number_density_salt_1 = [salt_1_conc[i] * Avogadros_num for i in range(ndata)]
        number_density_salt_2 = [salt_2_conc[i] * Avogadros_num for i in range(ndata)]

        # valency of cations
        valency_cation_1 = valency[0]
        valency_cation_2 = valency[1]

        # evaluating the number density of individual species
        number_density_cat_1 = number_density_salt_1
        number_density_cat_2 = number_density_salt_2
        number_density_an = [valency_cation_1 * number_density_salt_1[i] + valency_cation_2 * number_density_salt_2[i]
                             for i in range(ndata)]
    else:
        number_density_salt_1 = [salt_1_conc[i] * Avogadros_num for i in range(ndata)]
        number_density_salt_2 = [salt_2_conc[i] * Avogadros_num for i in range(ndata)]
        number_density_salt_3 = [salt_3_conc[i] * Avogadros_num for i in range(ndata)]

        # valency of cations
        valency_cation_1 = valency[0]
        valency_cation_2 = valency[1]
        valency_cation_3 = valency[2]

        # evaluating the number density of individual species
        number_density_cat_1 = number_density_salt_1
        number_density_cat_2 = number_density_salt_2
        number_density_cat_3 = number_density_salt_3
        number_density_an = [valency_cation_1 * number_density_salt_1[i] + valency_cation_2 * number_density_salt_2[i]
                             + valency_cation_3 * number_density_salt_3[i] for i in range(ndata)]

    # calculate the bulk conductivity of the salt solution
    all_cond_calc = np.zeros(ndata)
    for n in range(ndata):  # loop through the salts concentration
        if salt_2_conc is None and salt_3_conc is None:
            # number density of species
            n_cat = number_density_cat_1[n]
            n_an = number_density_an[n]
            number_density = [n_cat, n_an]
        elif salt_2_conc is not None and salt_3_conc is None:
            # number density of species
            n_cat_1 = number_density_cat_1[n]
            n_cat_2 = number_density_cat_2[n]
            n_an = number_density_an[n]
            number_density = [n_cat_1, n_cat_2, n_an]
        else:
            # number density of species
            n_cat_1 = number_density_cat_1[n]
            n_cat_2 = number_density_cat_2[n]
            n_cat_3 = number_density_cat_3[n]
            n_an = number_density_an[n]
            number_density = [n_cat_1, n_cat_2, n_cat_3, n_an]

        # calculate the species ionic mobility (omega) and the relative ionic strength (mew)
        omega = np.zeros(n_species)
        mew = np.zeros(n_species)
        for a in range(n_species):
            omega[a] = diff_coeff[a] / (k_B * temp)
            mew[a] = number_density[a] * ion_charge[a] ** 2 / sum(
                number_density[j] * ion_charge[j] ** 2 for j in range(n_species))

        # calculate the Debye length from equation (5)
        kappa = math.sqrt(
            sum(number_density[l] * ion_charge[l] ** 2 / (epsilon * epsilon_0 * k_B * temp) for l in range(n_species)))

        # calculate the mean mobility from equation (7)
        omega_bar = sum(mew[j] * omega[j] for j in range(n_species))

        # calculate the transport number of all species from equation (8)
        transport_num = np.zeros(n_species)
        for j in range(n_species):
            transport_num[j] = mew[j] * omega[j] / omega_bar

        # calculate the value of alpha from the function containing alpha
        alpha_scaled = np.zeros(n_species)
        def func_alpha(alpha_scaled):
            """
            Evaluates residuals of equation 12

            Arguments:
                alpha_scaled = alpha/omega_bar

            Returns:
                residual of equation 12
            """
            return ((alpha_scaled) * sum(
                transport_num[k] / ((omega[k] / omega_bar) ** 2 - (alpha_scaled) ** 2) for k in range(n_species)))

        # since alpha is scaled by omega_bar, omega need to be scaled by same factor
        omega_scaled = np.zeros(n_species)
        for k in range(n_species):
            omega_scaled[k] = omega[k] / omega_bar

        # updating the values of alpha_scaled
        for p in range(n_species):
            if p == 0:
                alpha_scaled[p] = bisect(func_alpha, 0, 0.991 * omega_scaled[p])
            else:
                alpha_scaled[p] = bisect(func_alpha, 1.001 * omega_scaled[p - 1], 0.991 * omega_scaled[p])

        # re-scaling alpha_scaled back to alpha
        alphas = [alpha_scaled[p] * omega_bar for p in range(n_species)]

        # calculate the q and N values of all species
        q = np.zeros(n_species)
        N = np.zeros(n_species)
        for p in range(n_species):
            q[p] = sum(omega_bar * transport_num[r] / (omega[r] + alphas[p]) for r in range(n_species))
            N[p] = math.sqrt(1 / sum(
                transport_num[r] * omega[r] ** 2 / (omega[r] ** 2 - alphas[p] ** 2) ** 2 for r in range(n_species)))

        # ksi matrix
        ksi = np.zeros((n_species, n_species))

        # calculate the hydrodynamic correction
        delta_vel_divide_vel = np.zeros(n_species)
        for i in range(n_species):
            # update the hydrodynamic correction; equation (16) divide equation (2)
            prefactor = -(Faraday ** 2) * np.abs(valency[i]) / (
                        12 * math.pi * epsilon_0 * epsilon * eta * k_B * temp * Avogadros_num * lambda_0[i])
            delta_vel_divide_vel[i] = prefactor * sum(number_density[j] * ion_charge[j] ** 2 *
                                                      ((np.exp(kappa * (diameters[i] - 0.5 * (
                                                                  diameters[i] + diameters[j]))) / (
                                                                    kappa * (1 + kappa * diameters[i]))) +
                                                       (np.exp(kappa * (diameters[j] - 0.5 * (
                                                                   diameters[i] + diameters[j]))) / (
                                                                    kappa * (1 + kappa * diameters[j]))))
                                                      for j in range(n_species))

            # update the ksi matrix
            for p in range(n_species):
                ksi[i][p] = N[p] * omega[i] / (omega[i] ** 2 - alphas[p] ** 2)

        # calculate the conductivity of the species
        cond_species = np.zeros(n_species)

        # calculate the relaxation correction
        delta_k_divide_k = np.zeros(n_species)
        for i in range(n_species):
            delta_k_divide_k_sum = 0.0
            for p in range(n_species):
                inner_delta_k_divide_k_sum = 0.0
                for j in range(n_species):
                    for a in range(n_species):
                        # calculate the average diameter
                        sigma_aj = 0.5 * (diameters[a] + diameters[j])

                        # first fraction of summation term
                        term1 = transport_num[j] * ksi[j][p] * mew[a] * (
                                ion_charge[a] * omega[a] - ion_charge[j] * omega[j]) / (
                                        ion_charge[a] * ion_charge[j] * (omega[a] + omega[j]))

                        # second fraction of summation term
                        term2 = math.sinh(kappa * math.sqrt(q[p]) * sigma_aj) / (kappa * math.sqrt(q[p]) * sigma_aj)

                        # evaluate the integral term
                        integral_prefactor = - (charge ** 2 * valency[a] * valency[j]) / (
                                    8 * math.pi * epsilon_0 * epsilon * k_B * temp)

                        integral_term1 = np.exp(kappa * diameters[j]) / (1 + kappa * diameters[j])
                        integral_term2 = np.exp(kappa * diameters[a]) / (1 + kappa * diameters[a])

                        integral_exp_term = np.exp(-kappa * (1 + np.sqrt(q[p])) * sigma_aj) / (
                                    kappa * (1 + np.sqrt(q[p])))

                        integral = integral_prefactor * (integral_term1 + integral_term2) * integral_exp_term

                        inner_delta_k_divide_k_sum += term1 * term2 * integral

                delta_k_divide_k_sum += ksi[i][p] * inner_delta_k_divide_k_sum

            # update the relaxation correction
            delta_k_divide_k[i] = -(kappa ** 2) * ion_charge[i] * delta_k_divide_k_sum / 3

            # update the conductivity of the species
            cond_species[i] = (charge ** 2) * number_density[i] * diff_coeff[i] * (valency[i] ** 2) * (
                    1 + delta_vel_divide_vel[i]) * (1 + delta_k_divide_k[i]) / (k_B * temp)

        # calculate the bulk conductivity of the solution in S/m
        all_cond_calc[n] = sum(cond_species)

    # convert the calculated conductivity from S/m to milli.S/cm
    cond_calc_con = [all_cond_calc[i] * 10 ** 3 / 100 for i in range(ndata)]

    print(f"Relaxation correction: {delta_k_divide_k}")
    print(f"Hydrodynamic correction: {delta_vel_divide_vel}")
    print(f"omega_bar: {omega_bar}")
    print(f"Alpha: {alphas}")
    print(f"Mew: {mew}")
    print(f"Transport number: {transport_num}")
    print(f"Omega: {omega}")
    print(f"ksi: {ksi}")

    return cond_calc_con
