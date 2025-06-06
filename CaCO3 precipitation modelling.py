#%%
"""
Calcium Carbonate Precipitation Model
-------------------------------------
This model simulates the precipitation of calcium carbonate polymorphs (calcite, aragonite, vaterite)
in CO2-H2O solutions at various temperatures, following CO2 degassing.

The model includes thermodynamic equilibria and kinetic processes for the CaCO3-CO2-H2O system.

Features:
- Multiple CaCO3 polymorphs (calcite, aragonite, vaterite)
- Temperature-dependent thermodynamic equilibria
- CO2 degassing kinetics with pH-dependent pathways
- Precipitation kinetics with nucleation and growth
- Experimental data comparison and visualization

Applications:
- CO2 geological sequestration
- Ocean acidification studies
- Material synthesis optimization
- Water treatment processes
"""

"""
QUICK START GUIDE
================

Basic Usage:
>>> model = CaCO3PrecipitationModel()
>>> model.set_polymorph(calcite=True)  # Enable calcite precipitation
>>> results = model.simulate(initial_Ca=0.01, P_CO2=1.0)  # Run simulation
>>> model.plot_results(results)  # Visualize results

Parameter Modification:
>>> model.set_parameters(**{
...     'degassing.k_degas': 1e-3,      # Adjust degassing rate
...     'precipitation.SI_critical': 1.2,  # Change nucleation threshold
...     'system_vars.Temperature': 35    # Set different temperature
... })

Load Experimental Data:
>>> model.load_experimental_data('experiment.csv', 'Exp_25C', temperature=25)

Common Applications:
- CO2 sequestration: High P_CO2, variable temperature
- Ocean chemistry: Low P_CO2, 0-30°C, aragonite/calcite
- Industrial crystallization: Controlled conditions, rapid kinetics
- Natural systems: Atmospheric P_CO2, seasonal temperature variation
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, curve_fit
import pandas as pd
import os

class CaCO3PrecipitationModel:
    """
    Model of calcium carbonate precipitation in aqueous systems
    including CO2 degassing and CaCO3 precipitation kinetics
 
    This model simulates the coupled processes of:
    1. CO2 degassing from supersaturated solutions
    2. Carbonate speciation equilibria (CO2 ⇌ HCO3- ⇌ CO3²-)
    3. CaCO3 precipitation kinetics and nucleation
    
    The model solves a system of ODEs representing mass balance and
    kinetic rate equations, with thermodynamic constraints.
    
    Quick Start:
    -----------
    >>> model = CaCO3PrecipitationModel()
    >>> model.set_polymorph(calcite=True)
    >>> results = model.simulate(initial_Ca=0.01, P_CO2=1.0)
    >>> model.plot_results(results)
    """
    
    def __init__(self):
        # Initialize the model with default settings
        self.setup_variables()
        self.setup_equilibria()
        self.exp_data = {}
    
    def setup_variables(self):
        """
        1. Define solution variables, their meanings, and reasonable ranges
        """
        # System Variables - These are modified during simulation
        self.system_vars = {
            'pH': 6.0,            # Solution pH (range: 0-14)
            'Ca': 0.01,           # Calcium ion concentration (mol/L) 
            'DIC': 0.01,          # Dissolved inorganic carbon (mol/L) 
            'CO2_aq': 0.005,      # Dissolved CO2 (mol/L) 
            'HCO3': 0.005,        # Bicarbonate ion (mol/L) 
            'CO3': 0.0001,        # Carbonate ion (mol/L) 
            'Temperature': 25,    # Solution temperature (°C) 
            'P_CO2': 0.0004       # CO2 partial pressure (atm) 
        }
        
        # Thermodynamic Parameters - Constants for each polymorph
        """
        Thermodynamic Parameters for CaCO3 Polymorphs
        Based on Plummer & Busenberg (1982) solubility data
        Each polymorph has different stability and solubility:
        - Calcite: Most stable, lowest solubility
        - Aragonite: Metastable, intermediate solubility  
        - Vaterite: Least stable, highest solubility
        """
        self.polymorph_params = {
            'calcite': {
                'enabled': True,                # Enable/disable this polymorph
                'name': 'Calcite',              # Display name
                'color': 'blue',                # Color for plotting
                # Parameters for log K calculation (log K = A + B*T + C/T + D*log(T))
                'eq_const_params': {
                    'A': -171.9065,
                    'B': -0.077993,
                    'C': 2839.319,
                    'D': 71.595
                },
                'K_sp_25C': 10**(-8.48)        # Solubility product at 25°C
            },
            'aragonite': {
                'enabled': False,
                'name': 'Aragonite',
                'color': 'green',
                'eq_const_params': {
                    'A': -171.9773,
                    'B': -0.077993,
                    'C': 2903.293,
                    'D': 71.595
                },
                'K_sp_25C': 10**(-8.336)
            },
            'vaterite': {
                'enabled': False,
                'name': 'Vaterite',
                'color': 'red',
                'eq_const_params': {
                    'A': -172.1295,
                    'B': -0.077993,
                    'C': 3074.688,
                    'D': 71.595
                },
                'K_sp_25C': 10**(-7.913)
            }
        }
        
        # Kinetic Parameters - These control the rates of processes
        """
        Kinetic Parameters - Control rates of chemical and physical processes
        These parameters determine how fast reactions occur (not just equilibrium)
        """
        self.kinetic_params = {
            # CO2 degassing parameters
            'degassing': {
                'enabled': True,               # Enable/disable CO2 degassing
                'k_degas': 2e-3,               # Base degassing rate constant (s⁻¹) (range: 1e-5 - 1e-2)
                
                # pH-dependent degassing pathways
                    # CO2 degassing occurs via multiple pH-dependent pathways:
                    # 1. Acid catalysis: H⁺ + HCO3⁻ → CO2 + H2O (fast at low pH)
                    # 2. Neutral pathway: HCO3⁻ → CO2 + OH⁻ (moderate at neutral pH)  
                    # 3. Base catalysis: OH⁻ + CO2 → HCO3⁻ (fast at high pH)
                'acid_catalysis': {
                    'enabled': True,           # Enable/disable acid catalysis
                    'k_acid': 5.0e-4           # Acid-catalyzed rate constant (range: 1e-5 - 1e-2)
                },
                'neutral_pathway': {
                    'enabled': True,           # Enable/disable neutral pathway
                    'k_neutral': 0.05          # Neutral pH rate constant (range: 0.001 - 1.0)
                },
                'base_catalysis': {
                    'enabled': True,           # Enable/disable base catalysis
                    'k_base': 5.0e-2           # Base-catalyzed rate constant (range: 1e-4 - 1e-1)
                },
                
                'non_equilibrium': {
                    'enabled': True,           # Enable/disable non-equilibrium effects
                    'k_factor': 0.005          # Non-equilibrium factor (range: 0.001 - 0.1)
                }
            },
            
            # Precipitation parameters
            'precipitation': {
                'enabled': True,               # Enable/disable precipitation
                'k_precip': 7e-9,              # Precipitation rate constant (L mol⁻¹ s⁻¹) (range: 1e-10 - 1e-6)
                'n_reaction': 1.0,             # Reaction order (range: 0.5 - 3.0)
                'SI_critical': 1.1,            # Critical saturation index for nucleation (range: 1.0 - 2.0)
                'induction_time': 45000,       # Induction period before precipitation (s) (range: 0 - 100000)
                
                # pH response
                'pH_lag': {
                    'enabled': True,           # Enable/disable pH lag effect
                    'tau': 1                   # pH lag time constant (s) (range: 0.1 - 10000)
                }
            }
        }

    def setup_equilibria(self):
        """
        2. Define equilibrium equations and temperature-dependent constants
        """
        T = self.system_vars['Temperature']
        T_K = T + 273.15  # Convert to Kelvin
        """
        Calculate temperature-dependent equilibrium constants
        All equilibria are fundamental to carbonate chemistry:

        Chemical Reactions:
        1. CO2(g) ⇌ CO2(aq)           [Henry's Law]
        2. CO2(aq) + H2O ⇌ H⁺ + HCO3⁻ [First dissociation]  
        3. HCO3⁻ ⇌ H⁺ + CO3²⁻         [Second dissociation]
        4. Ca²⁺ + CO3²⁻ ⇌ CaCO3(s)    [Precipitation]
        5. H2O ⇌ H⁺ + OH⁻             [Water autoionization]
        """
        # CO2-H2O Equilibria
        # ------------------
        # CO2(g) <--> CO2(aq)
        # First dissociation: K_a1 = [H⁺][HCO3⁻]/[CO2(aq)] [mol/L]
        self.constants = {}
        self.constants['K_H'] = 0.034 * np.exp(-0.0235 * (T - 25))
        
        # CO2(aq) + H2O <--> H+ + HCO3-
        # First dissociation constant
        pK1 = 6.35 - 0.0054 * (T - 25)
        self.constants['K_a1'] = 10**(-pK1)
        
        # HCO3- <--> H+ + CO3(2-)
        # Second dissociation constant
        pK2 = 10.33 - 0.009 * (T - 25)
        self.constants['K_a2'] = 10**(-pK2)
        
        # Water self-ionization: H2O <--> H+ + OH-
        pKw = 14.0 - 0.03 * (T - 25)
        self.constants['K_w'] = 10**(-pKw)
        
        # Calculate temperature-dependent Ksp for each polymorph
        for polymorph, params in self.polymorph_params.items():
            if params['enabled']:
                # Calculate log K based on temperature using parameters from Plummer & Busenberg
                A = params['eq_const_params']['A']
                B = params['eq_const_params']['B']
                C = params['eq_const_params']['C']
                D = params['eq_const_params']['D']
                
                log_K = A + B*T_K + C/T_K + D*np.log10(T_K)
                self.constants[f'K_sp_{polymorph}'] = 10**(log_K)
        
        # Temperature effects on kinetic parameters
        degas_temp_factor = np.exp(0.02 * (T - 25))
        precip_temp_factor = np.exp(0.035 * (T - 25))
        
        # Apply temperature factors to kinetic rates
        self.constants['k_degas_T'] = self.kinetic_params['degassing']['k_degas'] * degas_temp_factor
        self.constants['k_precip_T'] = self.kinetic_params['precipitation']['k_precip'] * precip_temp_factor
        
        # Ion association constants (CaHCO3+ and CaCO3(aq))
        # From Plummer & Busenberg data
        self.constants['K_CaHCO3'] = 10**(1.11)  # at 25°C
        self.constants['K_CaCO3aq'] = 10**(3.22)  # at 25°C
        
        # Print constants for reference
        print(f"Temperature set to {T}°C")
        print(f"K_H (Henry's constant) = {self.constants['K_H']:.6e}")
        print(f"K_a1 (First dissociation) = {self.constants['K_a1']:.6e}")
        print(f"K_a2 (Second dissociation) = {self.constants['K_a2']:.6e}")
        
        # Print Ksp values for enabled polymorphs
        for polymorph, params in self.polymorph_params.items():
            if params['enabled']:
                print(f"K_sp_{polymorph} = {self.constants[f'K_sp_{polymorph}']:.6e}")

    def set_polymorph(self, calcite=True, aragonite=False, vaterite=False):
        """
        Select which calcium carbonate polymorphs to enable in the simulation
        """
        self.polymorph_params['calcite']['enabled'] = calcite
        self.polymorph_params['aragonite']['enabled'] = aragonite
        self.polymorph_params['vaterite']['enabled'] = vaterite
        
        # Update equilibrium constants after changing polymorphs
        self.setup_equilibria()
        print(f"Enabled polymorphs: " + 
              (f"Calcite " if calcite else "") +
              (f"Aragonite " if aragonite else "") +
              (f"Vaterite" if vaterite else ""))

    def toggle_mechanism(self, mechanism, enabled=None):
        """
        Enable or disable specific mechanisms in the model
        
        Parameters:
        -----------
        mechanism : str
            The mechanism to toggle ('acid_catalysis', 'neutral_pathway', 'base_catalysis',
            'non_equilibrium', 'pH_lag')
        enabled : bool or None
            True to enable, False to disable, None to toggle
        """
        # Degassing mechanisms
        degassing_mechanisms = ['acid_catalysis', 'neutral_pathway', 'base_catalysis', 'non_equilibrium']
        if mechanism in degassing_mechanisms:
            current = self.kinetic_params['degassing'][mechanism]['enabled']
            new_value = not current if enabled is None else enabled
            self.kinetic_params['degassing'][mechanism]['enabled'] = new_value
            print(f"{mechanism} is now {'enabled' if new_value else 'disabled'}")
        
        # pH lag mechanism
        elif mechanism == 'pH_lag':
            current = self.kinetic_params['precipitation']['pH_lag']['enabled']
            new_value = not current if enabled is None else enabled
            self.kinetic_params['precipitation']['pH_lag']['enabled'] = new_value
            print(f"{mechanism} is now {'enabled' if new_value else 'disabled'}")
        
        # Master switches
        elif mechanism == 'degassing':
            current = self.kinetic_params['degassing']['enabled']
            new_value = not current if enabled is None else enabled
            self.kinetic_params['degassing']['enabled'] = new_value
            print(f"Degassing is now {'enabled' if new_value else 'disabled'}")
        
        elif mechanism == 'precipitation':
            current = self.kinetic_params['precipitation']['enabled']
            new_value = not current if enabled is None else enabled
            self.kinetic_params['precipitation']['enabled'] = new_value
            print(f"Precipitation is now {'enabled' if new_value else 'disabled'}")
        
        else:
            print(f"Unknown mechanism: {mechanism}")

    def set_parameters(self, **kwargs):
        """
        Set model parameters
        
        Parameters:
        -----------
        **kwargs : dict
            Parameters to set, using dot notation like 'degassing.k_degas' or 'precipitation.SI_critical'
        """
        parameters_updated = False
        temperature_updated = False
        
        for param_path, value in kwargs.items():
            path_parts = param_path.split('.')
            
            # Update system variables
            if path_parts[0] == 'system_vars' and len(path_parts) == 2:
                if path_parts[1] in self.system_vars:
                    old_value = self.system_vars[path_parts[1]]
                    self.system_vars[path_parts[1]] = value
                    print(f"Set system variable {path_parts[1]} = {value} (was {old_value})")
                    parameters_updated = True
                    
                    # Mark if temperature was updated
                    if path_parts[1] == 'Temperature':
                        temperature_updated = True
            
            # Update kinetic parameters
            elif path_parts[0] in ['degassing', 'precipitation'] and len(path_parts) >= 2:
                if len(path_parts) == 2:
                    # Direct parameter like 'degassing.k_degas'
                    if path_parts[1] in self.kinetic_params[path_parts[0]]:
                        old_value = self.kinetic_params[path_parts[0]][path_parts[1]]
                        self.kinetic_params[path_parts[0]][path_parts[1]] = value
                        print(f"Set {path_parts[0]}.{path_parts[1]} = {value} (was {old_value})")
                        parameters_updated = True
                
                elif len(path_parts) == 3:
                    # Nested parameter like 'degassing.acid_catalysis.k_acid'
                    if (path_parts[1] in self.kinetic_params[path_parts[0]] and 
                        path_parts[2] in self.kinetic_params[path_parts[0]][path_parts[1]]):
                        old_value = self.kinetic_params[path_parts[0]][path_parts[1]][path_parts[2]]
                        self.kinetic_params[path_parts[0]][path_parts[1]][path_parts[2]] = value
                        print(f"Set {path_parts[0]}.{path_parts[1]}.{path_parts[2]} = {value} (was {old_value})")
                        parameters_updated = True
            
            else:
                print(f"Warning: Unknown parameter {param_path}")
        
        # If temperature was updated, update equilibrium constants
        if temperature_updated:
            self.setup_equilibria()
        # If any degassing or precipitation parameter was updated, make sure constants are updated
        elif parameters_updated:
            # Update temperature-dependent kinetic parameters
            T = self.system_vars['Temperature']
            degas_temp_factor = np.exp(0.02 * (T - 25))
            precip_temp_factor = np.exp(0.035 * (T - 25))
            
            # Apply temperature factors to kinetic rates
            self.constants['k_degas_T'] = self.kinetic_params['degassing']['k_degas'] * degas_temp_factor
            self.constants['k_precip_T'] = self.kinetic_params['precipitation']['k_precip'] * precip_temp_factor
            
            print(f"Updated k_degas_T = {self.constants['k_degas_T']:.6e}")
            print(f"Updated k_precip_T = {self.constants['k_precip_T']:.6e}")

    def calculate_species_distribution(self, pH, DIC, Ca):
        """
        Calculate equilibrium distribution of carbonate species at given pH
        
        Parameters:
        -----------
        pH : float
            Solution pH
        DIC : float
            Total dissolved inorganic carbon (mol/L)
        Ca : float
            Calcium ion concentration (mol/L)
            
        Returns:
        --------
        dict
            Concentrations of all species
        """
        h_plus = 10**(-pH)
        oh_minus = self.constants['K_w'] / h_plus
        
        # Alpha factors (fraction of each carbonate species)
        denominator = (h_plus**2 + self.constants['K_a1']*h_plus + 
                      self.constants['K_a1']*self.constants['K_a2'])
        
        alpha0 = h_plus**2 / denominator
        alpha1 = self.constants['K_a1']*h_plus / denominator
        alpha2 = self.constants['K_a1']*self.constants['K_a2'] / denominator
        
        # Calculate species concentrations
        co2_aq = alpha0 * DIC
        hco3 = alpha1 * DIC
        co3 = alpha2 * DIC
        
        # Calculate saturation indices for each enabled polymorph
        SI = {}
        for polymorph, params in self.polymorph_params.items():
            if params['enabled']:
                SI[polymorph] = Ca * co3 / self.constants[f'K_sp_{polymorph}']
        
        return {
            'H': h_plus,
            'OH': oh_minus,
            'CO2': co2_aq,
            'HCO3': hco3,
            'CO3': co3,
            'SI': SI,
            'alpha0': alpha0
        }

    def calculate_pH(self, DIC, Ca):
        """
        Calculate pH from DIC and Ca based on charge balance
        
        Parameters:
        -----------
        DIC : float
            Total dissolved inorganic carbon (mol/L)
        Ca : float
            Calcium ion concentration (mol/L)
            
        Returns:
        --------
        float
            Calculated pH
        """
        def charge_balance(pH):
            h_plus = 10**(-pH)
            oh_minus = self.constants['K_w'] / h_plus
            
            # Alpha factors
            denominator = (h_plus**2 + self.constants['K_a1']*h_plus + 
                          self.constants['K_a1']*self.constants['K_a2'])
            
            alpha1 = self.constants['K_a1']*h_plus / denominator
            alpha2 = self.constants['K_a1']*self.constants['K_a2'] / denominator
            
            # Charge balance: [H+] + 2[Ca2+] = [HCO3-] + 2[CO32-] + [OH-]
            return h_plus + 2*Ca - alpha1*DIC - 2*alpha2*DIC - oh_minus
        
        # Initial guess at neutral pH
        pH_guess = 7.0
        
        try:
            pH = fsolve(charge_balance, pH_guess)[0]
            # Ensure reasonable range
            if pH < 0 or pH > 14:
                pH = 7.0
        except:
            pH = 7.0
        
        return pH

    def calculate_initial_conditions(self, initial_Ca, P_CO2):
        """
        Calculate initial equilibrium conditions for a given Ca and P_CO2
        
        Parameters:
        -----------
        initial_Ca : float
            Initial calcium concentration (mol/L)
        P_CO2 : float
            CO2 partial pressure (atm)
            
        Returns:
        --------
        dict
            Initial system state
        """
        def equations(vars):
            h_plus, ca, dic = vars
            
            # Calculate alpha factors
            denominator = (h_plus**2 + self.constants['K_a1']*h_plus + 
                          self.constants['K_a1']*self.constants['K_a2'])
            
            alpha0 = h_plus**2 / denominator
            alpha1 = self.constants['K_a1']*h_plus / denominator
            alpha2 = self.constants['K_a1']*self.constants['K_a2'] / denominator
            
            # CO2(aq) from Henry's Law
            co2_aq = self.constants['K_H'] * P_CO2
            
            # DIC based on CO2(aq) and pH
            dic_calc = co2_aq / alpha0 if alpha0 > 0 else co2_aq / 1e-6
            
            # CO3 concentration
            co3 = alpha2 * dic
            
            # Get Ksp for the active polymorph (use calcite as default)
            k_sp = None
            for polymorph, params in self.polymorph_params.items():
                if params['enabled']:
                    k_sp = self.constants[f'K_sp_{polymorph}']
                    break
            
            if k_sp is None:
                k_sp = self.constants['K_sp_calcite']  # Default to calcite
            
            # Ca concentration based on equilibrium
            if h_plus > 1.0e-7:  # pH < 7, dissolution dominates
                ca_calc = min(initial_Ca, k_sp / co3) if co3 > 0 else initial_Ca
            else:  # pH > 7, precipitation possible
                ca_calc = k_sp / co3 if co3 > 0 else initial_Ca
            
            # Charge balance
            charge_balance = h_plus + 2*ca - alpha1*dic - 2*alpha2*dic - self.constants['K_w']/h_plus
            
            return [dic - dic_calc, ca - ca_calc, charge_balance]
        
        # Initial guess: pH = 6, Ca = initial_Ca, DIC = CO2(aq) * 10
        co2_guess = self.constants['K_H'] * P_CO2
        initial_guess = [1.0e-6, initial_Ca, co2_guess * 10]
        
        try:
            solution = fsolve(equations, initial_guess)
            h_plus, ca, dic = solution
            
            # Calculate pH
            pH = -np.log10(h_plus)
            
            # Calculate species distribution
            species = self.calculate_species_distribution(pH, dic, ca)
            
            # Store initial state
            self.initial_pH = pH
            self.initial_Ca = ca
            self.initial_DIC = dic
            
            return {
                'pH': pH,
                'Ca': ca,
                'DIC': dic,
                'H': h_plus,
                'CO2': species['CO2'],
                'HCO3': species['HCO3'],
                'CO3': species['CO3'],
                'SI': species['SI']
            }
        except Exception as e:
            print(f"Error calculating initial equilibrium state: {e}")
            print("Using default values...")
            
            # Default values
            self.initial_pH = 6.0
            self.initial_Ca = initial_Ca
            self.initial_DIC = 0.01
            self.initial_H = 10**(-self.initial_pH)
            self.initial_CO2 = 0.005
            
            return {
                'pH': self.initial_pH,
                'Ca': self.initial_Ca,
                'DIC': self.initial_DIC,
                'H': self.initial_H,
                'CO2': self.initial_CO2,
                'HCO3': 0.005,
                'CO3': 0.0001,
                'SI': {'calcite': 1.0}
            }
        
    def degassing_precipitation_ode(self, t, y):
        """
        Coupled ODE system for CO2 degassing and CaCO3 precipitation
        
        State Variables (y):
        - y[0]: DIC [mol/L] - Total dissolved inorganic carbon
        - y[1]: Ca [mol/L] - Calcium ion concentration  
        - y[2]: pH_lag [-] - Lagged pH (accounts for kinetic pH response)
        
        Rate Processes:
        1. CO2 degassing: Removes CO2 from solution → decreases DIC
        2. CaCO3 precipitation: Removes Ca²⁺ and CO3²⁻ → decreases both Ca and DIC
        3. pH buffering: pH responds to changing DIC with time lag
        
        Returns:
        - [dDIC/dt, dCa/dt, dpH_lag/dt] [mol/(L·s), mol/(L·s), s⁻¹]
        """
        DIC, Ca, pH_lag = y
        
        # Calculate instantaneous pH
        pH_instant = self.calculate_pH(DIC, Ca)
        
        # Calculate pH lag rate
        if self.kinetic_params['precipitation']['pH_lag']['enabled']:
            dpH_lag_dt = (pH_instant - pH_lag) / self.kinetic_params['precipitation']['pH_lag']['tau']
        else:
            dpH_lag_dt = 0
            pH_lag = pH_instant
        
        # Calculate species using lagged pH
        species = self.calculate_species_distribution(pH_lag, DIC, Ca)
        co2_aq = species['CO2']
        
        # Degassing rate calculation
        if self.kinetic_params['degassing']['enabled']:
            # pH factor for degassing
            h_plus = 10**(-pH_lag)
            oh_minus = self.constants['K_w'] / h_plus
            
            pH_factor = 0
            if self.kinetic_params['degassing']['acid_catalysis']['enabled']:
                pH_factor += self.kinetic_params['degassing']['acid_catalysis']['k_acid'] * h_plus
            
            if self.kinetic_params['degassing']['neutral_pathway']['enabled']:
                pH_factor += self.kinetic_params['degassing']['neutral_pathway']['k_neutral']
            
            if self.kinetic_params['degassing']['base_catalysis']['enabled']:
                pH_factor += self.kinetic_params['degassing']['base_catalysis']['k_base'] * oh_minus
            
            # Non-equilibrium factor
            non_equil_factor = 1.0
            if self.kinetic_params['degassing']['non_equilibrium']['enabled']:
                alpha0 = species['alpha0']
                actual_fraction = co2_aq / DIC if DIC > 0 else alpha0
                deviation = abs(actual_fraction - alpha0)
                non_equil_factor += self.kinetic_params['degassing']['non_equilibrium']['k_factor'] * deviation
            
            # Calculate degassing rate - making sure to use the temperature-adjusted constant
            k_degas_T = self.constants['k_degas_T']
            
            # Debug print every 1000 seconds if needed
            if t % 10000 < 1:
                print(f"t={t:.1f}s, pH={pH_lag:.2f}, k_degas_T={k_degas_T:.6e}, pH_factor={pH_factor:.6e}")
                print(f"co2_aq={co2_aq:.6e}, K_H*P_atm={self.constants['K_H']*0.0004:.6e}")
                print(f"Degassing rate: {k_degas_T * (co2_aq - self.constants['K_H'] * 0.0004) * pH_factor * non_equil_factor:.6e}")
            
            r_degas = k_degas_T * (co2_aq - self.constants['K_H'] * 0.0004) * pH_factor * non_equil_factor
        else:
            r_degas = 0
        
        # Precipitation rate calculation remains the same
        if self.kinetic_params['precipitation']['enabled']:
            active_polymorph = None
            for polymorph, params in self.polymorph_params.items():
                if params['enabled']:
                    active_polymorph = polymorph
                    break
            
            if active_polymorph and active_polymorph in species['SI']:
                SI = species['SI'][active_polymorph]
                
                if SI > self.kinetic_params['precipitation']['SI_critical'] and t > self.kinetic_params['precipitation']['induction_time']:
                    # Calculate precipitation rate
                    r_precip = self.constants['k_precip_T'] * (SI - 1)**self.kinetic_params['precipitation']['n_reaction']
                    
                    # Limit precipitation at low Ca concentrations
                    if Ca < 0.001:  # 1 mmol/L threshold
                        r_precip *= (Ca / 0.001)**0.5
                else:
                    r_precip = 0
            else:
                r_precip = 0
        else:
            r_precip = 0
        
        # Calculate derivatives
        dDIC_dt = -r_degas - r_precip
        dCa_dt = -r_precip
        
        return [dDIC_dt, dCa_dt, dpH_lag_dt]

    def simulate(self, initial_Ca=0.01, P_CO2=1.0, time_span=(0, 150*3600), num_points=1000):
        """
        3. Simulate the system with user-defined parameters
        
        Parameters:
        -----------
        initial_Ca : float
            Initial calcium concentration (mol/L)
        P_CO2 : float
            Initial CO2 partial pressure (atm)
        time_span : tuple
            Simulation time range (start, end) in seconds
        num_points : int
            Number of output points
            
        Returns:
        --------
        dict
            Simulation results
        """
        # Calculate initial conditions
        initial_state = self.calculate_initial_conditions(initial_Ca, P_CO2)
        
        # Initial conditions for ODE solver [DIC, Ca, pH_lag]
        y0 = [initial_state['DIC'], initial_state['Ca'], initial_state['pH']]
        
        # Time points for output
        t_eval = np.linspace(time_span[0], time_span[1], num_points)
        
        # Solve ODE system
        try:
            solution = solve_ivp(
                self.degassing_precipitation_ode,
                time_span,
                y0,
                method='BDF',  # Backward Differentiation Formula for stiff problems
                t_eval=t_eval,
                rtol=1e-6,
                atol=1e-8
            )
            
            time = solution.t
            DIC_values = solution.y[0]
            Ca_values = solution.y[1]
            pH_lag_values = solution.y[2]
            
            # Calculate species for each time point
            CO2_values = np.zeros_like(time)
            HCO3_values = np.zeros_like(time)
            CO3_values = np.zeros_like(time)
            
            # Create dictionaries for SI values for each polymorph
            SI_values = {}
            for polymorph, params in self.polymorph_params.items():
                if params['enabled']:
                    SI_values[polymorph] = np.zeros_like(time)
            
            for i in range(len(time)):
                species = self.calculate_species_distribution(
                    pH_lag_values[i], DIC_values[i], Ca_values[i]
                )
                
                CO2_values[i] = species['CO2']
                HCO3_values[i] = species['HCO3']
                CO3_values[i] = species['CO3']
                
                for polymorph in SI_values:
                    if polymorph in species['SI']:
                        SI_values[polymorph][i] = species['SI'][polymorph]
            
            return {
                'time': time,
                'pH': pH_lag_values,
                'Ca': Ca_values,
                'DIC': DIC_values,
                'CO2': CO2_values,
                'HCO3': HCO3_values,
                'CO3': CO3_values,
                'SI': SI_values
            }
        except Exception as e:
            print(f"Simulation error: {str(e)}")
            # Return empty results on error
            dummy_array = np.zeros(num_points)
            empty_SI = {}
            for polymorph, params in self.polymorph_params.items():
                if params['enabled']:
                    empty_SI[polymorph] = dummy_array.copy()
                
            return {
                'time': np.linspace(time_span[0], time_span[1], num_points),
                'pH': dummy_array.copy(),
                'Ca': dummy_array.copy(),
                'DIC': dummy_array.copy(),
                'CO2': dummy_array.copy(),
                'HCO3': dummy_array.copy(),
                'CO3': dummy_array.copy(),
                'SI': empty_SI
            }

    def load_experimental_data(self, filepath, dataset_name=None, temperature=None):
        """
        4. Load experimental data for comparison
        
        Parameters:
        -----------
        filepath : str
            Path to CSV file containing experimental data
        dataset_name : str
            Name for the dataset (default: derived from filepath)
        temperature : float
            Temperature of the experimental data (°C)
            
        Returns:
        --------
        bool
            Success/failure
        """
        if not dataset_name:
            import os
            dataset_name = os.path.splitext(os.path.basename(filepath))[0]
            
        try:
            data = pd.read_csv(filepath)
            print(f"Successfully loaded data from {filepath}: {len(data)} points")
            
            # Check for required columns
            required_columns = ['Hours', 'Potential (pH)']
            missing_columns = [col for col in required_columns if col not in data.columns]
            
            if missing_columns:
                print(f"Warning: Missing required columns: {', '.join(missing_columns)}")
                return False
            
            # Store the data
            self.exp_data[dataset_name] = {
                'data': data,
                'temperature': temperature,
                'color': 'black',  # Default color
                'marker': 'o'      # Default marker
            }
            return True
            
        except Exception as e:
            print(f"Error loading experimental data: {str(e)}")
            return False
    
    def plot_results(self, results, title="Calcium Carbonate Precipitation Simulation", 
                    show_experimental=True, save_path=None, precip_stats=None):
        """
        Create plots of simulation results and experimental data
        
        Parameters:
        -----------
        results : dict
            Simulation results from simulate()
        title : str
            Plot title
        show_experimental : bool
            Whether to include experimental data in the plot
        save_path : str
            Path to save the figure, if None, only display
        precip_stats : dict
            Precipitation statistics to mark on the plot
        """
        # Create figure with subplots
        fig, axs = plt.subplots(2, 2, figsize=(14, 10))
        
        # Convert time to hours for plotting
        time_hours = results['time'] / 3600
        
        # Plot 1: pH vs Time
        axs[0, 0].plot(time_hours, results['pH'], 'k-', linewidth=2, label='Simulation')
        
        # Add precipitation start marker and average pH line if stats are provided
        if precip_stats and precip_stats['start_index'] is not None:
            start_hour = precip_stats['start_time_hours']
            start_pH = results['pH'][precip_stats['start_index']]
            
            # Mark precipitation start point
            axs[0, 0].scatter([start_hour], [start_pH], color='red', s=100, 
                            label=f'Precipitation Start (pH={start_pH:.2f})')
            axs[0, 0].axvline(x=start_hour, color='red', linestyle='--', alpha=0.5)
            
            # Add average pH line during precipitation
            axs[0, 0].axhline(y=precip_stats['average_pH'], color='green', linestyle='-', 
                            label=f'Avg pH={precip_stats["average_pH"]:.2f}')
        
        axs[0, 0].set_xlabel('Time (hours)')
        axs[0, 0].set_ylabel('pH')
        axs[0, 0].set_title('pH Evolution')
        axs[0, 0].grid(True, linestyle='--', alpha=0.7)
        
        # Add experimental pH data if available
        if show_experimental:
            for name, dataset in self.exp_data.items():
                if 'Hours' in dataset['data'].columns and 'Potential (pH)' in dataset['data'].columns:
                    axs[0, 0].plot(
                        dataset['data']['Hours'], 
                        dataset['data']['Potential (pH)'],
                        marker=dataset.get('marker', 'o'),
                        linestyle='',
                        color=dataset.get('color', 'blue'),
                        alpha=0.7,
                        label=f'Exp: {name}'
                    )
        
        # Plot 2: Calcium concentration vs Time
        axs[0, 1].plot(time_hours, results['Ca']*1000, 'g-', linewidth=2)  # Convert to mmol/L
        axs[0, 1].set_xlabel('Time (hours)')
        axs[0, 1].set_ylabel('Ca²⁺ Concentration (mmol/L)')
        axs[0, 1].set_title('Calcium Ion Evolution')
        axs[0, 1].grid(True, linestyle='--', alpha=0.7)
        
        # Add experimental calcium data if available
        if show_experimental:
            for name, dataset in self.exp_data.items():
                if 'Hours' in dataset['data'].columns and 'Calcium (mmol/L)' in dataset['data'].columns:
                    axs[0, 1].plot(
                        dataset['data']['Hours'], 
                        dataset['data']['Calcium (mmol/L)'],
                        marker=dataset.get('marker', 'o'),
                        linestyle='',
                        color=dataset.get('color', 'green'),
                        alpha=0.7,
                        label=f'Exp: {name}'
                    )
        
        # Plot 3: Carbonate Species vs Time
        # axs[1, 0].plot(time_hours, results['CO2']*1000, 'r-', linewidth=2, label='CO₂')
        axs[1, 0].plot(time_hours, results['HCO3']*1000, 'b-', linewidth=2, label='HCO₃⁻')
        axs[1, 0].plot(time_hours, results['CO3']*1000, 'g-', linewidth=2, label='CO₃²⁻')
        axs[1, 0].set_xlabel('Time (hours)')
        axs[1, 0].set_ylabel('Concentration (mmol/L)')
        axs[1, 0].set_title('Carbonate Species Evolution')
        axs[1, 0].grid(True, linestyle='--', alpha=0.7)
        axs[1, 0].legend()
        
        # Plot 4: Saturation Index vs Time
        for polymorph, si_values in results['SI'].items():
            params = self.polymorph_params[polymorph]
            axs[1, 1].plot(
                time_hours, 
                si_values, 
                linewidth=2, 
                color=params.get('color', 'blue'),
                label=f'{params["name"]} SI'
            )
        
        axs[1, 1].axhline(y=1.0, color='k', linestyle='--', alpha=0.7)
        axs[1, 1].text(time_hours[-1]*0.02, 1.05, 'Saturation (SI=1)', fontsize=9)
        
        # Add precipitation start marker on SI plot if stats are provided
        if precip_stats and precip_stats['start_index'] is not None:
            start_hour = precip_stats['start_time_hours']
            
            # Mark precipitation start point on SI plot
            axs[1, 1].axvline(x=start_hour, color='red', linestyle='--', alpha=0.5)
            
            # If we have SI values and start index is valid
            for polymorph, si_values in results['SI'].items():
                if precip_stats['start_index'] < len(si_values):
                    start_si = si_values[precip_stats['start_index']]
                    # Calculate average SI during precipitation
                    avg_si = np.mean(si_values[precip_stats['start_index']:])
                    
                    # Add average SI line
                    axs[1, 1].axhline(y=avg_si, color='green', linestyle='-', alpha=0.7)
                    axs[1, 1].text(time_hours[-1]*0.02, avg_si+0.05, 
                                f'Avg {polymorph} SI={avg_si:.2f}', fontsize=9, color='green')
        
        axs[1, 1].set_xlabel('Time (hours)')
        axs[1, 1].set_ylabel('Saturation Index (SI)')
        axs[1, 1].set_title('Saturation Index Evolution')
        axs[1, 1].grid(True, linestyle='--', alpha=0.7)
        axs[1, 1].legend()
        
        # Add horizontal line for critical SI if precipitation is enabled
        if self.kinetic_params['precipitation']['enabled']:
            si_crit = self.kinetic_params['precipitation']['SI_critical']
            axs[1, 1].axhline(y=si_crit, color='r', linestyle='--', alpha=0.7)
            axs[1, 1].text(time_hours[-1]*0.02, si_crit+0.05, f'Critical SI={si_crit}', fontsize=9, color='r')
        
        # Add legends where needed
        axs[0, 0].legend()
        if any('Calcium (mmol/L)' in dataset['data'].columns for dataset in self.exp_data.values()):
            axs[0, 1].legend()
        
        # Add super title
        plt.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # Save figure if requested
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {save_path}")
        
        plt.show()
    
    def calculate_average_pH_post_induction(self, results):
        """
        Calculate the average pH from the induction time (beginning of precipitation) until the end of simulation
        
        Parameters:
        -----------
        results : dict
            Simulation results from simulate()
            
        Returns:
        --------
        dict
            Dictionary containing average pH, start time index, and other relevant information
        """
        # Get induction time
        induction_time = self.kinetic_params['precipitation']['induction_time']
        
        # Find the first time index after the induction time
        time_array = results['time']
        induction_index = np.argmax(time_array >= induction_time)
        
        # If precipitation is not enabled or induction time is beyond simulation time, return None
        if (not self.kinetic_params['precipitation']['enabled'] or 
            induction_index >= len(time_array) or 
            induction_index == 0):
            return {
                'average_pH': None,
                'start_time_seconds': None,
                'start_time_hours': None,
                'start_index': None,
                'precipitation_enabled': self.kinetic_params['precipitation']['enabled'],
                'message': "No valid precipitation period found in the simulation"
            }
        
        # Extract pH values after induction time
        post_induction_pH = results['pH'][induction_index:]
        
        # Calculate average pH
        average_pH = np.mean(post_induction_pH)
        
        # Get start time in both seconds and hours
        start_time_seconds = time_array[induction_index]
        start_time_hours = start_time_seconds / 3600
        
        # Return results
        return {
            'average_pH': average_pH,
            'start_time_seconds': start_time_seconds,
            'start_time_hours': start_time_hours,
            'start_index': induction_index,
            'precipitation_enabled': self.kinetic_params['precipitation']['enabled'],
            'message': f"Average pH after induction time ({start_time_hours:.2f} hours): {average_pH:.4f}"
        }

    def get_precipitation_statistics(self, results):
        """
        Calculate comprehensive statistics about the precipitation phase
        
        Parameters:
        -----------
        results : dict
            Simulation results from simulate()
            
        Returns:
        --------
        dict
            Dictionary containing various precipitation statistics
        """
        # Get base information about precipitation phase
        precip_info = self.calculate_average_pH_post_induction(results)
        
        # If no valid precipitation period, return basic info
        if precip_info['average_pH'] is None:
            return precip_info
        
        # Get the start index for precipitation phase
        start_idx = precip_info['start_index']
        
        # Calculate calcium consumption (initial minus final)
        initial_Ca = results['Ca'][start_idx] 
        final_Ca = results['Ca'][-1]
        calcium_consumed = initial_Ca - final_Ca
        calcium_consumed_percent = (calcium_consumed / initial_Ca) * 100 if initial_Ca > 0 else 0
        
        # Calculate precipitation rate (average)
        time_elapsed = results['time'][-1] - results['time'][start_idx]
        average_precipitation_rate = calcium_consumed / (time_elapsed / 3600)  # mmol/L/hour
        
        # Calculate pH range during precipitation
        min_pH = np.min(results['pH'][start_idx:])
        max_pH = np.max(results['pH'][start_idx:])
        pH_range = max_pH - min_pH
        
        # Calculate maximum saturation index reached
        SI_values = list(results['SI'].values())[0]  # Get first polymorph's SI values
        max_SI = np.max(SI_values[start_idx:])
        
        # Return comprehensive statistics
        return {
            **precip_info,  # Include all information from the basic calculation
            'calcium_initial': initial_Ca,
            'calcium_final': final_Ca,
            'calcium_consumed': calcium_consumed,
            'calcium_consumed_percent': calcium_consumed_percent,
            'average_precipitation_rate': average_precipitation_rate,
            'min_pH': min_pH,
            'max_pH': max_pH,
            'pH_range': pH_range,
            'max_SI': max_SI,
            'message': (
                f"Precipitation Statistics (after {precip_info['start_time_hours']:.2f} hours):\n"
                f"- Average pH: {precip_info['average_pH']:.4f} (range: {min_pH:.4f} - {max_pH:.4f})\n"
                f"- Ca²⁺ consumed: {calcium_consumed*1000:.2f} mmol/L ({calcium_consumed_percent:.1f}%)\n"
                f"- Average precipitation rate: {average_precipitation_rate*1000:.4f} mmol/L/hour\n"
                f"- Maximum SI: {max_SI:.2f}"
            )
        }
    def plot_CO2_pH_relation(self, pH_range=(5.5, 8.5), DIC=0.01, Ca=0.01, title="CO₂ Concentration vs pH", save_path=None):
        """
        Plot the relationship between CO2 concentration and pH
        
        Parameters:
        -----------
        pH_range : tuple
            Range of pH values to plot
        DIC : float
            Dissolved inorganic carbon (mol/L)
        Ca : float
            Calcium concentration (mol/L)
        title : str
            Plot title
        save_path : str
            Path to save the figure, if None, only display
        """
        # Generate pH values
        pH_values = np.linspace(pH_range[0], pH_range[1], 100)
        
        # Calculate CO2 at each pH
        co2_values = np.zeros_like(pH_values)
        
        for i, pH in enumerate(pH_values):
            species = self.calculate_species_distribution(pH, DIC, Ca)
            co2_values[i] = species['CO2']
        
        # Convert to mg/L for plotting (CO2 molar mass = 44.01 g/mol)
        co2_mg_L = co2_values * 44.01 * 1000
        
        # Create figure
        plt.figure(figsize=(10, 6))
        plt.plot(pH_values, co2_mg_L, 'k-', linewidth=2.5)
        
        # Add exponential fit
        from scipy.optimize import curve_fit
        
        def exp_func(x, a, b):
            return a * np.exp(-b * x)
        
        try:
            popt, _ = curve_fit(exp_func, pH_values, co2_mg_L)
            plt.plot(pH_values, exp_func(pH_values, *popt), 'r--', 
                     label=f'y = {popt[0]:.1e}e^(-{popt[1]:.4f}x)')
            
            # Calculate R²
            residuals = co2_mg_L - exp_func(pH_values, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((co2_mg_L - np.mean(co2_mg_L))**2)
            r_squared = 1 - (ss_res / ss_tot)
            
            plt.text(0.05, 0.9, f'R² = {r_squared:.4f}', transform=plt.gca().transAxes)
        except:
            print("Could not fit exponential curve")
        
        plt.xlabel('pH')
        plt.ylabel('CO₂ Concentration (mg/L)')
        plt.title(title)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        
        # Save figure if requested
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {save_path}")
        
        plt.show()
    
    def customize_plot_appearance(self, dataset_name, color=None, marker=None):
        """
        Customize appearance of experimental data in plots
        
        Parameters:
        -----------
        dataset_name : str
            Name of the dataset to customize
        color : str
            Color for the dataset points
        marker : str
            Marker style for the dataset points
        """
        if dataset_name in self.exp_data:
            if color:
                self.exp_data[dataset_name]['color'] = color
            if marker:
                self.exp_data[dataset_name]['marker'] = marker
            print(f"Updated appearance for dataset '{dataset_name}'")
        else:
            print(f"Dataset '{dataset_name}' not found")
    
    def print_parameter_summary(self):
        """
        5. Print a summary of current model parameters
        """
        print("\n===== CaCO3 Precipitation Model Parameters =====")
        
        # Print enabled polymorphs
        print("\nEnabled Polymorphs:")
        for polymorph, params in self.polymorph_params.items():
            if params['enabled']:
                print(f"  - {params['name']}")
        
        # Print system variables
        print("\nSystem Variables:")
        for var, value in self.system_vars.items():
            if isinstance(value, float) and (value < 0.01 or value > 1000):
                print(f"  - {var}: {value:.2e}")
            else:
                print(f"  - {var}: {value}")
        
        # Print degassing parameters
        print("\nDegassing Parameters:")
        print(f"  - Enabled: {self.kinetic_params['degassing']['enabled']}")
        for key, value in self.kinetic_params['degassing'].items():
            if key != 'enabled' and not isinstance(value, dict):
                if isinstance(value, float) and (value < 0.01 or value > 1000):
                    print(f"  - {key}: {value:.2e}")
                else:
                    print(f"  - {key}: {value}")
        
        # Print sub-mechanisms
        print("\n  Sub-mechanisms:")
        for mech in ['acid_catalysis', 'neutral_pathway', 'base_catalysis', 'non_equilibrium']:
            print(f"  - {mech}: {'Enabled' if self.kinetic_params['degassing'][mech]['enabled'] else 'Disabled'}")
            for param, val in self.kinetic_params['degassing'][mech].items():
                if param != 'enabled':
                    if isinstance(val, float) and (val < 0.01 or val > 1000):
                        print(f"      {param}: {val:.2e}")
                    else:
                        print(f"      {param}: {val}")
        
        # Print precipitation parameters
        print("\nPrecipitation Parameters:")
        print(f"  - Enabled: {self.kinetic_params['precipitation']['enabled']}")
        for key, value in self.kinetic_params['precipitation'].items():
            if key != 'enabled' and not isinstance(value, dict):
                if isinstance(value, float) and (value < 0.01 or value > 1000):
                    print(f"  - {key}: {value:.2e}")
                else:
                    print(f"  - {key}: {value}")
        
        # Print pH lag info
        print(f"\n  pH Lag: {'Enabled' if self.kinetic_params['precipitation']['pH_lag']['enabled'] else 'Disabled'}")
        print(f"  - tau: {self.kinetic_params['precipitation']['pH_lag']['tau']}")
        
        # Print loaded experimental datasets
        if self.exp_data:
            print("\nLoaded Experimental Datasets:")
            for name, dataset in self.exp_data.items():
                points = len(dataset['data']) if 'data' in dataset else 0
                temp = dataset.get('temperature', 'Unknown')
                print(f"  - {name}: {points} points, Temperature: {temp}°C")


# Example usage
def main():
    """
    Example simulation demonstrating model capabilities
    
    This example shows:
    1. Model setup and parameter configuration
    2. Experimental data loading and comparison
    3. Simulation execution with realistic conditions
    4. Results analysis and visualization
    5. Statistical analysis of precipitation phase
    
    Modify parameters here to explore different scenarios:
    - Different temperatures (0-90°C)
    - Various initial Ca concentrations (0.001-0.1 mol/L)
    - Different CO2 pressures (atmospheric to high pressure)
    - Enable/disable different polymorphs and mechanisms
    """
    # Create model instance
    model = CaCO3PrecipitationModel()
    
    # Set polymorph to calcite
    model.set_polymorph(calcite=True, aragonite=False, vaterite=False)
    
    # Set parameters 
    model.set_parameters(**{
        'system_vars.Temperature': 25,
        'system_vars.P_CO2': 0.0001,
        'degassing.enabled': True,
        'degassing.k_degas': 1.03e-3,
        'degassing.acid_catalysis.enabled': True,
        'degassing.acid_catalysis.k_acid': 5.0e-8,
        'degassing.neutral_pathway.enabled': True,
        'degassing.neutral_pathway.k_neutral': 0.1,
        'degassing.base_catalysis.enabled': True,
        'degassing.base_catalysis.k_base': 5.0e-2,
        'degassing.non_equilibrium.enabled': True,
        'degassing.non_equilibrium.k_factor': 0.005,
        'precipitation.enabled': True,
        'precipitation.k_precip': 7e-9,
        'precipitation.n_reaction': 0.9,
        'precipitation.SI_critical': 1.0,
        'precipitation.induction_time': 45000,
        'precipitation.pH_lag.enabled': True,
        'precipitation.pH_lag.tau': 1
    })
    
    # Load experimental data (if available)
    model.load_experimental_data(
        filepath="25 No L.csv",
        dataset_name="Experimental Data",
        temperature=25
    )
    
    # Print parameter summary
    model.print_parameter_summary()
    
    ## Run simulation
    print("\nRunning simulation...")
    results = model.simulate(
        initial_Ca=0.7,
        P_CO2=1,
        time_span=(0, 200*3600),
        num_points=1000
    )
    
    # Calculate precipitation phase statistics
    precip_stats = model.get_precipitation_statistics(results)
    
    # Print precipitation statistics
    print("\n==== Precipitation Phase Statistics ====")
    print(precip_stats['message'])
    print(f"Precipitation start time: {precip_stats['start_time_hours']:.2f} 小时")
    
    # Output pH range during precipitation
    if precip_stats['start_index'] is not None:
        print(f"\npH range during precipitation:")
        print(f"- Minimum pH: {precip_stats['min_pH']:.4f}")
        print(f"- Maximum pH: {precip_stats['max_pH']:.4f}")
        print(f"- pH range: {precip_stats['pH_range']:.4f}")
        print(f"- Average pH: {precip_stats['average_pH']:.4f}")
    
   # Calculate and output saturation index range during precipitation
    if precip_stats['start_index'] is not None:
        start_idx = precip_stats['start_index']
        for polymorph, si_values in results['SI'].items():
            if model.polymorph_params[polymorph]['enabled']:
                # Calculate SI range during precipitation
                si_during_precip = si_values[start_idx:]
                min_si = np.min(si_during_precip)
                max_si = np.max(si_during_precip)
                avg_si = np.mean(si_during_precip)
                si_range = max_si - min_si
                
                print(f"\nSaturation index range during precipitation ({polymorph}):")
                print(f"- Minimum SI: {min_si:.4f}")
                print(f"- Maximum SI: {max_si:.4f}")
                print(f"- SI range: {si_range:.4f}")
                print(f"- Average SI: {avg_si:.4f}")
    # Calculate and output DIC range during precipitation
    if precip_stats['start_index'] is not None:
        start_idx = precip_stats['start_index']
        # Calculate DIC range during precipitation
        dic_during_precip = results['DIC'][start_idx:] * 1000  # 转换为mmol/L
        min_dic = np.min(dic_during_precip)
        max_dic = np.max(dic_during_precip)
        avg_dic = np.mean(dic_during_precip)
        dic_range = max_dic - min_dic
        
        # Calculate initial and final DIC values
        initial_dic = dic_during_precip[0]
        final_dic = dic_during_precip[-1]
        dic_change = final_dic - initial_dic
        dic_change_percent = (dic_change / initial_dic) * 100 if initial_dic > 0 else 0
        
        print(f"\n沉淀期间DIC范围:")
        print(f"- Minimum DIC: {min_dic:.4f} mmol/L")
        print(f"- Maximum DIC: {max_dic:.4f} mmol/L")
        print(f"- DIC range: {dic_range:.4f} mmol/L")
        print(f"- AverageDIC: {avg_dic:.4f} mmol/L")
        print(f"- Initial DIC: {initial_dic:.4f} mmol/L")
        print(f"- Final DIC: {final_dic:.4f} mmol/L")
        print(f"- DIC change: {dic_change:.4f} mmol/L ({dic_change_percent:.2f}%)")
    
    
    # Plot integrated results
    model.plot_results(
        results,
        title="CaCO3 Precipitation at 25°C",
        show_experimental=True,
        save_path="simulation_results.png",
        precip_stats=precip_stats  # Pass precipitation statistics to plotting function
    )
    
    print("\nSimulation completed.")

# Execute the main function if this script is run directly
if __name__ == "__main__":
    main()

#%%
