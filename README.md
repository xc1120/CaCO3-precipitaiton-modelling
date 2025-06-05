# CaCO3-precipitaiton-modelling

A comprehensive Python model for simulating CaCO₃ precipitation in CO₂-H₂O systems with advanced thermodynamic and kinetic capabilities.

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)

## 🌟 Overview

This model simulates the coupled processes of CO₂ degassing and calcium carbonate precipitation in aqueous systems. It's designed for researchers and engineers working in geochemistry, materials science, environmental engineering, and related fields.

### Key Features

- **Multiple CaCO₃ Polymorphs**: Calcite, aragonite, and vaterite with temperature-dependent solubilities
- **Advanced Kinetics**: pH-dependent CO₂ degassing pathways and precipitation kinetics
- **Temperature Effects**: Full temperature dependence (0-90°C) for all equilibrium and kinetic parameters
- **Experimental Integration**: Load and compare with experimental data
- **Comprehensive Visualization**: Multi-panel plots with statistical analysis
- **Flexible Parameter Control**: Easy modification of all model parameters

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/caco3-precipitation-model.git
cd caco3-precipitation-model

```

### Basic Usage

```python
from caco3_model import CaCO3PrecipitationModel

# Create model instance
model = CaCO3PrecipitationModel()

# Configure for calcite precipitation
model.set_polymorph(calcite=True)

# Run simulation
results = model.simulate(
    initial_Ca=0.01,    # 10 mM calcium
    P_CO2=1.0,          # 1 atm CO2 pressure
    time_span=(0, 48*3600),  # 48 hours
    num_points=1000
)

# Visualize results
model.plot_results(results)

# Get precipitation statistics
stats = model.get_precipitation_statistics(results)
print(stats['message'])
```

## 🔬 Scientific Background

### Chemical System

The model solves the coupled equilibria:

1. **CO₂ Solubility**: CO₂(g) ⇌ CO₂(aq)
2. **Carbonic Acid**: CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻
3. **Bicarbonate**: HCO₃⁻ ⇌ H⁺ + CO₃²⁻
4. **Precipitation**: Ca²⁺ + CO₃²⁻ ⇌ CaCO₃(s)
5. **Water**: H₂O ⇌ H⁺ + OH⁻

### Kinetic Processes

- **CO₂ Degassing**: Multiple pH-dependent pathways (acid, neutral, base catalysis)
- **Precipitation**: Nucleation and growth with saturation index control
- **pH Buffering**: Kinetic pH response with adjustable time constants

### Temperature Dependence

All equilibrium constants and kinetic parameters include empirical temperature corrections based on literature data (Plummer & Busenberg, 1982).

## 📋 Parameters Guide

### System Variables
| Parameter | Unit | Description | 
|-----------|------|-------------|
| `pH` | [-] | Solution pH | 0-14 |
| `Ca` | [mol/L] | Calcium concentration |
| `DIC` | [mol/L] | Dissolved inorganic carbon | 
| `Temperature` | [°C] | Solution temperature |
| `P_CO2` | [atm] | CO₂ partial pressure |

### Key Kinetic Parameters
| Parameter | Unit | Description | Typical Range |
|-----------|------|-------------|---------------|
| `degassing.k_degas` | [s⁻¹] | Base degassing rate | 1e-5 to 1e-2 |
| `precipitation.k_precip` | [L/(mol·s)] | Precipitation rate constant | 1e-10 to 1e-6 |
| `precipitation.SI_critical` | [-] | Critical saturation index | 1.0-2.0 |
| `precipitation.induction_time` | [s] | Nucleation induction time | 0-100000 |

## 📈 Example Results

The model produces comprehensive output including:

- **pH Evolution**: Time-dependent pH changes during degassing and precipitation
- **Species Distribution**: CO₂, HCO₃⁻, CO₃²⁻ concentrations over time
- **Saturation State**: Saturation indices for all enabled polymorphs
- **Precipitation Statistics**: Rates, consumed Ca²⁺, pH ranges during precipitation

## 🧪 Experimental Data Integration

Load experimental data for model validation:

```python
# Load experimental data
model.load_experimental_data(
    filepath="experiment_25C.csv",
    dataset_name="Lab_Experiment_1",
    temperature=25
)

# Customize appearance
model.customize_plot_appearance(
    dataset_name="Lab_Experiment_1",
    color='red',
    marker='s'
)

# Plot with experimental comparison
model.plot_results(results, show_experimental=True)
```

### Required CSV Format
```csv
Hours,Potential (pH),Calcium (mmol/L)
0.0,6.85,10.2
0.5,7.12,10.1
1.0,7.34,9.8
...
```

## The example provided in the paper

The default setting in the code and the document'25 No L.csv' (25°C, no ion impurities, active degassing method) was used to produce Fig. S3 in the paper. The resulting figure shows the (1) comparsion of pH vs time between the real experimental data; (2) calcium evolution; (3) DIC speciation evolution;(4) saturation index evolution.
![image](https://github.com/user-attachments/assets/868f1f40-b087-453a-9621-9a6451540201)


## 🙏 Acknowledgments

- Thermodynamic data based on Plummer & Busenberg (1982)
- Kinetic formulations adapted from classical precipitation theory
- Temperature corrections from multiple literature sources

---

**Keywords**: calcium carbonate, precipitation, geochemistry, CO2 sequestration, crystallization, Python modeling
