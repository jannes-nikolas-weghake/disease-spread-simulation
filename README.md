# Disease Spread Simulation

## Overview
This project implements computational models for simulating disease spread using SIS (Susceptible-Infected-Susceptible) and SIR (Susceptible-Infected-Recovered) models. It features both basic simulations and advanced multi-city scenarios with factors like mobility, vaccination, and mortality.

## Features
- Implementation of SIS and SIR epidemiological models
- Numerical simulations using Euler and 4th order Runge-Kutta methods
- Comparative analysis of algorithm accuracy and performance
- Advanced multi-city simulation incorporating:
  - Inter-city mobility
  - Vaccination programs
  - Mortality rates
- Visualization of simulation results

## Project Structure
- `SIS_model.py`: Implementation of the SIS model
- `SIR_model.py`: Implementation of the SIR model
- `euler_method.py`: Euler method for numerical integration
- `runge_kutta_method.py`: 4th order Runge-Kutta method
- `error_analysis.py`: Comparison of numerical errors between methods
- `multi_city_simulation.py`: Advanced simulation with multiple cities
- `visualization.py`: Functions for plotting results

## Requirements
- Python 3.x
- NumPy
- Matplotlib

## Results


The project includes various analyses:

Comparison of SIS and SIR model dynamics
Error analysis between Euler and Runge-Kutta methods
Performance benchmarks of numerical methods
Visualization of disease spread in multi-city scenarios

## Future Work

Implement stochastic variations of the models
Develop a user-friendly interface for parameter adjustment
