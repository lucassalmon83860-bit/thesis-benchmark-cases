# thesis-benchmark-cases

# Abaqus Cohesive Zone Benchmark Library

**Abaqus benchmark cases and UMAT cohesive zone models for thermo-mechanical crack propagation and healing. Input files and subroutines are provided, with reference results reported in the associated publication (https://doi.org/10.1007/s00466-025-02675-5) for reproducibility.**

This repository provides benchmark cases developed during the Ph.D. thesis:

**“Multiphysics modeling of crack propagation and healing in a thermo-elastic medium: application to nuclear fuel behavior”**  

It focuses on **thermo-mechanical crack propagation and healing** using **cohesive zone models implemented in Abaqus**. The goal is to provide tools enabling reproducible simulations of damage–healing cycles in cohesive interfaces.

## Content

- Automated procedures for inserting cohesive zones
- Abaqus input files for benchmark cases with associated UMAT subroutines (Fortran 90) implementing thermo-mechanical cohesive laws


## Benchmarks

The repository currently includes:

| Name | Description | Features |
|-----|-------------|----------|
| Notched plate (tension–compression) | Three cycles of tensile–compressive loading on a notched plate | Crack opening/closure and partial stiffness recovery |
| Notched plate (shear) | Three cycles of shear loading on a notched plate | Crack opening/closure and partial stiffness recovery |
| Fuel pellet quarter | Representative irradiation cycle on a quarter of a nuclear fuel pellet | Crack initiation, propagation, and healing under thermo-mechanical conditions |



