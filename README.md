# Resonant Tunneling Through a Double-Barrier Structure (NEGF + Tight Binding)

This repository contains a MATLAB implementation of electron transport through a 1D double-barrier resonant tunneling structure using a tight-binding Hamiltonian and NEGF formalism.

## Structure

Lead | Barrier (10 nm, 25 eV) | Well (20 nm) | Barrier (10 nm, 25 eV) | Lead

Default parameters:
- Lattice spacing: `a = 0.25 nm`
- Barrier width: `Lb = 10 nm` (each)
- Well width: `Lw = 20 nm`
- Lead length (simulation padding): `Lc = 10 nm` (each side)
- Barrier height: `Vb = 25 eV`
- Effective mass: `m = m0`

## Outputs

1. Potential profile `V(x)`
2. Transmission spectrum `T(E)`
3. First resonant-state electron density `|psi1(x)|^2`
4. Second resonant-state electron density `|psi2(x)|^2`
5. Combined plot with potential, resonant energies, and densities

## Files

- `double_barrier_negf.m` — main simulation script

## How to run

1. Open MATLAB.
2. Navigate to the repository folder.
3. Run:

```matlab
run('double_barrier_negf.m')
```

The script will:
- print the hopping parameter and first two resonant energies,
- compute and plot transmission,
- plot resonant-state densities,
- plot potential profile and combined visualizations.

## Results

### 1) Combined: Potential Profile + Resonant-State Densities

![Combined plot](figures/figure1_combined.png)

### 2) Double Barrier Potential Profile

![Potential profile](figures/figure2_potential.png)

### 3) First and Second Resonant-State Densities

![Resonant state densities](figures/figure3_resonant_densities.png)

### 4) Transmission Spectrum

![Transmission spectrum](figures/figure4_transmission.png)

## Notes

- Transmission is computed from:
  \[
  T(E) = \mathrm{Tr}[\Gamma_1 G^r \Gamma_2 G^a]
  \]
- Resonant energies are estimated from the isolated well sub-Hamiltonian.
- Density plots use idealized sinusoidal bound-state shapes for report clarity.

If desired, this can be extended to compute densities directly from full-device scattering states / Green's functions at transmission peaks.
