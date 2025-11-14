# Delay-and-Sum Passive Acoustic Mapping (`das-pam`)

**MRes Neurotechnology Baseline Beamforming & k-Wave Simulation Module (March to April 2025)**  
Department of Bioengineering, Imperial College London

---

## Overview

**Passive Acoustic Mapping (PAM)** is a non-invasive ultrasound therapy monitoring technique used during Focused Ultrasound (FUS) interventions, such as targeted blood-brain barrier opening and thermal ablation. Unlike active B-mode ultrasound, which transmits pulses and measures backscattered echoes, PAM operates entirely in receive-only mode: a linear or 2D transducer array records the acoustic emissions radiated by oscillating or collapsing microbubbles (acoustic cavitation) and reconstructs a spatial map of cavitation power across the treatment domain.

This repository contains the foundational **2D k-Wave forward acoustic simulation** (`kWave_simulation.m`), **numerical grid convergence study** (`convergence_test.m`), and **Time-Domain Delay-and-Sum (DAS) beamformer** (`generate_pam.m`, `run_pam.m`) developed during the initial stage of my MRes Neurotechnology thesis (*Visualising Sound: A Physics-Guided Deep Learning Framework for Passive Acoustic Mapping*).

---

## Physical & Mathematical Formulation

### 1. 2D Acoustic Wave Propagation (`kWave_simulation.m`)

Acoustic wave propagation from a microbubble point source inside a homogeneous, lossless, isotropic soft-tissue medium is simulated using the **k-Wave** pseudo-spectral time-domain solver (`kspaceFirstOrder2D`). The solver integrates the coupled first-order conservation equations for momentum, mass, and the linearised equation of state:

$$\frac{\partial \mathbf{u}}{\partial t} = -\frac{1}{\rho_0} \nabla p, \qquad \frac{\partial \rho}{\partial t} = -\rho_0 \nabla \cdot \mathbf{u}, \qquad p = c_0^2 \rho$$

where $\mathbf{u}(\mathbf{r}, t)$ is the acoustic particle velocity, $p(\mathbf{r}, t)$ is the acoustic pressure, $\rho(\mathbf{r}, t)$ is the acoustic density perturbation, $\rho_0$ is the ambient medium density, and $c_0$ is the isentropic speed of sound. Combining the time derivative of the state equation with the divergence of the momentum equation yields the equivalent second-order wave equation driven by a point source $S(\mathbf{r}, t)$:

$$\nabla^2 p(\mathbf{r}, t) - \frac{1}{c_0^2} \frac{\partial^2 p(\mathbf{r}, t)}{\partial t^2} = -S(\mathbf{r}, t)$$

Fixed environment and discretisation parameters (`kWave_simulation.m`):

| Parameter | Symbol | Value |
| :--- | :--- | :--- |
| Computational domain | $L_x \times L_y$ | $100\text{ mm} \times 100\text{ mm}$ |
| Source Region of Interest (ROI) | $x \times y$ | $x \in [10, 90]\text{ mm},\; y \in [10, 80]\text{ mm}$ ($80\text{ mm} \times 70\text{ mm}$) |
| Medium sound speed | $c_0$ | $1500\text{ m/s}$ ($1.5\text{ mm/}\mu\text{s}$) |
| Medium ambient density | $\rho_0$ | $1040\text{ kg/m}^3$ |
| Source excitation | $f_c,\; N_c$ | Single point source, $1\text{ MHz}$ ($\lambda = 1.5\text{ mm}$), $5$-cycle sinusoidal tone burst |
| Production spatial resolution | $\text{PPW},\; N_x \times N_y$ | $15\text{ points per wavelength}$ ($1000 \times 1000$ grid, $\Delta x = \Delta y = 0.1\text{ mm}$) |
| Courant-Friedrichs-Lewy number | $\text{CFL}$ | $0.5$ ($\Delta t = \text{CFL} \cdot \Delta x / c_0 = 33.3\text{ ns}$) |
| Sensor array geometry | $M$ | $64$-element linear array along $x \in [10, 90]\text{ mm}$ at $y = 90\text{ mm}$ ($1.27\text{ mm}$ pitch) |
| Boundary absorption | $\text{PML}$ | $20$-grid-point Perfectly Matched Layer ($\alpha = 2\text{ Np/grid point}$) |

### 2. Numerical Discretisation & Convergence Study (`convergence_test.m`)

Spatial and temporal discretisation are coupled through the Courant-Friedrichs-Lewy ($\text{CFL}$) stability number:

$$\text{CFL} = \frac{c_0 \Delta t}{\Delta x}$$

To determine the minimum spatial resolution required to suppress numerical dispersion and preserve the full wavefront bandwidth across the $100\text{ mm} \times 100\text{ mm}$ domain, `convergence_test.m` sweeps grid dimensions from $N_x = 100$ ($1.5\text{ PPW}$) to $N_x = 1000$ ($15\text{ PPW}$) across $\text{CFL} \in \{0.3, 0.5\}$ and evaluates four metrics:
1. **Peak Positive Pressure:** Maximum recorded pressure amplitude $\max_{m, n} p_m(n\Delta t)$ across the 64-element array.
2. **Time-Integrated Total Energy:** Discrete acoustic energy $\sum_{m=1}^{M} \sum_{n=0}^{N_t - 1} p_m(n\Delta t)^2 \, \Delta t$ integrated across all channels.
3. **Waveform Similarity:** Pearson correlation coefficient $r$ between the central sensor waveform ($m = 32$) and the high-resolution reference waveform ($N_x = 1000$, $15\text{ PPW}$) interpolated onto a common physical time axis.
4. **Computational Runtime:** Wall-clock execution time per simulation as a function of grid size $N_x$ and $\text{CFL}$.

Summary of convergence regimes across the spatial resolution sweep:

| Grid Size ($N_x \times N_y$) | Grid Spacing ($\Delta x$) | Points Per Wavelength ($\text{PPW}$) | Time Step at $\text{CFL}=0.5$ ($\Delta t$) | Numerical Convergence Behaviour |
| :---: | :---: | :---: | :---: | :--- |
| $100 \times 100$ | $1.00\text{ mm}$ | $1.5$ | $333.3\text{ ns}$ | Under-resolved below Nyquist limit ($2\text{ PPW}$); severe phase dispersion and unstable amplitude/energy |
| $200 \times 200$ | $0.50\text{ mm}$ | $3.0$ | $166.7\text{ ns}$ | Waveform shape correlation exceeds $r > 0.99$ relative to $N_x = 1000$, but peak pressure and total energy remain transitional |
| $400 \times 400$ | $0.25\text{ mm}$ | $6.0$ | $83.3\text{ ns}$ | Onset of numerical convergence; peak positive pressure and time-integrated total energy plateau ($r > 0.99$) |
| $667 \times 667$ | $0.15\text{ mm}$ | $10.0$ | $50.0\text{ ns}$ | Fully converged intermediate resolution suitable for rapid desktop DAS prototyping (`run_pam.m`) |
| **$1000 \times 1000$** | **$0.10\text{ mm}$** | **$15.0$** | **$33.3\text{ ns}$** | **Production high-fidelity regime:** eliminates residual numerical dispersion for training the 10,000-simulation deep learning dataset |

Comparing $\text{CFL} = 0.3$ and $\text{CFL} = 0.5$ showed identical pressure, energy, and waveform correlation curves for $N_x \ge 400$ while reducing computational runtime by approximately $40\%$ at $\text{CFL} = 0.5$. Consequently, $N_x = 1000$ ($15\text{ PPW}$) and $\text{CFL} = 0.5$ were adopted for full dataset generation.

### 3. Delay-and-Sum (DAS) Beamforming (`generate_pam.m`)

Given recorded sensor time series $s_m(t)$ across $M = 64$ receivers at coordinates $\mathbf{r}_m = (x_m, y_m)$, each candidate pixel $\mathbf{r}_p = (x_p, y_p)$ on the reconstruction grid is focussed by computing the geometric time-of-flight propagation delay:

$$\tau_m(\mathbf{r}_p) = \frac{\|\mathbf{r}_p - \mathbf{r}_m\|_2}{c_0} = \frac{\sqrt{(x_p - x_m)^2 + (y_p - y_m)^2}}{c_0}$$

Each channel signal is advanced by $\tau_m(\mathbf{r}_p)$ via linear interpolation (`interp1`) to align wavefront arrivals originating from $\mathbf{r}_p$, coherently summed across the array aperture, and integrated over the acquisition window $T$ (Time Exposure Acoustics, [Gyöngy & Coussios, 2010](https://doi.org/10.1109/TBME.2009.2026907)) to estimate the spatial cavitation power $I_{\text{DAS}}(\mathbf{r}_p)$:

$$q(\mathbf{r}_p, t) = \sum_{m=1}^{M} w_m \, s_m\!\left(t + \tau_m(\mathbf{r}_p)\right), \qquad I_{\text{DAS}}(\mathbf{r}_p) = \sum_{n=0}^{N_t - 1} \left|q(\mathbf{r}_p, n\Delta t)\right|^2$$

with uniform sensor weighting $w_m = 1$.

---

## Classical Limitations & Progression to Physics-Guided Deep Learning

While Delay-and-Sum beamforming provides a transparent, closed-form physical baseline for localising acoustic cavitation, the one-sided $80\text{ mm}$ linear array aperture and diffraction limits introduce characteristic spatial artefacts:

1. **Axial Elongation ("Tail" Artefact):** Point spread functions are tightly focussed laterally ($1.29\text{ mm}$ median lateral FWHM at $0\text{ dB}$ SNR) but heavily elongated along the propagation axis ($3.15\text{ mm}$ median axial FWHM), degrading depth resolution.
2. **Cross-Shaped Sidelobes & Low-SNR Breakdown:** Coherent summation without adaptive spatial filtering produces prominent off-axis sidelobes ($\text{PAR} = -0.19\text{ dB}$) and fails below $-15\text{ dB}$ SNR ($0\%$ valid PSF at $-25\text{ dB}$ SNR).
3. **Lossy Summation Bottleneck:** Collapsing the aligned multi-channel spatio-temporal wavefront into a single scalar energy value per pixel discards the fine inter-channel coherence structure needed to distinguish true cavitation sources from off-axis constructive interference.

In the companion MRes repositories, this foundation is extended to:
- **Delay-and-Classify (`ResDAC`):** Retaining the physics-based time-of-flight alignment step while replacing the lossy summation with a ResNet-18 binary classifier trained with hard negative sampling, weighted cross-entropy, and a correlated noise curriculum (see [`ResNet_Beamformer`](https://github.com/uam225/ResNet_Beamformer)).
- **Physics-Informed Neural Networks (`Wave_PINN`):** Enforcing the 2D acoustic wave PDE residual ($\nabla^2 p - c_0^{-2}\partial_{tt}p = 0$), quiescent initial velocity conditions, and k-Wave PML boundary decay directly in a Fourier-feature neural network loss (see [`Wave_PINN`](https://github.com/uam225/Wave_PINN)).

---

## Repository Structure

```text
├── kWave_simulation.m   # 2D k-Wave forward simulation (1000x1000 grid, 15 PPW, 64-element array, 1 MHz burst)
├── convergence_test.m   # Grid resolution (Nx = 100..1000, 1.5..15 PPW) and CFL (0.3, 0.5) convergence study
├── generate_pam.m       # Time-domain Delay-and-Sum (DAS) Passive Acoustic Mapping beamformer function
├── run_pam.m            # Reconstruction driver script and 2D cavitation power map visualisation
├── .gitignore           # Excludes large .mat simulation dumps and MATLAB autosave files
└── LICENSE              # MIT Licence
```

---

## Usage

### Prerequisites

- **MATLAB** (R2023a or newer recommended)
- [**k-Wave MATLAB Toolbox**](http://www.k-wave.org/) added to your MATLAB path

### 1. Run Grid Discretisation Convergence Study

To reproduce the spatial and temporal discretisation sweep across $N_x \in [100, 1000]$ ($1.5$ to $15\text{ PPW}$) and $\text{CFL} \in \{0.3, 0.5\}$:

```matlab
convergence_test
```

This prints the convergence table (`Nx`, `PPW`, `Max_Pressure`, `Total_Energy`, `Waveform_Correlation`, `Runtime_s`) and plots the four-panel convergence figure.

### 2. Run Forward k-Wave Simulation

Open `kWave_simulation.m` in MATLAB, uncomment the `save(...)` line at the bottom of the script to export `sensor_data` to a `.mat` file, and run:

```matlab
kWave_simulation
```

### 3. Reconstruct Passive Acoustic Map

Ensure the `load(...)` path in `run_pam.m` points to your generated `.mat` simulation output, then execute:

```matlab
run_pam
```

This computes `generate_pam(sensor_data, sensor_positions, grid_x, grid_y, c, dt)` and plots the 2D spatial cavitation power map ($X, Y \in [0, 100]\text{ mm}$).

---

## Licence

Released under the [MIT Licence](LICENSE).
