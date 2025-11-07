# Whipple Bicycle LQG Figures Generator

MATLAB code to reproduce all figures (2-8) from the paper:
> High-Rate Discrete-Time LQG Control of the Whipple Bicycle: Practical ZOH Modeling and Frequency-Domain Robustness

## Setup
- Requires MATLAB Control System Toolbox.
- Place Whipple model files (e.g., `build_whipple_ct_aug.m`) in the same directory.
- Run: `make_all_figs()`
- Outputs: `./figs/` with PDF/PNG figures.

## Usage
Edit parameters in the script (e.g., `Vgrid = [4 5 6];` for speeds).

## Citation
Cite the paper: [Link to paper PDF or DOI].

License: MIT
