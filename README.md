# Compressible Interface Solver

This MATLAB code models the evolution of a gas-liquid interface and the associated pressure field during the injection of a compressible gas into a confined porous layer. The interface is bounded by two dynamic contact lines: Xₗ at the lower boundary and Xᵤ at the upper boundary.

- **Region I (0 < x < Xₗ):** Here, the pressure field `P` is solved for directly based on the interface shape.
- **Region II (Xₗ < x < Xᵤ):** The equations are rewritten in terms of `B` and `F`, which are evolved in time.

When the lower contact line reaches the origin (i.e., `Xₗ = 0`), only Region II remains, and both `B` and `F` are solved across the full domain.

The discretization uses a second-order accurate method-of-lines (MOL) scheme with central differencing in space. Time integration is handled using MATLAB's ODE solvers (e.g., `ode15s`), making use of a custom Jacobian sparsity pattern for efficiency.

This code was developed to accompany the paper:  
**"Dynamic injection of a compressible gas into a confined porous layer"**  
*Peter Castellucci, Radha Boya, Lin Ma, Igor Chernyavsky, Oliver Jensen. (2025)*  

## 🛠️ Requirements

- **MATLAB R2024b** or later
- Tested on **macOS Ventura**
- No third-party toolboxes required (except core functions like `numjac`)

---

---

## ▶️ How to Run

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/compressible-interface-solver.git
    cd compressible-interface-solver/code
    ```

2. Open MATLAB in this directory.

3. To simulate the system, run the main script:

    ```matlab
    MOL
    ```

This will simulate the system using default parameters and save the data to a file called "Data.txt".

🔧 You can modify key physical parameters (e.g. `M`, `zeta`, `L`) as well as the time dependence of the source `Q(t)` directly in `MOL.m` to explore different dynamics.
