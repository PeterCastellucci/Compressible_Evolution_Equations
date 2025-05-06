# Compressible Interface Dynamics Solver

This repository contains MATLAB code for simulating the evolution of a compressible interface in a two-region domain, as presented in our paper:

**[Paper Title]**  
*Author1, Author2, Author3*  
*Journal of Fluid Mechanics*, 2025 (forthcoming)  
[DOI if available]

---

## 📌 Overview

The solver implements a second-order accurate **method-of-lines (MOL)** scheme to simulate the evolution of:
- **F(x, t)**: Interface height
- **B(x, t)**: Depth-integrated density, where B = (1 − F)·P

The domain is split into **two regions** depending on the position of the lower contact line. The script handles both **two-region** and **single-region** cases automatically, based on the contact line dynamics.

---

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

3. Run the main script:
    ```matlab
    run_main
    ```

This will simulate the system with default parameters and generate outputs including interface evolution and sample plots.

---

## 📣 Citation

If you use this code, please cite the following paper:

> Author1, Author2, Author3. *Title*. Journal of Fluid Mechanics, 2025. [DOI]

BibTeX:
```bibtex
@article{author2025jfm,
  title={Title},
  author={Author1 and Author2 and Author3},
  journal={Journal of Fluid Mechanics},
  year={2025}
}
