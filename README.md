<h2 align="center">🚀 DIANO: DIfferentiable Autoencoding Neural Operator</h2>

<p align="center">
  <a href="https://www.arxiv.org/abs/2510.00233"><img src="https://img.shields.io/badge/arXiv-2510.00233-b31b1b.svg" alt="arXiv"></a>
  <a href="https://github.com/siva-viknesh/DIfferentiable_Autoencoding_Neural_Operator"><img src="https://img.shields.io/badge/GitHub-DIANO-blue?logo=github" alt="GitHub"></a>
  <img src="https://img.shields.io/badge/PyTorch-Framework-orange?logo=pytorch" alt="PyTorch">
  <img src="https://img.shields.io/badge/License-MIT-green" alt="License">
</p>

**DIANO** is a deterministic, physics-integrated framework for interpretable latent-space modeling of high-dimensional spatio-temporal flow fields. It unifies autoencoding, neural operator learning, and differentiable PDE solvers into a single end-to-end differentiable ML architecture — enabling compact, physically meaningful latent representations and robust generalization across different Reynolds numbers, without sacrificing predictive accuracy.

<p align="center">
  <img src="https://raw.githubusercontent.com/siva-viknesh/DIfferentiable_Autoencoding_Neural_Operator/main/0_Dimension_Reduction_Static_Mapping/DIANO.png" 
       alt="DIANO Framework Overview" 
       width="900"/>
</p>

## ✨ Key Contributions

- 🧩 **Deterministic Autoencoding Neural Operator**  
  Mesh-invariant operator framework that maps high-dimensional flow fields (N×N) to physically interpretable coarse-grid latent representations (M×M, N > M) and reconstructs them via Fourier-basis encoder/decoder operators.

- 🔍 **Interpretable Coarse-Grid Latent Representation**  
  The latent space is defined as a coarsened physical grid, enabling direct visualization and physical interpretation of flow structures (e.g., vortex shedding, pressure fields) without additional post-processing.

- ⚡ **Latent-Space PDE Integration**  
  Fully differentiable PDE solvers (finite difference, RK4, Point-Jacobi) are embedded directly in the latent space, ensuring temporal evolution is governed by physics rather than black-box ML approximations.

- 🎯 **Flexible Solver-Accuracy Trade-offs**  
  Supports a spectrum of PDE fidelities — from full nonlinear VTE to Stokes flow and 1D simplified variants — enabling principled trade-offs between computational efficiency, reconstruction accuracy, and latent interpretability.

- 📐 **Geometrical Reduction with Operator Learning**  
  Learns data-driven geometrical reductions (e.g., 2D → 1D), solves a geometrically reduced PDE in the latent space, and decodes back to the original high-dimensional geometry.

---

## 🏗️ Framework Overview

DIANO consists of three tightly coupled components:

```
High-dim Input (N×N)
        │
   ┌────▼────────────────────┐
   │   Encoder Neural        │   Fourier layers + AvgPool2D downsampling
   │   Operator              │   → coarse-grid latent space (M×M)
   └────────────────────┬────┘
                        │
             ┌──────────▼──────────┐
             │  Differentiable     │   FDM (OUCS2 upwind + central diff)
             │  PDE Solver         │   RK4 (parabolic) / Point-Jacobi (elliptic)
             └──────────┬──────────┘
                        │
   ┌────────────────────▼────┐
   │   Decoder Neural        │   ConvTranspose2D upsampling + Fourier layers
   │   Operator              │   → reconstructed high-dim output (N×N)
   └─────────────────────────┘
```

---

## 🔬 Modeling Scenarios

DIANO is evaluated across four representative configurations:

| Scenario | Description | PDE Solver |
|---|---|---|
| **Static Mapping** | Nonlinear dimensionality reduction at a single time instant | — |
| **Temporal Marching** | Latent-space time advancement via PDE solver | 2D VTE variants |
| **Geometrical Reduction** | 2D → 1D geometric compression + temporal marching | 1D linearized VTE |
| **Many-to-One Mapping** | Multiple inputs (u, v, w velocity) → single output (pressure) | 3D Pressure Poisson |

### PDE Formulations Supported

**2D Vorticity Transport Equation (VTE)** and variants:
- Full 2D linearized VTE
- 2D Stokes flow (no convection)
- 2D inviscid linearized VTE
- 1D linearized VTE (streamwise or normal direction)

**3D Pressure Poisson Equation (PPE):**
- Solved iteratively via Point-Jacobi on the coarse latent grid

---

## 📊 Benchmark Problems

| Case | Description | Scenario |
|---|---|---|
| **Flow past a 2D cylinder** | Re = 100–225, periodic vortex shedding (von Kármán street) | Static mapping, temporal marching |
| **2D symmetric stenosed artery** | Pulsatile flow, 50% blockage, post-stenotic recirculation | Geometrical reduction |
| **3D patient-specific LAD artery** | Patient-specific coronary artery, stenosed, 3D velocity → pressure | Many-to-one mapping |

---

## 📈 Results Highlights

- DIANO achieves reconstruction MSE on the order of **10⁻⁶–10⁻⁷** for static mapping, competitive with or better than CNN-AE, NN-AE, and CNO baselines, using **~1.2M parameters** vs. 31M for CNO.
- Latent vorticity fields faithfully reproduce the **k⁻³ energy cascade** spectrum of the ground-truth vorticity field.
- Autoregressive rollout errors **saturate at O(10⁻²)** across seen, interpolated (Re = 180), and extrapolated (Re = 225) Reynolds numbers — no unbounded error growth.
- DIANO outperforms PPNN and LaSDI on long-horizon rollout stability while maintaining interpretable latent dynamics.
- Robust to input noise up to **15%** (Gaussian), with graceful degradation at 25%.

---

## 🔭 Future Directions

- **Unstructured grid support** via geometry-aware neural operators
- **Turbulent flow modeling** and multiscale latent representations
- **Latent PDE discovery** — treating PDE parameters as trainable quantities
- **Attention-based autoencoding operators** for long-range spatial dependencies

---

## 📄 Citation

If you use DIANO in your research, please cite:

```bibtex
@article{viknesh2026diano,
  title   = {Differentiable Autoencoding Neural Operator for Interpretable 
             and Integrable Latent Space Modeling},
  author  = {Viknesh, Siva and Arzani, Amirhossein},
  journal = {arXiv preprint arXiv:2510.00233},
  year    = {2026}
}
```

---
