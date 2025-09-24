<h2 align="center">🚀 DIANO: Differentiable Autoencoding Neural Operator</h2>

**DIANO** is a framework for interpretable, physics-consistent latent-space modeling of high-dimensional spatio-temporal fields. It methodologically combines autoencoding, operator learning, and differentiable PDE solvers to enable efficient and interpretable simulations.

---

## ✨ Key Contributions

- 🧩 **Deterministic Autoencoding Neural Operator**  
  Mesh-invariant operator mapping of dimensionality reduction from high-dimensional fields to coarse latent representations and reconstructing them.

- 🔍 **Interpretable Coarse-Grid Latent Representation**  
  Coarsened latent space enabling visualization and interpretation of physics in the latent space.

- ⚡ **Latent-Space PDE Integration**  
  Differentiable PDE solvers embedded in the latent space for efficient, physically consistent temporal evolution.

- 🎯 **Flexible Solver-Accuracy Trade-offs**  
  Allows inexpensive lower-fidelity PDE solutions in the latent space with a trade-off between efficiency, accuracy, and interpretability.

- 📐 **Geometrical Reduction with Operator Learning**  
  Learns geometrical reductions (e.g., 2D → 1D), solves reduced PDEs, and maps back to high-dimensional fields.

---

## 🧪 Problem Scenarios

- 🖼️ **Nonlinear Dimensionality Reduction (Static Mapping)**  
  Encodes and decodes high-dimensional fields without temporal evolution:    
  ![Equation 1](https://latex.codecogs.com/png.latex?\mathbf{u}(t^n)%20\xrightarrow{\text{Encoder}}%20\mathbf{z}(t^n)%20\xrightarrow{\text{Decoder}}%20\hat{\mathbf{u}}(t^n))


- ⏩ **Nonlinear Dimensionality Reduction with Temporal Marching**  
  Latent evolution via differentiable PDE solver:  
  $$
  \mathbf{u}(t^n)
  \;\xrightarrow{\text{Encoder}}\;
  \mathbf{z}(t^n)
  \;\xrightarrow{\text{PDE Evolution}}\;
  \mathbf{z}(t^{n+1})
  \;\xrightarrow{\text{Decoder}}\;
  \hat{\mathbf{u}}(t^{n+1})
  $$

- 📏 **Geometrical Reduction with Temporal Marching**  
  Compresses geometric dimensions, evolves latent state, reconstructs original field:  
  $$
  \mathbf{u}_{D_h}(t^n)
  \;\xrightarrow{\text{Encoder}}\;
  \mathbf{z}_{D_\ell}(t^n)
  \;\xrightarrow{\text{PDE Evolution}}\;
  \mathbf{z}_{D_\ell}(t^{n+1})
  \;\xrightarrow{\text{Decoder}}\;
  \hat{\mathbf{u}}_{D_h}(t^{n+1})
  $$

- 🔗 **Many-to-One Functional Mapping via Latent Fusion**  
  Fuses multiple input fields in latent space for complex interactions:  
  $$
  \{\mathbf{u}^i(t^n)\}_{i=1}^m
  \;\xrightarrow{\text{Encoder}}\;
  \{\mathbf{z}^i(t^n)\}_{i=1}^m
  \;\xrightarrow{\text{PDE Mapping}}\;
  \mathbf{p}(t^n)
  \;\xrightarrow{\text{Decoder}}\;
  \hat{\mathbf{P}}(t^n)
  $$


---

