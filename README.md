# Methodology

## 1. Modelling Approach

The objective of this study is to compute the flow distribution through the terminal branches of an arterial tree and to quantify the total pressure required drop to drive flow through the system. Blood is modelled as a steady, incompressible, Newtonian fluid ($`\rho = 1 \times 10^{-3}`$ g/mm³, $`\mu = 4 \times 10^{-3}`$ g/(mm·s)) flowing through a two-dimensional cross-sectional representation. This approach is appropriate for analysing time-averaged hemodynamics in resistance-dominated vascular regions, such as arterioles and capillaries, where viscous effects dominate the inertial effects.

Considering that blood behaves as a viscous, incompressible fluid, its motion through the cardiovascular system must respect the laws of conservation of momentum and conservation of mass, which are described by a group of governing equations. The Finite Element Method (FEM) is used to numerically solve these governing equations in complex geometries. A steady-state formulation is adopted because the imposed inlet flow is constant and the primary interest lies in the spatial distribution of velocity and pressure rather than transient dynamics.

## 2. Strong Form Equation

Blood flow within the vascular domain $`\Omega`$ is governed by the incompressible Navier-Stokes equations:

**Momentum Conservation**

```math
\rho \mathbf{v} \cdot \nabla \mathbf{v} - \mu \Delta \mathbf{v} + \nabla p = 0 \quad \text{on } \Omega
```

**Mass Conservation**

```math
\nabla \cdot \mathbf{v} = 0 \quad \text{on } \Omega
```

where $`\mathbf{v}`$ is the velocity vector field, $`p`$ the pressure, $`\rho`$ the fluid density, and $`\mu`$ the viscosity of the fluid.

The non-linear convective term $`\rho \mathbf{v} \cdot \nabla \mathbf{v}`$ represents the momentum transport achieved by advection, the Laplacian term $`\mu \Delta \mathbf{v}`$ models the viscous dissipation, and $`\nabla p`$ the pressure gradient provides the driving force for flow.

## 3. Boundary Conditions

Using physiological constraints, the boundary conditions are imposed on the domain boundaries, which are divided into three distinct regions as illustrated in Figure 1:

<div align="center">
  <img width="242" height="210" alt="image" src="https://github.com/user-attachments/assets/0b107943-caf0-4b58-ba09-b338efa67f7c" />
  <br/>
  <b>Figure 1:</b> Branching Vessel Tree with color-coded boundary regions
</div>


**Inlet boundary ($`\Gamma_1`$)**, where a velocity Dirichlet condition is imposed:
```math
\mathbf{v} = \mathbf{v}_0 = (-50, 0)^T \text{ mm s}^{-1} \quad \text{on } \Gamma_1
```

**Outlet boundaries ($`\Gamma_n = \Gamma_2 \cup \Gamma_3 \cup \Gamma_4 \cup \Gamma_5`$)**, with zero normal traction, Neumann conditions are imposed:
```math
(\mu \nabla \mathbf{v} - pI) \cdot \mathbf{n} = 0 \quad \text{on } \Gamma_n
```

**Vessel walls ($`\Gamma_6`$)**, with a no-slip condition applied, represent the stationary walls:
```math
\mathbf{v} = 0 \quad \text{on } \Gamma_6
```

The traction-free condition allows natural outflow without the need for prescribing specific pressure or velocity values.

## 4. Weak Form Equation

### 4.1 Momentum Equation Weak Form

To achieve the weak form, the momentum equation is multiplied by a vector test function $`\mathbf{w} \in [H^1(\Omega)]^2`$ and integrated over the domain:

```math
\int_\Omega \left( \rho (\mathbf{v} \cdot \nabla) \mathbf{v} - \mu \Delta \mathbf{v} + \nabla p \right) \cdot \mathbf{w} \, d\Omega = 0
```

Applying integration by parts to the viscous and pressure gradient terms gives:

```math
\int_\Omega \rho (\mathbf{v} \cdot \nabla \mathbf{v}) \cdot \mathbf{w} \, d\Omega - \int_\Gamma [(\mu \nabla \mathbf{v} - pI) \cdot \mathbf{n}] \cdot \mathbf{w} \, d\Gamma + \int_\Omega \mu \nabla \mathbf{v} : \nabla \mathbf{w} \, d\Omega - \int_\Omega p \nabla \cdot \mathbf{w} \, d\Omega = 0
```

Looking at the boundary flux terms, which can also be written as:

```math
\int_{\Gamma_1} (\mu \nabla \mathbf{v} - pI) \cdot \mathbf{n} \cdot \mathbf{w} \, d\Gamma_1 + \int_{\Gamma_6} (\mu \nabla \mathbf{v} - pI) \cdot \mathbf{n} \cdot \mathbf{w} \, d\Gamma_6 + \int_{\Gamma_n} (\mu \nabla \mathbf{v} - pI) \cdot \mathbf{n} \cdot \mathbf{w} \, d\Gamma_n
```

This boundary integral term vanishes because the test function $`\mathbf{w} = 0`$ on $`\Gamma_1`$ and $`\Gamma_6`$ (Dirichlet boundary conditions) and the traction $`(\mu \nabla \mathbf{v} - pI) \cdot \mathbf{n} = 0`$ on $`\Gamma_2 \cup \Gamma_3 \cup \Gamma_4 \cup \Gamma_5`$ (Neumann boundary conditions).

### 4.2 Continuity Equation Weak Form

The weak form of the conservation of mass is achieved by multiplying the continuity equation by a scalar test function $`q \in P = L^2(\Omega)`$ and integrating over the domain:

```math
\int_\Omega (\nabla \cdot \mathbf{v}) \cdot q \, d\Omega = 0
```

### 4.3 Continuous Weak Form

Find $`(\mathbf{v}, p) \in V \times P`$ such that,

```math
a(\mathbf{v}, \mathbf{w}) + b(p, \mathbf{w}) + b(q, \mathbf{v}) = 0 \quad \forall (\mathbf{w}, q) \in V_0 \times P
```

With function spaces:

```math
V_0 = \{\mathbf{w} \in [H^1(\Omega)]^2 : \mathbf{w} = 0 \text{ on } \Gamma_1 \cup \Gamma_6\}
```

```math
V = \{\mathbf{v} \in [H^1(\Omega)]^2 : \mathbf{v} = \mathbf{v}_0 \text{ on } \Gamma_1, \mathbf{v} = 0 \text{ on } \Gamma_6\}
```

### 4.4 Bilinear Operators

**Momentum term:**

```math
a(\mathbf{v}, \mathbf{w}) = \int_\Omega \rho (\mathbf{v} \cdot \nabla \mathbf{v}) \cdot \mathbf{w} \, d\Omega + \int_\Omega \mu \nabla \mathbf{v} : \nabla \mathbf{w} \, d\Omega
```

**Pressure-velocity term:**

```math
b(p, \mathbf{w}) = -\int_\Omega p (\nabla \cdot \mathbf{w}) \, d\Omega
```

**Continuity constraint:**

```math
b(q, \mathbf{v}) = \int_\Omega (\nabla \cdot \mathbf{v}) \cdot q \, d\Omega
```

## 5. Spatial Discretisation

In terms of mesh generation, the domain is discretized using an unstructured triangular mesh to accommodate the complex branching geometry of the arterial tree.

For the element spaces, a Taylor-Hood (P2-P1) mixed finite element formulation is employed, with velocity ($`\mathbf{v}`$) being continuous piecewise quadratic (P2) polynomials and pressure ($`p`$) being continuous piecewise linear (P1) polynomials. This satisfies the inf-sup stability condition.

Writing the function in the discrete space gives:

Find $`(\mathbf{v}_h, p_h) \in V^h \times P^h`$ such that,

```math
a(\mathbf{v}_h, \mathbf{w}_h) + b(p_h, \mathbf{w}_h) + b(q_h, \mathbf{v}_h) = 0 \quad \forall (\mathbf{w}_h, q_h) \in V_0^h \times P^h
```

**Discrete velocity space:**

```math
V^h := \left\{ \mathbf{v}_h \in [H^1(\Omega)]^2 \mid \mathbf{v} = \sum_{i=1}^{N_v} \sum_{k=1}^2 v_k^i \hat{\phi}_i, \text{ for some } \{v_k^i\} \in \mathbb{R} \right\}
```

and $`V_0^h = \{\mathbf{w}_h \in V^h : \mathbf{w}_h = 0 \text{ on Dirichlet boundaries}\}`$

**Discrete pressure space:**

```math
P^h := \left\{ q_h \in L^2(\Omega) \mid q = \sum_{i=1}^{N_p} q^i \hat{\psi}_i, \text{ for some } \{q^i\} \in \mathbb{R} \right\}
```

After assembly and linearisation, the discrete system at each iteration takes the block matrix form:

```math
\begin{bmatrix} A & B \\ B^T & 0 \end{bmatrix} \begin{bmatrix} U_v \\ U_p \end{bmatrix} = \begin{bmatrix} R_v \\ R_p \end{bmatrix}
```

where $`U_v`$ and $`U_p`$ are the global velocity and pressure degrees of freedom, $`R_v`$ is the momentum residual vector, and $`R_p`$ is the continuity residual vector.

## 6. Numerical Solution

The non-linear system resulting from the discretized equations is solved using a Newton-Raphson iterative scheme. This involves assembling the momentum and continuity residuals, after which the boundary conditions are enforced on both residuals and the corresponding Jacobian blocks. The resulting linearised system is then solved at each iteration to update the velocity and pressure fields. Iterations continue until the residual norm falls below a tolerance of $`10^{-6}`$.

## 7. Numerical Implementation and Workflow

The FEM solver is implemented in MATLAB following a structured workflow illustrated in Figure 2. The first step of this solver begins with loading the mesh and defining the vascular domain $`\Omega`$, following with the generation of P2-P1 basis functions for Taylor Hood elements. Subsequently, function spaces for velocity and pressure are initialized with physical parameters, and boundary conditions are assigned. The Newton-Raphson iteration loop assembles the momentum and residual vectors, applies the Dirichlet boundary conditions at inlets and walls, and finally checks the convergence. If not yet converged, the Jacobian Blocks are assembled, boundary conditions are enforced on the global system, and the linear system is solved to update the velocity and pressure fields. Post-processing steps include computing outflow rates through terminal boundaries ($`\Gamma_2`$-$`\Gamma_5`$) by averaging velocity at edge endpoints and integrating normal components using the midpoint rule, verifying mass conservation between inlet and total outlet flow, and calculating total pressure drop using length-weighted averaging of nodal pressures along boundary edges.

<p align="center">
  <img src="https://github.com/user-attachments/assets/04bdd25a-438a-41bd-9d67-d425da634f64" alt="image" width="411"/>
  <br/>
  <b>Figure 2:</b> Diagram of computational workflow for FEM solution of arterial flow.
</p>


----

# Results

## 1. Convergence and Numerical Accuracy

The nonlinear system that arises from the convective term, was solved using the Newton-Raphson algorithm. Figure 3 presents the convergence history showing systematic residual norm reduction across four progressively refined meshes.

The trend of this algorithm resembles that of a quadratic convergence, with residual norms decreasing roughly by an order of magnitude at each iteration in the final stages. This type of rapid convergence is characteristic of Newton-Raphson methods near the solution. The sparser Mesh 0 required an additional iteration to converge.

<p align="center">
  <img src="https://github.com/user-attachments/assets/89c4fcef-925b-4429-b420-d1c8bdc46bcd" alt="Convergence history" width="600"/>
  <br/>
  <b>Figure 3:</b> Convergence history for Newton-Raphson. All meshes achieved convergence tolerance of 10<sup>-6</sup> in 5 to 6 iterations.
</p>

## 2. Mesh Convergence

The four progressively refined meshes used were analysed to assess solution convergence. Table 1 summarizes the main mesh characteristics and highlights how the spatial resolution improves through refinement. As can be seen in this table, the maximum element size systematically reduces by approximately a factor of 2 between the different meshes. Figure 4 visualises the refinement of Mesh 0 and Mesh 2, clearly showing the improved spatial resolution everywhere, but most importantly in the more complex regions, such as near branching points and vessel walls, where the velocity gradients are the steepest.

**Table 1:** Mesh characteristics across refinement levels

| Mesh | Elements | Nodes | DOF | Min element size (mm) | Max element size (mm) |
|------|----------|-------|-----|----------------------|----------------------|
| 0 | 448 | 1031 | 2354 | 0.785 | 8.964 |
| 1 | 705 | 1598 | 3643 | 0.745 | 4.662 |
| 2 | 1906 | 4087 | 9265 | 0.746 | 2.411 |
| 3 | 7491 | 15514 | 35040 | 0.676 | 1.177 |

![image](https://github.com/user-attachments/assets/19f87d36-52e5-46cf-952c-cfbbdd82572a)

**Figure 4:** Visual comparison of Mesh 0 and Mesh 2 showing increased resolution.

## 3. Mass Conservation Error Analysis

A critical validation metric is global mass conservation. Mass should be conserved between inflow 1 and outflows 2 to 5; we should not lose any blood flow. Therefore, the continuity equation requires:

```math
\int_{\Gamma_1} \mathbf{v} \cdot \mathbf{n} \, d\Gamma = \sum_{i=2}^5 \int_{\Gamma_i} \mathbf{v} \cdot \mathbf{n} \, d\Gamma
```

Table 2 represents the flow rates for each mesh and the conservation error, illustrating an overall convergence towards an acceptable conservation error.

As can be seen in the table, Mesh 0 exhibits a severe mass conservation violation (55.5%), indicating inadequate resolution. This error systematically reduces with refinement, reaching 0.4% with mesh 3, confirming numerical convergence. An acceptable error is considered below 2%; therefore, only meshes 2 and 3 achieve acceptable errors.

**Table 2:** Conservation of mass analysis across meshes

| Mesh | Inlet Flow (mm³/s) | Total Outlet Flow $`\Gamma_n = \Gamma_2 \cup \Gamma_3 \cup \Gamma_4 \cup \Gamma_5`$ (mm³/s) | Conservation Error (%) |
|------|-------------------|--------------------------------------------------------------------------------------|----------------------|
| 0 | 985.2 | 438.1 | 55.5 |
| 1 | 1056.1 | 922.5 | 12.7 |
| 2 | 1054.7 | 1037.1 | 1.7 |
| 3 | 1107.4 | 1102.6 | 0.4 |

## 4. Velocity and Pressure Field Distribution

<table>
<tr>
<td align="center">
  <img src="https://github.com/user-attachments/assets/bdcd0768-757f-407c-89e6-0ee9d8954914" alt="Velocity Magnitude Field" width="400"/>
  <br/>
  <b>Figure 5:</b> Velocity Magnitude Field (Mesh 3)
</td>
<td align="center">
  <img src="https://github.com/user-attachments/assets/7e8f2d2e-accf-47f0-a0a8-e6792b1d41c4" alt="Pressure Field Distribution" width="400"/>
  <br/>
  <b>Figure 6:</b> Pressure Field Distribution (Mesh 3)
</td>
</tr>
</table>


Figure 5 shows the converged velocity magnitude field for the most refined mesh, Mesh 3. This figure shows the spatial distribution of flow through the branching arterial network. It can be noted that the defined velocities of 50 mm/s correctly occur at the inlet boundary. Additionally, asymmetric flow patterns at bifurcations, with some branches having more flow.

Figure 6 shows the corresponding pressure field distribution. The pressure decreases monotonically from the inlet to the outlets, with the steepest gradients occurring in regions of highest flow resistance, which are the branching points.

## 5, Flow Distribution Through Terminal Branches

Table 3 presents the converged flow distribution for Mesh 3, directly addressing the primary research question about how the inlet flow is distributed along the terminals.

It can be noticed that the distribution of flow is highly asymmetric with $`\Gamma_2`$ receiving 45% of the total flow, which is 2.5 times more than other outlets, despite all boundaries having a comparable amount of boundary edges. This pattern is present across all meshes, signifying it's a feature of our model rather than a numerical artefact.

**Table 3:** Converged outlet flow distribution for mesh 3 with inlet flow = 1107.4 mm³/s

| Outlet | Flow Rate (mm³/s) | Percentage of Inlet | Number of edges |
|--------|------------------|-------------------|----------------|
| $`\Gamma_2`$ | 498.0 | 45.0% | 8 |
| $`\Gamma_3`$ | 195.7 | 17.7% | 8 |
| $`\Gamma_4`$ | 193.7 | 17.5% | 8 |
| $`\Gamma_5`$ | 215.1 | 19.4% | 8 |
| **Total Outflow** | **1102.6** | **99.6%** | **32** |

## 6. Total Pressure Drop and Pressure Field

The total pressure drop across the arterial network can be calculated using:

```math
\Delta p_{\text{total}} = p_{\text{inlet}} - \frac{1}{4} \sum_{i=2}^5 p_i
```

Table 4 represents the output of this calculation. The value of the pressure drop is 2 to 3 orders of magnitude lower than typical physiological arteriolar pressure drops (40-60 mmHg). This assumption can be verified using Poiseuille's Law:

```math
\Delta p = \frac{8 \mu L Q}{\pi r^4} \approx 4.1 \text{ Pa}
```

These two values agree closely, confirming that FEM solution correctly captures local pressure gradients for the modeled geometry.

**Table 4:** Pressure analysis for mesh 3

| Parameter | Value in Pa | Value in mmHg |
|-----------|------------|---------------|
| Mean Inlet pressure | 4.4182 | 0.0331 |
| Mean Outlet pressure | 0.0001 | 0.0001 |
| Pressure Drop $`\Delta p`$ | 4.4181 | 0.0331 |

## 7. Error Analysis

Following the finite element error theory, total error decomposes as:

```math
\varepsilon = \|u_m - u_h\| \leq \|u_m - u\| + \|u - u_h\|
```

where $`u_m`$ is the true physical solution, $`u`$ is the continuous weak form solution, and $`u_h`$ is the discrete FEM approximation. The first term represents modelling error, which arises from mathematical simplifications, and physical assumptions. The second term represents the simulation errors appearing from discretization, numerical integration and iterative convergence.

In terms of simulation errors, the mass conservation error of 0.43% for mesh 3 satisfies the incompressibility constraint within an acceptable tolerance below 1%. And for modelling errors, the Newton-Raphson convergence tolerance $`10^{-6}`$ ensures iterative solution error is negligible compared to other error sources. The main source of error is the traction-free boundary condition, which is inconsistent with the physics.
