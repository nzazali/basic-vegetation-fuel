import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, solve, Eq

# Define symbols
V, F, lam = symbols('V F lambda')

# Define the functions
def f(V):
    """Function f(V) = V*(1 - V)"""
    return V * (1 - V)

def g_lambda(V, F, lam):
    """Function g_λ(V, F) = λ(V - F)"""
    return lam * (V - F)

# Set up the equations for equilibrium points
# f[V] == 0  and  g_λ[V, F] == 0
eq1 = Eq(f(V), 0)
eq2 = Eq(g_lambda(V, F, lam), 0)

print("Equations to solve:")
print(f"f(V) = 0: {eq1}")
print(f"g_λ(V, F) = 0: {eq2}")
print()

# Solve the system of equations
solutions = solve([eq1, eq2], [V, F])

print("Solutions (equilibrium points):")
for i, sol in enumerate(solutions):
    print(f"Solution {i+1}: V = {sol[0]}, F = {sol[1]}")

# Alternative: solve each equation separately to see all solutions
print("\nSeparate analysis:")
print("Solutions to f(V) = 0:")
v_solutions = solve(eq1, V)
print(f"V = {v_solutions}")

print("\nSolutions to g_λ(V, F) = 0:")
f_solutions = solve(eq2, F)
print(f"F = {f_solutions}")

# For specific value of lambda, you can substitute
print("\nExample with λ = 1:")
lam_val = 1
eq2_specific = Eq(g_lambda(V, F, lam_val), 0)
solutions_specific = solve([eq1, eq2_specific], [V, F])
print(f"Solutions: {solutions_specific}")

# =============================================================================
# PHASE PORTRAIT PLOT (converted from Mathematica Show[] command)
# =============================================================================

# Set lambda value
lambda_val = 0.636

# Define numerical versions of the functions
def f_num(v):
    return v * (1 - v)

def g_num(v, f_val, lam_val):
    return lam_val * (v - f_val)

# Create meshgrid for plotting
v_range = np.linspace(-0.02, 1.5, 30)
f_range = np.linspace(-0.02, 1.5, 30)
V_mesh, F_mesh = np.meshgrid(v_range, f_range)

# Calculate vector field components
dV_dt = f_num(V_mesh)
dF_dt = g_num(V_mesh, F_mesh, lambda_val)

# Create the plot
fig, ax = plt.subplots(figsize=(10, 8))

# 1. Vector field (quiver plot) - equivalent to VectorPlot
ax.quiver(V_mesh, F_mesh, dV_dt, dF_dt,
          alpha=0.6, color='lightgray', scale=50, width=0.002)

# 2. Stream plot - equivalent to StreamPlot
ax.streamplot(V_mesh, F_mesh, dV_dt, dF_dt,
              color='black', density=1.5, linewidth=0.8)

# 3. Nullclines - equivalent to ContourPlot
v_fine = np.linspace(-0.02, 1.5, 200)
f_fine = np.linspace(-0.02, 1.5, 200)
V_fine, F_fine = np.meshgrid(v_fine, f_fine)

# f(V) = 0 nullcline (blue dashed)
dV_dt_fine = f_num(V_fine)
ax.contour(V_fine, F_fine, dV_dt_fine, levels=[0],
           colors=['blue'], linestyles=['dashed'], linewidths=2)

# g_λ(V,F) = 0 nullcline (red dashed)
dF_dt_fine = g_num(V_fine, F_fine, lambda_val)
ax.contour(V_fine, F_fine, dF_dt_fine, levels=[0],
           colors=['red'], linestyles=['dashed'], linewidths=2)

# 4. Equilibrium points - equivalent to ListPlot
equilibrium_points = [(0, 0), (1, 1)]
for point in equilibrium_points:
    ax.plot(point[0], point[1], 'ko', markersize=8, markerfacecolor='black')

# 5. Special point with white center - equivalent to Epilog
ax.plot(0, 0, 'ko', markersize=6, markerfacecolor='white',
        markeredgecolor='black', markeredgewidth=2)

# Formatting - equivalent to Frame, Axes, TicksStyle, FrameLabel
ax.set_xlim(-0.02, 1.5)
ax.set_ylim(-0.02, 1.5)
ax.set_xlabel('Vegetation', fontsize=16, color='black')
ax.set_ylabel('Fuels', fontsize=16, color='black')
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=12)

# Add title
ax.set_title(f'Phase Portrait (λ = {lambda_val})', fontsize=14)

plt.tight_layout()
plt.show()

# Print information about the plot
print(f"\nPhase Portrait Analysis (λ = {lambda_val}):")
print("- Gray arrows: Vector field showing flow direction")
print("- Black lines: Streamlines (trajectories)")
print("- Blue dashed: f(V) = 0 nullcline (V = 0 and V = 1)")
print("- Red dashed: g_λ(V,F) = 0 nullcline (V = F)")
print("- Black dots: Equilibrium points at (0,0) and (1,1)")
print("- White-centered dot: Special emphasis on origin")
