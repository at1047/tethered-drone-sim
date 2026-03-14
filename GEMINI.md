# Capstone Simulation: Quadcopter with Slung Payload

This project is a MATLAB-based simulation environment for a quadcopter carrying a cable-suspended payload. It focuses on deriving the system dynamics using Lagrangian mechanics and implementing various control strategies to manage the payload's swing.

## Rules
- Don't delete commented out code
- Don't change variable names unless explicitly requested
- Don't change how the code flows unless explicitly requested
- No need to execute files after making changes, I will execute and provide error messages if exists

## Project Overview

- **Purpose**: To simulate and control the dynamics of a quadcopter-payload system, specifically addressing vibration/swing reduction.
- **Core Technology**: MATLAB, Symbolic Math Toolbox.
- **Architecture**:
    - **Dynamics Derivation**: Symbolic scripts derive the Equations of Motion (EoM) using the Lagrangian method. Both generalized coordinates (angles) and redundant coordinates (Lagrange multipliers for cable constraints) are supported.
    - **Simulation**: Numerical solvers (`ode45` or discrete-time loops) integrate the dynamics over time.
    - **Control**: Includes PD-based position control and Input Shaping techniques to minimize residual oscillations.
    - **Optimization**: Symbolic matrices (Mass, Coriolis, Gravity, Input) are exported to optimized MATLAB functions (`auto_*.m`) for high-performance numerical simulation.
