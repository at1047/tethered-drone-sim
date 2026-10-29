import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from mpc_controller import MPCController

# --- Parameters ---
m_q = 0.7
m_p = 0.7
l = 0.5
g = 9.81

# Create the controller instance (MPC)
ctrl = MPCController(m_q, m_p, l, g)

def drone_velocity_dynamics(t, s):
    # Extract states (Convention: [xq, yq, zq, dxq, dyq, dzq, xp, yp, zp, dxp, dyp, dzp])
    xq, yq, zq = s[0:3]
    dxq, dyq, dzq = s[3:6]
    xp, yp, zp = s[6:9]
    dxp, dyp, dzp = s[9:12]

    # --- Desired Setpoint ---
    # target_xq = 1.0, target_yq = 0.0, target_zq = 2.0
    s_ref = [1.0, 1.0, 2.0, 0, 0, 0, 1.0, 1.0, 1.5, 0, 0, 0]
    
    # --- Controller (MPC) ---
    f1, f2, f3 = ctrl.compute_force(t, s, s_ref)

    # ==========================================================

    # PASTE YOUR 6 EoM HERE (from lagrangian_drone_dynamics.py)
    # ==========================================================
    
    ddxq = 1.0*(f1*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) + m_p*(xp - xq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddyq = 1.0*(f2*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) + m_p*(yp - yq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddzq = 1.0*(m_p*(zp - zq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq) + (f3 - g*m_q)*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddxp = -1.0*(xp - xq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq)/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    ddyp = -1.0*(yp - yq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq)/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    ddzp = 1.0*(-g*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) - (zp - zq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    # ==========================================================

    # --- Air Resistance (Drag) ---
    Cd_q = 0.1286 # Drone drag coefficient
    Cd_p = 0.05   # Payload drag coefficient (approximate)

    # Apply drag accelerations to drone
    ddxq -= (Cd_q * dxq) / m_q
    ddyq -= (Cd_q * dyq) / m_q
    ddzq -= (Cd_q * dzq) / m_q

    # Apply drag accelerations to payload
    ddxp -= (Cd_p * dxp) / m_p
    ddyp -= (Cd_p * dyp) / m_p
    ddzp -= (Cd_p * dzp) / m_p

    # Assemble derivative
    return [dxq, dyq, dzq, ddxq, ddyq, ddzq, dxp, dyp, dzp, ddxp, ddyp, ddzp]

# --- Initial Conditions ---
xq0, yq0, zq0 = 0, 0, 2
xp0, yp0, zp0 = xq0, yq0, zq0 - l
s0 = [xq0, yq0, zq0, 0, 0, 0, xp0, yp0, zp0, 0, 0, 0]

t_span = (0, 10)
t_eval = np.linspace(0, 10, 500)

def drone_velocity_dynamics_with_forces(t, s, f1, f2, f3):
    # Extract states (Convention: [xq, yq, zq, dxq, dyq, dzq, xp, yp, zp, dxp, dyp, dzp])
    xq, yq, zq = s[0:3]
    dxq, dyq, dzq = s[3:6]
    xp, yp, zp = s[6:9]
    dxp, dyp, dzp = s[9:12]

    # --- Physics (Using auto-generated redundant EOM) ---
    # (Same as before but uses provided forces)
    
    # ==========================================================
    # PASTE YOUR 6 EoM HERE (from lagrangian_drone_dynamics.py)
    # ==========================================================
    ddxq = 1.0*(f1*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) + m_p*(xp - xq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddyq = 1.0*(f2*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) + m_p*(yp - yq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddzq = 1.0*(m_p*(zp - zq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq) + (f3 - g*m_q)*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddxp = -1.0*(xp - xq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq)/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    ddyp = -1.0*(yp - yq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq)/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    ddzp = 1.0*(-g*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) - (zp - zq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    # ==========================================================

    # --- Air Resistance (Drag) ---
    k_q = 0.1286 # Drone drag coefficient
    k_p = 0.05   # Payload drag coefficient (approximate)

    # Apply drag accelerations to drone
    ddxq -= (k_q * dxq) / m_q
    ddyq -= (k_q * dyq) / m_q
    ddzq -= (k_q * dzq) / m_q

    # Apply drag accelerations to payload
    ddxp -= (k_p * dxp) / m_p
    ddyp -= (k_p * dyp) / m_p
    ddzp -= (k_p * dzp) / m_p

    # Assemble derivative
    return [dxq, dyq, dzq, ddxq, ddyq, ddzq, dxp, dyp, dzp, ddxp, ddyp, ddzp]

# --- Run Simulation (Discrete Time Loop for Speed) ---
dt_sim = 0.05
t_end = 20.0
t_steps = np.arange(0, t_end, dt_sim)
s_hist = [s0]
t_hist = [0]

print("Starting simulation... (MPC is computationally heavy)")
s_current = np.array(s0)

for i in range(len(t_steps)-1):
    t_curr = t_steps[i]
    
    # --- Desired Setpoint ---
    s_ref = [1.0, 0.0, 2.0, 0, 0, 0, 1.0, 0.0, 1.5, 0, 0, 0]
    
    # Update MPC at fixed intervals
    f1, f2, f3 = ctrl.compute_force(t_curr, s_current, s_ref)
    
    # Small integration step for physics
    # We use a tighter tolerance here but it's only for a small dt
    sol_step = solve_ivp(drone_velocity_dynamics_with_forces, [t_curr, t_curr + dt_sim], 
                         s_current, args=(f1, f2, f3), rtol=1e-4, atol=1e-4)
    
    s_current = sol_step.y[:, -1]
    s_hist.append(s_current)
    t_hist.append(t_curr + dt_sim)
    
    if i % 20 == 0:
        print(f"Simulation Progress: {100*t_curr/t_end:.1f}%")

t = np.array(t_hist)
s = np.array(s_hist)

# --- Plotting ---
fig = plt.figure(figsize=(15, 5))

# Drone Position
ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(t, s[:, 0:3])
ax1.legend(['xq', 'yq', 'zq'])
ax1.set_title('Drone Position')
ax1.grid(True)

# Payload Position
ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(t, s[:, 6:9])
ax2.legend(['xp', 'yp', 'zp'])
ax2.set_title('Payload Position')
ax2.grid(True)

# 3D Trajectory
ax3 = fig.add_subplot(1, 3, 3, projection='3d')
ax3.plot(s[:, 0], s[:, 1], s[:, 2], 'b', label='Drone')
ax3.plot(s[:, 6], s[:, 7], s[:, 8], 'r--', label='Payload')
ax3.legend()
ax3.set_title('3D Trajectory')
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.set_zlabel('Z')
ax3.set_aspect('equal')

plt.tight_layout()

# --- Animation ---
fig_anim = plt.figure('Drone Animation')
ax_anim = fig_anim.add_subplot(111, projection='3d')

# Initialize the line object
line, = ax_anim.plot([], [], [], 'o-', lw=2)

def init():
    ax_anim.set_xlim([-2, 2])
    ax_anim.set_ylim([-2, 2])
    ax_anim.set_zlim([0, 3])
    ax_anim.set_xlabel('X (m)')
    ax_anim.set_ylabel('Y (m)')
    ax_anim.set_zlabel('Z (m)')
    return line,

def update(i):
    # Line between Drone (0:3) and Payload (6:9)
    x = [s[i, 0], s[i, 6]]
    y = [s[i, 1], s[i, 7]]
    z = [s[i, 2], s[i, 8]]
    line.set_data(x, y)
    line.set_3d_properties(z)
    ax_anim.set_title(f'Time: {t[i]:.2f}s')
    return line,

# Create animation
ani = animation.FuncAnimation(fig_anim, update, frames=range(0, len(t), 2), 
                              init_func=init, blit=False, interval=20)

plt.show()
