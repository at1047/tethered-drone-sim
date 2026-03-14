import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation

# --- Parameters ---
m_q = 0.7
m_p = 0.7
l = 0.5
g = 9.81

def drone_velocity_dynamics(t, s):
    # Extract states (Convention: [xq, yq, zq, dxq, dyq, dzq, xp, yp, zp, dxp, dyp, dzp])
    xq, yq, zq = s[0:3]
    dxq, dyq, dzq = s[3:6]
    xp, yp, zp = s[6:9]
    dxp, dyp, dzp = s[9:12]

    # --- Command Input (Circular Trajectory) ---
    R = 1.0
    omega = 1.5
    vx_cmd = -R * omega * np.sin(omega * t)
    vy_cmd =  R * omega * np.cos(omega * t)
    vz_cmd = 0
    
    ddxq_cmd = -R * omega**2 * np.cos(omega * t)
    ddyq_cmd = -R * omega**2 * np.sin(omega * t)
    ddzq_cmd = 0

    # --- Controller ---
    kv = 20
    f1 = m_q * (ddxq_cmd + kv * (vx_cmd - dxq))
    f2 = m_q * (ddyq_cmd + kv * (vy_cmd - dyq))
    f3 = m_q * (ddzq_cmd + kv * (vz_cmd - dzq)) + (m_q + m_p) * g

    # ==========================================================
    # PASTE YOUR 6 EoM HERE (from lagrangian_drone_dynamics.py)
    # ==========================================================
    # Note: Ensure you use numpy functions (e.g., np.sin instead of sin) if needed.
    # The current expressions from sympy are already in terms of symbols.
    
    ddxq = 1.0*(f1*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) + m_p*(xp - xq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddyq = 1.0*(f2*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) + m_p*(yp - yq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddzq = 1.0*(m_p*(zp - zq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq) + (f3 - g*m_q)*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))/(m_q*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2))
    ddxp = -1.0*(xp - xq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq)/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    ddyp = -1.0*(yp - yq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq)/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    ddzp = 1.0*(-g*(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2) - (zp - zq)*(dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq))/(m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
    # ==========================================================

    # Assemble derivative
    return [dxq, dyq, dzq, ddxq, ddyq, ddzq, dxp, dyp, dzp, ddxp, ddyp, ddzp]

# --- Initial Conditions ---
xq0, yq0, zq0 = 0, 0, 2
xp0, yp0, zp0 = xq0, yq0, zq0 - l
s0 = [xq0, yq0, zq0, 0, 0, 0, xp0, yp0, zp0, 0, 0, 0]

t_span = (0, 10)
t_eval = np.linspace(0, 10, 500)

# --- Run Simulation ---
sol = solve_ivp(drone_velocity_dynamics, t_span, s0, t_eval=t_eval, rtol=1e-8, atol=1e-8)

# --- Plotting ---
t = sol.t
s = sol.y.T

fig = plt.figure(figsize=(12, 8))

# Drone Position
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(t, s[:, 0:3])
ax1.legend(['xq', 'yq', 'zq'])
ax1.set_title('Drone Position')
ax1.grid(True)

# Payload Position
ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(t, s[:, 6:9])
ax2.legend(['xp', 'yp', 'zp'])
ax2.set_title('Payload Position')
ax2.grid(True)

# Pendulum Angles
ax3 = fig.add_subplot(2, 2, 3)
th = np.arccos(np.clip((s[:, 2] - s[:, 8]) / l, -1, 1))
phi = np.atan2(s[:, 7] - s[:, 1], s[:, 6] - s[:, 0])
ax3.plot(t, np.degrees(th), label='theta')
ax3.plot(t, np.degrees(phi), label='phi')
ax3.legend()
ax3.set_title('Pendulum Angles (deg)')
ax3.grid(True)

# 3D Trajectory
ax4 = fig.add_subplot(2, 2, 4, projection='3d')
ax4.plot(s[:, 0], s[:, 1], s[:, 2], 'b', label='Drone')
ax4.plot(s[:, 6], s[:, 7], s[:, 8], 'r--', label='Payload')
ax4.legend()
ax4.set_title('3D Trajectory')
ax4.set_xlabel('X')
ax4.set_ylabel('Y')
ax4.set_zlabel('Z')
ax4.set_aspect('equal')

plt.tight_layout()

# --- Animation ---
fig_anim = plt.figure()
ax_anim = fig_anim.add_subplot(111, projection='3d')

def update(i):
    ax_anim.clear()
    ax_anim.plot([s[i, 0], s[i, 6]], [s[i, 1], s[i, 7]], [s[i, 2], s[i, 8]], 'o-', lw=2)
    ax_anim.set_xlim([-2, 2])
    ax_anim.set_ylim([-2, 2])
    ax_anim.set_zlim([0, 3])
    ax_anim.set_xlabel('X')
    ax_anim.set_ylabel('Y')
    ax_anim.set_zlabel('Z')
    ax_anim.set_title(f'Time: {t[i]:.2f}s')

# Create animation (only showing first 100 frames to keep it fast)
# ani = animation.FuncAnimation(fig_anim, update, frames=range(0, len(t), 5), interval=50)

plt.show()
