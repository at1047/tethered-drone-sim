"""Lagrangian slung-payload simulation with a square trajectory.

Run with: python simulate_drone_velocity.py
Requires numpy, scipy, and matplotlib.
The drone is a force-actuated point mass and the payload hangs on a rigid,
massless tether whose fixed-length constraint is always active.
State: [xq, yq, zq, dxq, dyq, dzq, xp, yp, zp, dxp, dyp, dzp].
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation

# --- Parameters (same as the MATLAB baseline) ---
m_q = 0.7
m_p = 0.7
l = 0.5
g = 9.81
max_payload_force = 0.5  # N, magnitude of the horizontal force vector
kp_payload = 8.0  # 1/s^2, payload position feedback
kd_payload = 4.0  # 1/s, payload velocity feedback

# --- Square trajectory parameters ---
start_position = np.array([0.0, 0.0, 2.0])  # m
square_side = 2.0  # m, horizontal square; stop at each corner
target_position = start_position + np.array([square_side, 0.0, 0.0])
waypoints = np.array([start_position, target_position,
                      target_position + [0.0, square_side, 0.0],
                      start_position + [0.0, square_side, 0.0], start_position])
max_velocity = 0.5  # m/s, total reference and commanded speed
max_acceleration = 0.5  # m/s^2, reference acceleration


def segment_duration(start_position, target_position):
    if max_velocity <= 0 or max_acceleration <= 0:
        raise ValueError('max_velocity and max_acceleration must be positive')
    distance = np.linalg.norm(target_position - start_position)
    peak_velocity = min(max_velocity, np.sqrt(distance * max_acceleration))
    return distance / peak_velocity + peak_velocity / max_acceleration if distance else 0.0


def point_to_point_trajectory(t, start_position=start_position, target_position=target_position):
    """Return position, velocity, acceleration for a rest-to-rest straight move.

    Uses a trapezoidal speed profile, or a triangular one for short moves.
    Holds the destination after the move finishes.
    """
    if max_velocity <= 0 or max_acceleration <= 0:
        raise ValueError('max_velocity and max_acceleration must be positive')
    displacement = target_position - start_position
    distance = np.linalg.norm(displacement)
    if distance == 0:
        return start_position.copy(), np.zeros(3), np.zeros(3)
    direction = displacement / distance
    peak_velocity = min(max_velocity, np.sqrt(distance * max_acceleration))
    t_accel = peak_velocity / max_acceleration
    t_cruise = max(0.0, distance / peak_velocity - t_accel)
    t_finish = 2 * t_accel + t_cruise

    if t < 0:
        return start_position.copy(), np.zeros(3), np.zeros(3)
    if t < t_accel:
        position = 0.5 * max_acceleration * t**2
        velocity = max_acceleration * t
        acceleration = max_acceleration
    elif t < t_accel + t_cruise:
        position = 0.5 * peak_velocity * t_accel + peak_velocity * (t - t_accel)
        velocity = peak_velocity
        acceleration = 0.0
    elif t < t_finish:
        remaining = t_finish - t
        position = distance - 0.5 * max_acceleration * remaining**2
        velocity = max_acceleration * remaining
        acceleration = -max_acceleration
    else:
        return target_position.copy(), np.zeros(3), np.zeros(3)
    return (start_position + direction * position,
            direction * velocity, direction * acceleration)


def square_trajectory(t):
    """Traverse four rest-to-rest segments, then hold at the starting corner."""
    for start_position, target_position in zip(waypoints[:-1], waypoints[1:]):
        duration = segment_duration(start_position, target_position)
        if t < duration:
            return point_to_point_trajectory(t, start_position, target_position)
        t -= duration
    return waypoints[-1].copy(), np.zeros(3), np.zeros(3)


def payload_tracking_force(t, s):
    """Horizontal payload tracking force [Fx, Fy] in N, capped at 0.5 N."""
    position_ref, velocity_ref, acceleration_ref = square_trajectory(t)
    payload_force = m_p * (acceleration_ref[:2]
                           + kp_payload * (position_ref[:2] - s[6:8])
                           + kd_payload * (velocity_ref[:2] - s[9:11]))
    magnitude = np.linalg.norm(payload_force)
    if magnitude > max_payload_force:
        payload_force *= max_payload_force / magnitude
    return payload_force


def drone_velocity_dynamics(t, s, payload_force=None):
    """Coupled dynamics; optionally supply horizontal payload force [Fx, Fy] (N).

    Without an explicit input, use payload tracking feedback. Supplied inputs
    are also limited to max_payload_force; pass [0, 0] to disable assistance.
    """
    xq, yq, zq = s[0:3]
    dxq, dyq, dzq = s[3:6]
    xp, yp, zp = s[6:9]
    dxp, dyp, dzp = s[9:12]

    # --- Command Input (Square Trajectory) ---
    position_ref, velocity_ref, acceleration_ref = square_trajectory(t)
    # Position correction holds the endpoint despite velocity tracking errors.
    kp = 2.0
    velocity_cmd = velocity_ref + kp * (position_ref - s[0:3])
    speed_cmd = np.linalg.norm(velocity_cmd)
    if speed_cmd > max_velocity:
        velocity_cmd *= max_velocity / speed_cmd
    vx_cmd, vy_cmd, vz_cmd = velocity_cmd
    ddxq_cmd, ddyq_cmd, ddzq_cmd = acceleration_ref

    # Previous position setpoint, retained for reference:
    # target_xq = 1.0, target_yq = 0.0, target_zq = 2.0
    # --- Controller (velocity tracking with acceleration feedforward) ---
    kv = 20.0
    f1 = m_q * (ddxq_cmd + kv * (vx_cmd - dxq))
    f2 = m_q * (ddyq_cmd + kv * (vy_cmd - dyq))
    f3 = m_q * (ddzq_cmd + kv * (vz_cmd - dzq)) + (m_q + m_p) * g

    # --- Lagrangian dynamics with an always-active constraint ---
    # |r|^2 = l^2, where r points from the drone to the payload.
    # Eliminating the Lagrange multiplier using r·(a_p-a_q) + |v_rel|^2 = 0:
    # T = (|v_rel|^2 + r·(payload_force/m_p - force/m_q))
    #     / ((1/m_q + 1/m_p)*l).
    # The tether stays taut; there is no slack mode or length comparison.
    # PASTE YOUR 6 EoM HERE (from lagrangian_drone_dynamics.py)
    r = np.array([xp - xq, yp - yq, zp - zq])
    v_rel = np.array([dxp - dxq, dyp - dyq, dzp - dzq])
    force = np.array([f1, f2, f3])
    if payload_force is None:
        payload_force = payload_tracking_force(t, s)
    payload_force = np.array([payload_force[0], payload_force[1], 0.0])
    magnitude = np.linalg.norm(payload_force)
    if magnitude > max_payload_force:
        payload_force *= max_payload_force / magnitude
    gravity = np.array([0.0, 0.0, -g])
    tension = (np.dot(v_rel, v_rel) + np.dot(r, payload_force / m_p - force / m_q)) / ((1 / m_q + 1 / m_p) * l)
    ddxq, ddyq, ddzq = force / m_q + gravity + tension * r / (m_q * l)
    ddxp, ddyp, ddzp = payload_force / m_p + gravity - tension * r / (m_p * l)

    return [dxq, dyq, dzq, ddxq, ddyq, ddzq,
            dxp, dyp, dzp, ddxp, ddyp, ddzp]


def main():
    # --- Initial Conditions ---
    xq0, yq0, zq0 = start_position
    xp0, yp0, zp0 = xq0, yq0, zq0 - l
    s0 = [xq0, yq0, zq0, 0, 0, 0, xp0, yp0, zp0, 0, 0, 0]
    t_finish = sum(segment_duration(a, b) for a, b in zip(waypoints[:-1], waypoints[1:]))
    t_span = (0.0, t_finish + 5.0)
    t_eval = np.linspace(*t_span, int(np.ceil(t_span[1] / 0.02)) + 1)

    # Both cases use identical initial conditions, trajectory, solver, and sample times.
    position_ref = np.array([square_trajectory(ti)[0] for ti in t_eval])
    # Desired payload position is directly below the commanded drone position.
    # Error is simulated payload position minus this reference, not swing offset.
    payload_ref = position_ref + np.array([0.0, 0.0, -l])
    results = {}
    for label, assisted in [('Force disabled', False), ('Force enabled', True)]:
        sol = solve_ivp(drone_velocity_dynamics, t_span, s0, t_eval=t_eval,
                        args=(None if assisted else [0.0, 0.0],),
                        rtol=1e-8, atol=1e-10, max_step=0.02)
        if not sol.success:
            raise RuntimeError(f'{label}: {sol.message}')
        t = sol.t
        s = sol.y.T
        payload_error = s[:, 6:9] - payload_ref
        payload_error_norm = np.linalg.norm(payload_error, axis=1)
        # Summary metrics use the full run, except hold RMSE (after square completion).
        metrics = np.array([np.sqrt(np.mean(payload_error_norm**2)),
                            payload_error_norm.max(),
                            np.sqrt(np.mean(payload_error_norm[t >= t_finish]**2))])
        results[label] = dict(states=s, error=payload_error, error_norm=payload_error_norm,
                              metrics=metrics)
        print(f'{label}: 3D RMSE={metrics[0]:.6f} m, peak error={metrics[1]:.6f} m, '
              f'hold RMSE={metrics[2]:.6f} m')

    # Existing state plots and animation show the assisted run.
    s = results['Force enabled']['states']

    # All plotted states are simulated truth on the same t_eval grid.
    # Angles are derived from the relative position.
    xq_plot, yq_plot, zq_plot = s[:, 0:3].T
    xp_plot, yp_plot, zp_plot = s[:, 6:9].T
    th_plot = np.arccos(np.clip((zq_plot - zp_plot) / l, -1.0, 1.0))
    phi_plot = np.arctan2(yp_plot - yq_plot, xp_plot - xq_plot)
    # Azimuth is undefined when the tether is vertical.
    phi_plot[np.hypot(xp_plot - xq_plot, yp_plot - yq_plot) < 1e-8] = np.nan

    fig = plt.figure('Drone Velocity Control Simulation — Force Enabled', figsize=(12, 9))
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(t, s[:, 0:3])
    for i, label in enumerate(['x', 'y', 'z']):
        ax1.plot(t, position_ref[:, i], '--', color=f'C{i}', label=f'{label} reference')
    ax1.legend(['xq', 'yq', 'zq', 'x reference', 'y reference', 'z reference'])
    ax1.set(xlabel='Time (s)', ylabel='Position (m)', title='Drone Position')
    ax1.grid(True)

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(t, s[:, 6:9])
    for i in range(3):
        ax2.plot(t, payload_ref[:, i], '--', color=f'C{i}')
    ax2.legend(['xp', 'yp', 'zp', 'x reference', 'y reference', 'z reference'])
    ax2.set(xlabel='Time (s)', ylabel='Position (m)', title='Payload Position')
    ax2.grid(True)

    ax_angles = fig.add_subplot(2, 2, 3)
    ax_angles.plot(t, np.rad2deg(th_plot), label='Inclination θ')
    ax_angles.plot(t, np.rad2deg(phi_plot), label='Azimuth φ')
    ax_angles.set(xlabel='Time (s)', ylabel='Angle (deg)', title='Pendulum Angles')
    ax_angles.legend()
    ax_angles.grid(True)

    ax3 = fig.add_subplot(2, 2, 4, projection='3d')
    ax3.plot(s[:, 0], s[:, 1], s[:, 2], 'b-', label='Drone')
    ax3.plot(s[:, 6], s[:, 7], s[:, 8], 'r--', label='Payload (force enabled)')
    s_baseline = results['Force disabled']['states']
    ax3.plot(*s_baseline[:, 6:9].T, color='orange', linestyle='-.', label='Payload (force disabled)')
    ax3.plot(*position_ref.T, 'k:', label='Drone reference')
    ax3.plot(*payload_ref.T, color='gray', linestyle=':', label='Payload reference')
    ax3.legend()
    ax3.set(xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)', title='3D Trajectory')
    ax3.set_aspect('equal')
    ax3.grid(True)
    fig.tight_layout()

    # Compare simulated truth errors on the common t_eval grid (actual − reference).
    fig_error, axes_error = plt.subplots(4, 1, num='Payload Tracking Comparison',
                                         figsize=(10, 10), sharex=True)
    styles = [('Force disabled', 'C0', '--'), ('Force enabled', 'C1', '-')]
    for label, color, linestyle in styles:
        for i, axis in enumerate(['x', 'y', 'z']):
            axes_error[i].plot(t, results[label]['error'][:, i],
                               color=color, linestyle=linestyle, label=label)
            axes_error[i].set_ylabel(f'{axis} error (m)')
        axes_error[3].plot(t, results[label]['error_norm'],
                           color=color, linestyle=linestyle, label=label)
    axes_error[0].set_title('Payload Tracking Error (actual − reference)')
    axes_error[3].set(xlabel='Time (s)', ylabel='3D error norm (m)')
    for ax in axes_error:
        ax.axvline(t_finish, color='gray', linestyle=':', label='Square complete')
        ax.legend()
        ax.grid(True)
    fig_error.tight_layout()

    fig_metrics, ax_metrics = plt.subplots(num='Payload Performance Metrics', figsize=(10, 5))
    metric_names = ['3D position RMSE', 'Peak 3D error', 'Hold 3D RMSE']
    x = np.arange(len(metric_names))
    for i, (label, color, _) in enumerate(styles):
        bars = ax_metrics.bar(x + (i - 0.5) * 0.35, results[label]['metrics'],
                              width=0.35, color=color, label=label)
        ax_metrics.bar_label(bars, fmt='%.4g', padding=3)
    ax_metrics.set(xticks=x, xticklabels=metric_names, ylabel='Position error (m)',
                   title='Payload Performance: Force Disabled vs Enabled')
    ax_metrics.margins(y=0.2)
    ax_metrics.legend()
    ax_metrics.grid(True, axis='y')
    fig_metrics.tight_layout()

    # Applied payload force on the same simulation time grid.
    payload_force = np.array([payload_tracking_force(ti, si) for ti, si in zip(t, s)])
    fig_force, ax_force = plt.subplots(num='Payload Horizontal Force', figsize=(10, 4))
    ax_force.plot(t, payload_force[:, 0], 'r', label='Fx')
    ax_force.plot(t, payload_force[:, 1], 'g', label='Fy')
    ax_force.plot(t, np.linalg.norm(payload_force, axis=1), 'k', label='Magnitude')
    ax_force.plot(t, np.zeros_like(t), 'C0--', label='Force disabled (magnitude)')
    ax_force.axhline(max_payload_force, color='gray', linestyle='--', label='Force limit')
    ax_force.set(xlabel='Time (s)', ylabel='Force (N)', title='Applied Payload Force')
    ax_force.legend()
    ax_force.grid(True)
    fig_force.tight_layout()

    # --- Animation ---
    fig_anim = plt.figure('Drone Animation — Force Enabled')
    ax_anim = fig_anim.add_subplot(111, projection='3d')
    line, = ax_anim.plot([], [], [], 'o-', lw=2)
    ax_anim.set(xlim=(s[:, [0, 6]].min() - 0.3, s[:, [0, 6]].max() + 0.3),
                ylim=(s[:, [1, 7]].min() - 0.3, s[:, [1, 7]].max() + 0.3),
                zlim=(s[:, [2, 8]].min() - 0.3, s[:, [2, 8]].max() + 0.3),
                xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
    ax_anim.set_aspect('equal')

    def init():
        line.set_data([], [])
        line.set_3d_properties([])
        return line,

    def update(i):
        # Line between Drone (0:3) and Payload (6:9)
        line.set_data([s[i, 0], s[i, 6]], [s[i, 1], s[i, 7]])
        line.set_3d_properties([s[i, 2], s[i, 8]])
        ax_anim.set_title(f'Time: {t[i]:.2f} s')
        return line,

    # Keep the animation alive until the plot window closes.
    ani = animation.FuncAnimation(fig_anim, update, frames=range(0, len(t), 2),
                                  init_func=init, blit=False, interval=40)
    plt.show()
    return ani


if __name__ == '__main__':
    main()
