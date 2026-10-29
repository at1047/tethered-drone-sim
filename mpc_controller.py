import numpy as np
from scipy.optimize import minimize

class MPCController:
    def __init__(self, m_q, m_p, l, g):
        self.m_q = m_q
        self.m_p = m_p
        self.l = l
        self.g = g
        
        # MPC Parameters
        self.dt = 0.1       # Prediction step size
        self.N = 5          # Prediction horizon
        self.num_inputs = 3 # f1, f2, f3
        
        # Weights
        self.W_pos = np.array([10.0, 10.0, 10.0]) # Drone position
        self.W_vel = np.array([1.0, 1.0, 1.0])    # Drone velocity
        self.W_swing_y = 50.0                     # Lateral swing penalty (CRITICAL)
        self.W_swing_x = 1.0                      # Longitudinal swing penalty (LOW)
        self.W_control = np.array([0.01, 0.01, 0.01]) # Control effort

    def dynamics(self, s, u):
        """Discrete-time dynamics for prediction."""
        xq, yq, zq, dxq, dyq, dzq, xp, yp, zp, dxp, dyp, dzp = s
        f1, f2, f3 = u
        
        # Use the derived physics (simplified/copied for the controller's internal model)
        # These must match the lagrangian_drone_dynamics derivation
        m_p, m_q, g, l = self.m_p, self.m_q, self.g, self.l
        
        # Constraint: (xq-xp)^2 + (yq-yp)^2 + (zq-zp)^2 - l^2 = 0
        # We use a simplified version for the MPC internal model if needed, 
        # but here we'll use the full expressions for better accuracy.
        
        # Re-using the derived expressions (shortened for brevity in the internal model)
        den = (m_p*(xp - xq)**2 + m_p*(yp - yq)**2 + m_p*(zp - zq)**2 + m_q*(xp - xq)**2 + m_q*(yp - yq)**2 + m_q*(zp - zq)**2)
        term_common = (dxp**2*m_q - 2.0*dxp*dxq*m_q + dxq**2*m_q + dyp**2*m_q - 2.0*dyp*dyq*m_q + dyq**2*m_q + dzp**2*m_q - 2.0*dzp*dzq*m_q + dzq**2*m_q - f1*xp + f1*xq - f2*yp + f2*yq - f3*zp + f3*zq)
        
        ddxq = (f1*(den/m_p) + (xp - xq)*term_common) * (m_p / (m_q * den))
        ddyq = (f2*(den/m_p) + (yp - yq)*term_common) * (m_p / (m_q * den))
        ddzq = (m_p*(zp - zq)*term_common + (f3 - g*m_q)*den) / (m_q * den)
        
        ddxp = - (xp - xq)*term_common / den
        ddyp = - (yp - yq)*term_common / den
        ddzp = (-g*den - (zp - zq)*term_common) / den

        # Euler integration
        ds = np.array([dxq, dyq, dzq, ddxq, ddyq, ddzq, dxp, dyp, dzp, ddxp, ddyp, ddzp])
        return s + ds * self.dt

    def cost_function(self, u_flat, s0, s_ref):
        u_seq = u_flat.reshape((self.N, self.num_inputs))
        s = s0.copy()
        total_cost = 0
        
        for i in range(self.N):
            u = u_seq[i]
            s = self.dynamics(s, u)
            
            # Position error
            total_cost += np.sum(self.W_pos * (s[0:3] - s_ref[0:3])**2)
            # Velocity error
            total_cost += np.sum(self.W_vel * s[3:6]**2)
            
            # Control effort
            total_cost += np.sum(self.W_control * u**2)
            
            # Penalize swing only at the end of the horizon (Terminal Cost)
            if i == self.N - 1:
                # Lateral swing (Y axis) - We want yp approx yq
                # Indices: yq=1, dyq=4, yp=7, dyp=10
                total_cost += self.W_swing_y * (s[7] - s[1])**2 # (yp - yq)^2
                total_cost += self.W_swing_y * (s[10] - s[4])**2 # (dyp - dyq)^2
                
                # Longitudinal swing (X axis) - allow some swing but penalize at end
                # Indices: xq=0, dxq=3, xp=6, dxp=9
                total_cost += self.W_swing_x * (s[6] - s[0])**2
                total_cost += self.W_swing_x * (s[9] - s[3])**2
            
        return total_cost

    def compute_force(self, t, s, s_ref):
        # Initial guess: Hover force
        u0 = np.tile([0, 0, (self.m_q + self.m_p) * self.g], self.N)
        
        # Bounds (Optional, but good for stability)
        bounds = [(-20, 20), (-20, 20), (0, 40)] * self.N
        
        res = minimize(self.cost_function, u0, args=(s, s_ref), 
                       method='SLSQP', bounds=bounds,
                       options={'maxiter': 10, 'disp': False})
        
        u_optimal = res.x.reshape((self.N, self.num_inputs))
        return u_optimal[0] # Return the first control action
