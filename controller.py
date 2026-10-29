import numpy as np

class PDController:
    def __init__(self, m_q, m_p, l, g):
        self.m_q = m_q
        self.m_p = m_p
        self.l = l
        self.g = g
        
        # PD Gains for Drone Position
        self.kp_q = 15.0
        self.kd_q = 10.0
        
        # PD Gains for Payload Swing (Full State Feedback)
        # These help damp the oscillations of the payload
        self.kp_p = 5.0
        self.kd_p = 2.0

    def compute_force(self, t, s, s_ref):
        """
        Computes control forces f1, f2, f3 based on full state feedback.
        s: current state [xq, yq, zq, dxq, dyq, dzq, xp, yp, zp, dxp, dyp, dzp]
        s_ref: desired state (same format, usually velocities are 0)
        """
        xq, yq, zq = s[0:3]
        dxq, dyq, dzq = s[3:6]
        xp, yp, zp = s[6:9]
        dxp, dyp, dzp = s[9:12]
        
        xq_ref, yq_ref, zq_ref = s_ref[0:3]
        dxq_ref, dyq_ref, dzq_ref = s_ref[3:6]
        xp_ref, yp_ref, zp_ref = s_ref[6:9]
        dxp_ref, dyp_ref, dzp_ref = s_ref[9:12]

        # 1. Horizontal Control (X and Y)
        # We combine drone position error and payload position error
        # Effectively, we pull the drone towards the target while also pulling 
        # it in a way that damps the payload swing.
        
        ex = (xq_ref - xq)
        edx = (dxq_ref - dxq)
        exp = (xp_ref - xp)
        edxp = (dxp_ref - dxp)
        
        ey = (yq_ref - yq)
        edy = (dyq_ref - dyq)
        eyp = (yp_ref - yp)
        edyp = (dyp_ref - dyp)

        # Feedforward: none for horizontal
        f1 = self.m_q * (self.kp_q * ex + self.kd_q * edx + self.kp_p * exp + self.kd_p * edxp)
        f2 = self.m_q * (self.kp_q * ey + self.kd_q * edy + self.kp_p * eyp + self.kd_p * edyp)

        # 2. Vertical Control (Z)
        ez = (zq_ref - zq)
        edz = (dzq_ref - dzq)
        ezp = (zp_ref - zp)
        edzp = (dzp_ref - dzp)
        
        # Feedforward: gravity compensation for both masses
        f_gravity = (self.m_q + self.m_p) * self.g
        f3 = self.m_q * (self.kp_q * ez + self.kd_q * edz + self.kp_p * ezp + self.kd_p * edzp) + f_gravity

        return f1, f2, f3
