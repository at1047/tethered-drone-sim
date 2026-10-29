import sympy as sp

# 1. Define symbolic variables
print("Defining symbols...")
xq, yq, zq, xp, yp, zp = sp.symbols('xq yq zq xp yp zp', real=True)
dxq, dyq, dzq, dxp, dyp, dzp = sp.symbols('dxq dyq dzq dxp dyp dzp', real=True)
ddxq, ddyq, ddzq, ddxp, ddyp, ddzp = sp.symbols('ddxq ddyq ddzq ddxp ddyp ddzp', real=True)
f1, f2, f3 = sp.symbols('f1 f2 f3', real=True)
m_q, m_p, l, g, la = sp.symbols('m_q m_p l g la', real=True)

q = sp.Matrix([xq, yq, zq, xp, yp, zp])
dq = sp.Matrix([dxq, dyq, dzq, dxp, dyp, dzp])
ddq = sp.Matrix([ddxq, ddyq, ddzq, ddxp, ddyp, ddzp])
u = sp.Matrix([f1, f2, f3])
q_act = sp.Matrix([xq, yq, zq])

# 2. Energies
print("Calculating energies...")
T = 0.5 * m_q * (dxq**2 + dyq**2 + dzq**2) + 0.5 * m_p * (dxp**2 + dyp**2 + dzp**2)
V = m_q * g * zq + m_p * g * zp

# 3. Constraint: (xq-xp)^2 + (yq-yp)^2 + (zq-zp)^2 - l^2 = 0
c = (xq - xp)**2 + (yq - yp)**2 + (zq - zp)**2 - l**2

# 4. Dynamics Matrices (D, C, G, B)
print("Deriving D matrix...")
D = sp.Matrix([[sp.diff(sp.diff(T, dq[i]), dq[j]) for j in range(len(q))] for i in range(len(q))]).simplify()

print("Deriving C matrix (this may take a moment)...")
n = len(q)
C = sp.zeros(n, n)
for k in range(n):
    for j in range(n):
        for i in range(n):
            C[k, j] += 0.5 * (sp.diff(D[k, j], q[i]) + sp.diff(D[k, i], q[j]) - sp.diff(D[i, j], q[k])) * dq[i]
C = sp.simplify(C)

print("Deriving G and B matrices...")
G = sp.Matrix([sp.diff(V, q[i]) for i in range(len(q))]).simplify()
B = q_act.jacobian(q).T

# 5. Constraint Jacobian and Acceleration-level constraint (ddc = 0)
print("Calculating constraint acceleration...")
A = sp.Matrix([c]).jacobian(q) # 1x6
# dA/dt * dq
Adot = A.jacobian(q) * dq
ddc = (A * ddq)[0] + (Adot.dot(dq))

# 6. EOM with Lagrange Multiplier
# D*ddq = B*u - C*dq - G + A.T*la
print("Solving for ddq in terms of lambda...")
rhs_eom = B * u - C * dq - G + A.T * la
ddq_sol = D.LUsolve(rhs_eom)

# 7. Solve for lambda (la)
# Substitute ddq_sol into ddc = 0
print("Solving for lambda...")
eq_for_la = ddc.subs([(ddq[i], ddq_sol[i]) for i in range(len(q))])
la_sol = sp.solve(eq_for_la, la)[0]

# 8. Substitute lambda back into ddq_sol
print("Calculating final equations of motion...")
ddq_final = sp.simplify(ddq_sol.subs(la, la_sol))

# 9. Print the 6 Equations of Motion
print()
print("="*50)
print("FINAL 6 EQUATIONS OF MOTION (ddq_final)")
print("="*50)
labels = ['ddxq', 'ddyq', 'ddzq', 'ddxp', 'ddyp', 'ddzp']
for i in range(len(labels)):
    print(f"{labels[i]} = {ddq_final[i]}")
print("="*50)
