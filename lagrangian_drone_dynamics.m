clear
clc

% syms xq(t) yq(t) zq(t) xp(t) yp(t) zp(t) 
syms th phi real
syms f1 f2 f3 real
syms xq yq zq xp yp zp real 
syms dxq dyq dzq dxp dyp dzp real
syms ddxq ddyq ddzq ddxp ddyp ddzp real
syms dth ddth dphi ddphi real
syms m_q m_p l g real
syms la real
% 
% dxq = diff(xq, t);
% dyq = diff(yq, t);
% dzq = diff(zq, t);
% dxp = diff(xp, t);
% dyp = diff(yp, t);
% dzp = diff(zp, t);
% 
% ddxq = diff(dxq, t);
% ddyq = diff(dyq, t);
% ddzq = diff(dzq, t);
% ddxp = diff(dxp, t);
% ddyp = diff(dyp, t);
% ddzp = diff(dzp, t);

% xp = l*sin(th)*cos(phi)+xq;
% yp = l*sin(th)*sin(phi)+yq;
% zp = -l*cos(th)+zq;
% 
% dxp = dxq + l*dth*cos(th)*cos(phi) - l*dphi*sin(th)*sin(phi);
% dyp = dyq + l*dth*cos(th)*sin(phi) + l*dphi*sin(th)*cos(phi);
% dzp = dzq - l*dth*sin(th);

% T = 0.5*m_q*(dxq^2 + dyq^2 + dzq^2) + 0.5*m_p*(dxp^2 + dyp^2 + dzp^2);
% V = m_q*g*zq + m_p*g*zp;

T = 0.5*m_q*(dxq^2 + dyq^2 + dzq^2) + 0.5*m_p*(dxp^2 + dyp^2 + dzp^2);
V = m_q*g*zq + m_p*g*zp;

c = (xq-xp)^2 + (yq-yp)^2 + (zq-zp)^2 - l^2;

L = T - V;

q = [xq; yq; zq; xp; yp; zp];
dq = [dxq; dyq; dzq; dxp; dyp; dzp];
ddq = [ddxq; ddyq; ddzq; ddxp; ddyp; ddzp];
q_act = [xq; yq; zq];
u = [f1;f2;f3];

D = simplify( jacobian(jacobian(T,dq), dq) ) ;
for k=1:length(q)
    for j=1:length(q)
        C(k,j) = sym(0) ;
        for i=1:length(q)
            C(k,j) = C(k,j) + 1/2 * ( diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)) ) * dq(i) ;
        end
    end
end
C = simplify(C) ;
G = simplify( jacobian(V,q) )' ;
B = jacobian(q_act, q)' ;

% c = (xq-xp)^2 + (yq-yp)^2 + (zq-zp)^2 - l^2
A = jacobian(c, q); % This is dc/dq (1x6 row vector)

% To get the acceleration-level constraint (ddc = 0)
% We need: A*ddq + (partial A / partial q)*dq^2 = 0
% A simpler way in MATLAB symbolic:
syms t real
syms xqt(t) yqt(t) zqt(t) xpt(t) ypt(t) zpt(t) 
q_t = [xqt(t); yqt(t); zqt(t); xpt(t); ypt(t); zpt(t)];
c_t = subs(c, q, q_t);
ddc = diff(c_t, t, 2); % This is the d^2/dt^2 of the constraint

% Clean it back to symbols
ddc = subs(ddc, diff(q_t, t, 2), ddq);
ddc = subs(ddc, diff(q_t, t), dq);
ddc = subs(ddc, q_t, q);

% 1. Express ddq in terms of lambda from the Equations of Motion
% D*ddq = B*u - C*dq - G + A'*la
rhs_eom = B*u - C*dq - G + A.'*la;
ddq_sol = D \ rhs_eom; % This is ddq expressed with 'la' inside it

% 2. Substitute ddq_sol into the constraint acceleration equation
% We want to solve: ddc(ddq_sol) = 0
eq_for_lambda = subs(ddc, ddq, ddq_sol);

% 3. Solve for lambda (la)
la_sol = solve(eq_for_lambda, la);

% 4. Substitute lambda back into ddq_sol to get pure accelerations
ddq_final = simplify(subs(ddq_sol, la, la_sol));

% --- Export to optimized MATLAB functions ---
vars = {xq, yq, zq, xp, yp, zp, dxq, dyq, dzq, dxp, dyp, dzp, f1, f2, f3, m_q, m_p, l, g};
matlabFunction(ddq_final, 'File', 'auto_ddq_final', 'Vars', vars);