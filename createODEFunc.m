function createODEFunc(W1, W2, R)
    syms x y z psi_var theta phi alpha_var beta_var
    syms x_dot y_dot z_dot psi_dot theta_dot phi_dot alpha_dot beta_dot t_f
    U = sym('U',[4 1]);
    syms l M m_l I_psi I_theta I_phi I g
    consts = [l M m_l I_psi I_theta I_phi I g];
    syms lam [16 1]
    
    zeta = [x; y; z];
    zeta_dot = [x_dot; y_dot; z_dot];
    eta = [psi_var; theta; phi];
    eta_dot = [psi_dot; theta_dot; phi_dot];
    mu = [alpha_var; beta_var];
    mu_dot = [alpha_dot; beta_dot];
    % U = [U_1; U_2; U_3; U_4];
    
    q = [zeta; eta; mu];
    q_dot = [zeta_dot; eta_dot; mu_dot];
    
    c_alpha = cos(alpha_var);
    s_alpha = sin(alpha_var);
    c_beta = cos(beta_var);
    s_beta = sin(beta_var);
    c_psi = cos(psi_var);
    s_psi = sin(psi_var);
    c_theta = cos(theta);
    s_theta = sin(theta);
    c_phi = cos(phi);
    s_phi = sin(phi);
    m_11 = M + m_l;
    m_22 = m_11;
    m_33 = m_22;
    m_17 = m_l*l*cos(alpha_var)*cos(beta_var);
    m_71 = m_17;
    m_18 = -m_l*l*sin(alpha_var)*sin(beta_var);
    m_81 = m_18;
    m_27 = m_l*l*cos(alpha_var)*sin(beta_var);
    m_72 = m_27;
    m_28 = m_l*l*sin(alpha_var)*cos(beta_var);
    m_82 = m_28;
    m_37 = m_l*l*sin(alpha_var);
    m_73 = m_37;
    m_44 = I_psi*s_theta*s_theta + c_theta*c_theta*(I_theta*s_phi*s_phi + I_phi*c_phi*c_phi);
    m_45 = (I_theta - I_phi)*c_theta*s_phi*c_phi;
    m_54 = m_45;
    m_55 = I_theta*c_phi*c_phi + I_phi*s_phi*s_phi;
    m_77 = m_l*l*l + I;
    m_88 = m_l*l*l*s_alpha*s_alpha + I;
    
    M_q = [   m_11,       0,      0,      0,              0,      0,              m_17,       m_18;
            0,          m_22,   0,      0,              0,      0,              m_27,       m_28;
            0,          0,      m_33,   0,              0,      0,              m_37,       0;
            0,          0,      0,      m_44,           m_45,   -I_psi*s_theta, 0,          0;
            0,          0,      0,      m_54,           m_55,   0,              0,          0;
            0,          0,      0,      -I_psi*s_theta, 0,      I_psi,          0,          0;
            m_71,       m_72,   m_73,   0,              0,      0,              m_77,       0;
            m_81,       m_82,   0,      0,              0,      0,              0,          m_88    ];
    
    
    c_17 = -m_l*l*(c_alpha*s_beta*beta_dot + s_alpha*c_beta*alpha_dot);
    c_18 = -m_l*l*(c_alpha*s_beta*alpha_dot + s_alpha*c_beta*beta_dot);
    c_27 = m_l*l*(c_alpha*c_beta*beta_dot - s_alpha*s_beta*alpha_dot);
    c_28 = m_l*l*(c_alpha*c_beta*alpha_dot - s_alpha*s_beta*beta_dot);
    c_44 = I_psi*theta_dot*s_theta*c_theta - (I_theta + I_phi)*theta_dot*s_theta*c_theta*s_phi*s_phi + ...
            (I_theta - I_phi)*phi_dot*c_theta*c_theta*s_phi*c_phi;
    c_45 = I_psi*psi_dot*s_theta*c_theta - (I_theta - I_phi)*(theta_dot*s_theta*c_phi*s_phi + ...
            phi_dot*c_theta*s_phi*s_phi) - (I_theta + I_phi)*(psi_dot*s_theta*c_theta*c_phi*c_phi - ...
            phi_dot*c_theta*c_phi*c_phi);
    c_46 = -(I_psi*theta_dot*c_theta - (I_theta - I_phi)*(psi_dot*c_theta*c_theta*s_phi*c_phi));
    c_54 = psi_dot*s_theta*c_theta*(-I_psi + I_theta*s_phi*s_phi + I_phi*c_phi*c_phi);
    c_55 = -(I_theta - I_phi)*phi_dot*s_phi*c_phi;
    c_56 = I_psi*psi_dot*c_theta + (I_theta - I_phi)*(-theta_dot*s_theta*c_phi + psi_dot*c_theta*c_phi*c_phi - ...
            psi_dot*c_theta*s_phi*s_phi);
    c_64 = -(I_theta - I_phi)*psi_dot*c_theta*c_theta*s_phi*c_phi;
    c_65 = -I_psi*psi_dot*c_theta + (I_theta - I_phi)*(theta_dot*s_phi*c_phi + psi_dot*c_theta*s_phi*s_phi - ...
            psi_dot*c_theta*c_phi*c_phi);
    
    C_q = [   0,      0,      0,      0,      0,      0,      c_17,                               c_18;
            0,      0,      0,      0,      0,      0,      c_27,                               c_28;
            0,      0,      0,      0,      0,      0,      m_l*l*c_alpha*alpha_dot,            0;
            0,      0,      0,      c_44,   c_45,   c_46,   0,                                  0;
            0,      0,      0,      c_54,   c_55,   c_56,   0,                                  0;
            0,      0,      0,      c_64,   c_65,   0,      0,                                  0;
            0,      0,      0,      0,      0,      0,      0,                                  -m_l*l*l*s_alpha*c_alpha*beta_dot;
            0,      0,      0,      0,      0,      0,      m_l*l*l*s_alpha*c_alpha*beta_dot,   m_l*l*l*s_alpha*c_alpha*alpha_dot   ];
    
    G_q = [0; 0;  (M+m_l)*g;  0; 0; 0; m_l*l*g*s_alpha; 0];
    
    b_q = [   s_alpha*s_psi + c_phi*c_psi*s_theta,    0,  0,  0;
            c_phi*s_theta*s_psi - c_psi*s_phi,      0,  0,  0;
            c_theta*c_phi,                          0,  0,  0;
            0,                                      1,  0,  0;
            0,                                      0,  1,  0;
            0,                                      0,  0,  1;
            0,                                      0,  0,  0;
            0,                                      0,  0,  0   ];
    
    N = -inv(M_q)*(C_q*q_dot + G_q);
    Z = inv(M_q);
    F_1 = q_dot;
    F_2 = N + Z*b_q*U;
    
    % W1 = eye(8);
    % W2 = eye(8);
    % R = 0.1*eye(8);
    
    F = [F_1; F_2];
    X = transpose([x y z psi_var theta phi alpha_var beta_var x_dot y_dot z_dot psi_dot theta_dot phi_dot alpha_dot beta_dot]);
    
    W = [W1 zeros(8,8);zeros(8,8) W2];
    L = 0.5*transpose(X)*W*X+0.5*transpose(U)*R*U;
    H = L + transpose(lam)*F;
    eq1 = jacobian(H, U) == 0;

    U_sol = solve(eq1, U);

    U_star = [U_sol.U1; U_sol.U2; U_sol.U3; U_sol.U4];

    % U_star = -transpose(inv(R))*[zeros(4, 8) transpose(Z*b_q)]*lam;

    U_star_formula = matlabFunction(U_star, 'File', 'UStarFormula', 'Vars', {X, lam, consts});
    
    lamdot = -transpose(jacobian(H,X));
    lamdot = subs(lamdot,U,U_star);
    X_dot = subs(F,U,U_star);
    syms t
    X_aug = [X; lam];
    
    out = [X_dot; lamdot];
    bvpode = matlabFunction(out, 'File', 'bvpode', 'Vars', {t, X_aug, consts});

end