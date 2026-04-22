function [phi_pos, phi_vel , phi_acc , t_out] = simulate_flying_object(dt, k, b, J)
    tspan = 0:dt:20; 
    y0 = [0; 0]; 
    
    % Ορισμός της εισόδου
    u = @(t) 0.25 * sin(0.5 * pi * t);
    
    % Πίνακες State-Space
    A = [0,      1; 
         0,   -k/J];
    B = [0; 
         b/J];
         
    C = [1 , 0;
         0 , 1;
         0 , -k/J];
    D = [0; 
         0; 
         b/J];
         
    % y_dot = A*y + B*u(t)
    sys_odes = @(t, y) A * y + B * u(t);
    
    % Λύση με ODE45
    options = odeset('MaxStep', dt); 
    [t_out, y_out] = ode45(sys_odes, tspan, y0, options);
    
    
    u_vals = u(t_out); 
    

    z = C * y_out' + D * u_vals'; 
    
    % 3. Διαχωρισμός των αποτελεσμάτων και μετατροπή ξανά σε στήλες (Nx1)
    phi_pos = z(1, :)';
    phi_vel = z(2, :)';
    phi_acc = z(3, :)';
    
end