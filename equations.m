function dX = equations(t, X, J, k_real, b_real, gamma1, gamma2, d_func)
    
    u = 0.25 * sin(0.5 * pi * t);
    d = d_func(t); 
    
    % Ανάγνωση των 6 καταστάσεων
    phi         = X(1);
    phi_dot     = X(2);
    phi_hat     = X(3);
    phi_dot_hat = X(4);
    k_hat       = X(5);
    b_hat       = X(6);
    
    % Πραγματικό Σύστημα 
    phi_ddot = (1/J) * (-k_real * phi_dot + b_real * u + d);
    
    % Μοντέλο 
    phi_ddot_hat = (1/J) * (-k_hat * phi_dot_hat + b_hat * u);
    
    % Μέθοδος Κλίσης
    k_hat_dot = gamma1 * (-J*phi_ddot*phi_dot - k_hat*phi_dot^2 + b_hat*u*phi_dot);
    b_hat_dot = -gamma2 * (-J*phi_ddot*u - k_hat*phi_dot*u + b_hat*u^2);
    
    % Παράγωγοι
    dX = zeros(6,1);
    % x1 = φ , x2 = φ_dot

    dX(1) = phi_dot; %x1_dot = x2 = φ_dot
    
    dX(2) = phi_ddot;  
    
    % x1 = φ_hat , x2 = φ_hat_dot
    dX(3) = phi_dot_hat;   
    dX(4) = phi_ddot_hat;

    dX(5) = k_hat_dot;     
    dX(6) = b_hat_dot;     
end