clear; clc; close all;

% --- 1. Παράμετροι Συστήματος ---
g = 9.81;     % Επιτάχυνση βαρύτητας
l = 1.5;      % Πραγματικό μήκος (επιλογή στο [1, 2.5])
c = 0.5;      % Πραγματική απόσβεση (επιλογή στο [0.2, 0.8])

% --- 2. Παράμετροι Εκτιμητή ---
gamma1 = 10;  % Κέρδος προσαρμογής για το a (a = g/l)
gamma2 = 5;   % Κέρδος προσαρμογής για το c
L1 = 10;      % Κέρδος σφάλματος θέσης (Observer)
L2 = 10;      % Κέρδος σφάλματος ταχύτητας (Observer)

tspan = [0 150];
% Αρχικές συνθήκες: [theta, theta_dot, theta_hat, theta_dot_hat, a_hat, c_hat]
% Έστω ότι ξεκινάμε με μια αρχική εκτίμηση l_hat = 1.0 -> a_hat = 9.81
X0 = [0; 0; 0; 0; g/1.0; 0.2]; 

% =========================================================================
% ΣΥΝΑΡΤΗΣΗ ΔΥΝΑΜΙΚΟΥ ΣΥΣΤΗΜΑΤΟΣ (ΕΝΣΩΜΑΤΩΜΕΝΗ)
% =========================================================================
function dX = pendulum_sys(t, X, g, l_real, c_real, gamma1, gamma2, L1, L2, eta_0)
    % Είσοδος
    u = 0.5 * sin(t);
    
    % Θόρυβος Μετρήσεων
    eta = eta_0 * sin(20 * pi * t);
    
    % Καταστάσεις
    theta_true     = X(1);
    theta_dot_true = X(2);
    theta_hat      = X(3);
    theta_dot_hat  = X(4);
    a_hat          = X(5); 
    c_hat          = X(6);
    
    % Μετρήσεις (προσθήκη θορύβου)
    theta_m     = theta_true + eta;
    theta_dot_m = theta_dot_true + eta;
    u_m         = u + eta;
    
    % Πραγματικό Σύστημα
    theta_ddot_true = -(g/l_real) * sin(theta_true) - c_real * theta_dot_true + u;
    
    % Σφάλματα βασισμένα στις μετρήσεις
    e_theta     = theta_m - theta_hat;
    e_theta_dot = theta_dot_m - theta_dot_hat;
    
    % Εκτιμητής ταχύτητας & επιτάχυνσης 
    theta_ddot_hat = -a_hat * sin(theta_m) - c_hat * theta_dot_m + u_m;
    
    % Νόμοι Προσαρμογής Lyapunov
    a_hat_dot = -gamma1 * e_theta_dot * sin(theta_m);
    c_hat_dot = -gamma2 * e_theta_dot * theta_dot_m;
    
    % Παράγωγοι
    dX = zeros(6,1);
    dX(1) = theta_dot_true;
    dX(2) = theta_ddot_true;
    dX(3) = theta_dot_hat + L1 * e_theta;       
    dX(4) = theta_ddot_hat + L2 * e_theta_dot;  
    dX(5) = a_hat_dot;
    dX(6) = c_hat_dot;
end

% =========================================================================
% ΘΕΜΑ 2.α: ΧΩΡΙΣ ΘΟΡΥΒΟ
% =========================================================================
eta_a = 0;
[t_a, X_a] = ode45(@(t,X) pendulum_sys(t, X, g, l, c, gamma1, gamma2, L1, L2, eta_a), tspan, X0);

theta_a     = X_a(:,1);
theta_hat_a = X_a(:,3);
e_theta_a   = theta_a - theta_hat_a;

% Μετατροπή a_hat πίσω σε l_hat
l_hat_a = g ./ X_a(:,5);
c_hat_a = X_a(:,6);

e_l_a = l_hat_a - l;
e_c_a = c_hat_a - c;

figure('Name', 'Θέμα 2.α: Εκτίμηση με Μέθοδο Lyapunov (Χωρίς Θόρυβο)');
subplot(4,1,1); plot(t_a, theta_a, 'b', t_a, theta_hat_a, 'r--'); title('Γωνία \theta(t) vs \theta_{hat}(t)'); grid on;
subplot(4,1,2); plot(t_a, e_theta_a, 'g'); title('Σφάλμα Γωνίας e_\theta(t)'); grid on;
subplot(4,1,3); plot(t_a, e_l_a, 'k'); title('Σφάλμα Παραμέτρου e_l(t)'); grid on;
subplot(4,1,4); plot(t_a, e_c_a, 'm'); title('Σφάλμα Παραμέτρου e_c(t)'); xlabel('Χρόνος (s)'); grid on;

% =========================================================================
% ΘΕΜΑ 2.β: ΜΕ ΘΟΡΥΒΟ ΚΑΙ ΜΕΛΕΤΗ ΣΦΑΛΜΑΤΟΣ
% =========================================================================
% Δημιουργία γραφήματος για μία ενδεικτική τιμή θορύβου (π.χ. eta_0 = 0.05)
eta_b = 0.05;
[t_b, X_b] = ode45(@(t,X) pendulum_sys(t, X, g, l, c, gamma1, gamma2, L1, L2, eta_b), tspan, X0);

l_hat_b = g ./ X_b(:,5);
e_l_b = l_hat_b - l;
e_c_b = X_b(:,6) - c;

figure('Name', 'Θέμα 2.β: Ενδεικτικό Γράφημα με Θόρυβο (\eta_0 = 0.05)');
subplot(2,1,1); plot(t_b, e_l_b, 'k'); title('Σφάλμα e_l(t) παρουσία θορύβου'); grid on;
subplot(2,1,2); plot(t_b, e_c_b, 'm'); title('Σφάλμα e_c(t) παρουσία θορύβου'); xlabel('Χρόνος (s)'); grid on;

% Μελέτη σφάλματος ως προς το πλάτος του θορύβου
eta_vals = 0:0.01:0.1;
err_l_steady = zeros(length(eta_vals), 1);
err_c_steady = zeros(length(eta_vals), 1);

for i = 1:length(eta_vals)
    [~, X_temp] = ode45(@(t,X) pendulum_sys(t, X, g, l, c, gamma1, gamma2, L1, L2, eta_vals(i)), tspan, X0);
    
    l_temp = g ./ X_temp(:, 5);
    c_temp = X_temp(:, 6);
    
    % Υπολογισμός Mean Absolute Error στα τελευταία 30 sec (steady state)
    idx_steady = round(0.8 * length(X_temp)) : length(X_temp);
    err_l_steady(i) = mean(abs(l_temp(idx_steady) - l));
    err_c_steady(i) = mean(abs(c_temp(idx_steady) - c));
end

figure('Name', 'Θέμα 2.β: Σφάλμα vs Πλάτος Θορύβου');
plot(eta_vals, err_l_steady, 'ko-', 'LineWidth', 1.5); hold on;
plot(eta_vals, err_c_steady, 'ms-', 'LineWidth', 1.5);
title('Μέσο Απόλυτο Σφάλμα Εκτίμησης vs Πλάτος Θορύβου \eta_0');
xlabel('\eta_0'); ylabel('Σφάλμα (Μόνιμη Κατάσταση)');
legend('Σφάλμα e_l', 'Σφάλμα e_c'); grid on;