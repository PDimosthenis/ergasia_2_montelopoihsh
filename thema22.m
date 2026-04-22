clear; clc; close all;

% =========================================================================
% 1. ΠΑΡΑΜΕΤΡΟΙ ΠΡΑΓΜΑΤΙΚΟΥ ΣΥΣΤΗΜΑΤΟΣ
% =========================================================================
g = 9.81;        % Επιτάχυνση βαρύτητας
l_true = 1.5;    % Πραγματικό μήκος l (επιλογή [1, 2.5])
c_true = 0.5;    % Πραγματική απόσβεση c (επιλογή [0.2, 0.8])

% Οι πραγματικές "άγνωστες" παράμετροι (όπως τις ορίσαμε στη θεωρία)
theta1_star = g / l_true; 
theta2_star = c_true;

% =========================================================================
% 2. ΠΑΡΑΜΕΤΡΟΙ ΕΚΤΙΜΗΤΗ & LYAPUNOV
% =========================================================================
gamma1 = 15;     % Κέρδος προσαρμογής για το theta1 (μήκος)
gamma2 = 5;      % Κέρδος προσαρμογής για το theta2 (απόσβεση)
am = 10;         % Κέρδος Series-Parallel (για το σφάλμα ταχύτητας)
am1 = 10;        % Κέρδος διόρθωσης για το σφάλμα θέσης (προαιρετικό για το θ_hat)

tspan = [0 100]; % Χρόνος προσομοίωσης

% Αρχικές συνθήκες: [theta, theta_dot, theta_hat, theta_dot_hat, theta1_hat, theta2_hat]
% Έστω ότι ξεκινάμε με λάθος εκτίμηση l=1.0 και c=0.2
X0 = [0; 0; 0; 0; g/1.0; 0.2]; 

% =========================================================================
% 3. ΠΡΟΣΟΜΟΙΩΣΗ ΘΕΜΑ 2(α) - ΧΩΡΙΣ ΘΟΡΥΒΟ
% =========================================================================
noise_amp = 0; % Μηδενικός θόρυβος για το 2α
[t, X] = ode45(@(t,X) pendulum_lyapunov(t, X, theta1_star, theta2_star, gamma1, gamma2, am, am1, noise_amp), tspan, X0);

% Εξαγωγή δεδομένων
theta     = X(:,1);
theta_hat = X(:,3);
th1_hat   = X(:,5);
th2_hat   = X(:,6);

% Υπολογισμός πραγματικών μεγεθών και σφαλμάτων
l_hat = 9.81 ./ th1_hat; 
e_theta = theta - theta_hat;
e_l = l_hat - l_true;
e_c = th2_hat - c_true;

% --- Figure 1: Γραφικές Παραστάσεις 2(α) ---
figure('Name', 'Θέμα 2.α: Εκτίμηση Παραμέτρων (Χωρίς Θόρυβο)', 'Position', [100, 100, 800, 800]);
subplot(4,1,1);
plot(t, theta, 'b', 'LineWidth', 1.5); hold on;
plot(t, theta_hat, 'r--', 'LineWidth', 1.5);
title('Γωνία \theta(t) vs Εκτιμώμενη \theta_{hat}(t)'); legend('\theta', '\theta_{hat}'); grid on;
subplot(4,1,2);
plot(t, e_theta, 'g', 'LineWidth', 1.5);
title('Σφάλμα Γωνίας e_\theta(t)'); grid on;
subplot(4,1,3);
plot(t, e_l, 'k', 'LineWidth', 1.5);
title('Σφάλμα Εκτίμησης Μήκους e_l(t) (Χωρίς Θόρυβο)'); grid on;
subplot(4,1,4);
plot(t, e_c, 'm', 'LineWidth', 1.5);
title('Σφάλμα Εκτίμησης Απόσβεσης e_c(t) (Χωρίς Θόρυβο)'); xlabel('Χρόνος (sec)'); grid on;

% =========================================================================
% 4. ΠΡΟΣΟΜΟΙΩΣΗ ΘΕΜΑ 2(β) - ΜΕ ΘΟΡΥΒΟ (ΕΝΔΕΙΚΤΙΚΗ ΜΕΛΕΤΗ ΧΡΟΝΟΥ)
% =========================================================================
eta0_example = 0.05; % Ενδεικτικό πλάτος θορύβου
[t_n, X_n] = ode45(@(t,X) pendulum_lyapunov(t, X, theta1_star, theta2_star, gamma1, gamma2, am, am1, eta0_example), tspan, X0);

l_hat_n = 9.81 ./ X_n(:,5);
e_l_n = l_hat_n - l_true;
e_c_n = X_n(:,6) - c_true;

% --- Figure 2: Γραφικές Παραστάσεις 2(β) Χρονικά ---
figure('Name', 'Θέμα 2.β: Επίδραση Θορύβου στο Χρόνο (\eta_0 = 0.05)', 'Position', [150, 150, 800, 500]);
subplot(2,1,1);
plot(t_n, e_l_n, 'k', 'LineWidth', 1.2);
title('Σφάλμα Μήκους e_l(t) με Θόρυβο (\eta_0 = 0.05)'); grid on;
subplot(2,1,2);
plot(t_n, e_c_n, 'm', 'LineWidth', 1.2);
title('Σφάλμα Απόσβεσης e_c(t) με Θόρυβο (\eta_0 = 0.05)'); xlabel('Χρόνος (sec)'); grid on;

% =========================================================================
% 5. ΜΕΛΕΤΗ ΣΦΑΛΜΑΤΟΣ vs ΠΛΑΤΟΣ ΘΟΡΥΒΟΥ
% =========================================================================
eta0_vals = 0:0.01:0.1; % Δοκιμή για διάφορα πλάτη θορύβου
err_l_steady = zeros(length(eta0_vals), 1);
err_c_steady = zeros(length(eta0_vals), 1);

for i = 1:length(eta0_vals)
    [~, X_temp] = ode45(@(t,X) pendulum_lyapunov(t, X, theta1_star, theta2_star, gamma1, gamma2, am, am1, eta0_vals(i)), tspan, X0);
    
    l_temp = 9.81 ./ X_temp(:,5);
    c_temp = X_temp(:,6);
    
    % Υπολογισμός Mean Absolute Error στα τελευταία 20 sec (80% έως 100%)
    % (Μελετάμε τη μόνιμη κατάσταση, αγνοώντας το αρχικό μεταβατικό)
    idx_steady = round(0.8 * length(X_temp)) : length(X_temp);
    err_l_steady(i) = mean(abs(l_temp(idx_steady) - l_true));
    err_c_steady(i) = mean(abs(c_temp(idx_steady) - c_true));
end

% --- Figure 3: Διάγραμμα Σφάλματος σε σχέση με τον Θόρυβο ---
figure('Name', 'Θέμα 2.β: Σφάλμα vs Πλάτος Θορύβου', 'Position', [200, 200, 800, 400]);
plot(eta0_vals, err_l_steady, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(eta0_vals, err_c_steady, 'ms-', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
title('Μέσο Απόλυτο Σφάλμα Παραμέτρων (Μόνιμη Κατάσταση) vs Πλάτος Θορύβου \eta_0');
xlabel('Πλάτος Θορύβου \eta_0');
ylabel('Μέσο Απόλυτο Σφάλμα |e|');
legend('Σφάλμα Μήκους |e_l|', 'Σφάλμα Απόσβεσης |e_c|', 'Location', 'Best');
grid on;


% =========================================================================
% ΣΥΝΑΡΤΗΣΗ ΔΥΝΑΜΙΚΟΥ ΣΥΣΤΗΜΑΤΟΣ 
% =========================================================================
function dX = pendulum_lyapunov(t, X, th1_star, th2_star, gamma1, gamma2, am, am1, noise_amp)
    % Καταστάσεις
    x1 = X(1); x2 = X(2); x1_hat = X(3); x2_hat = X(4); th1_hat = X(5); th2_hat = X(6);
    
    u = 0.5 * sin(t); % Είσοδος (Συνθήκη Επιμένουσας Διέγερσης)
    
    % Θόρυβος στις μετρήσεις
    eta = noise_amp * sin(20 * pi * t);
    x1_m = x1 + eta; x2_m = x2 + eta; u_m  = u + eta;
    
    % Μη-Γραμμικές Συναρτήσεις (Series-Parallel τροφοδότηση με μετρήσεις)
    f1 = -sin(x1_m); 
    f2 = -x2_m;
    
    % ΠΡΑΓΜΑΤΙΚΟ ΣΥΣΤΗΜΑ
    dx1 = x2;
    dx2 = th1_star * f1 + th2_star * f2 + u; 
    
    % ΣΦΑΛΜΑ
    e = x2_m - x2_hat;
    
    % ΜΟΝΤΕΛΟ ΕΚΤΙΜΗΣΗΣ (Series-Parallel)
    dx1_hat = x2_hat + am1 * (x1_m - x1_hat); 
    dx2_hat = th1_hat * f1 + th2_hat * f2 + u_m + am * e;
    
    % ΝΟΜΟΙ ΠΡΟΣΑΡΜΟΓΗΣ LYAPUNOV
    dth1_hat = gamma1 * e * f1;
    dth2_hat = gamma2 * e * f2;
    
    dX = [dx1; dx2; dx1_hat; dx2_hat; dth1_hat; dth2_hat];
end