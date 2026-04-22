clear; clc; close all;

% Παράμετροι 
J = 0.025;      % Ροπή αδράνειας 
k = 0.3;        % Πραγματική τιμή k (επιλογή στο [0.1, 0.5]) 
b = 1.0;        % Πραγματική τιμή b (επιλογή στο [0.5, 1.5]) 
gamma1 = 500;  % Κέρδος προσαρμογής για k
gamma2 = 500;  % Κέρδος προσαρμογής για b

tspan = [0 100]; % Χρόνος Προσομοίωσης
X0 = zeros(6,1); % Αρχικές συνθήκες: [phi, phi_dot, phi_hat, phi_dot_hat, k_hat, b_hat]
X0(5) = 0.1;
X0(6) = 0.5;
% Θέμα 1.α : Χωρίς Διαταραχή (d = 0)
d_a = @(t) 0; 
[t_a, X_a] = ode45(@(t,X) equations(t, X, J, k, b, gamma1, gamma2, d_a), tspan, X0);

% Εξαγωγή αποτελεσμάτων 1.α
phi_a = X_a(:,1); % Πραγματικό φ(t)
phi_hat_a = X_a(:,3); % Εκτίμηση φ_hat(t)
e_phi_a = phi_a - phi_hat_a; % Σφάλμα Μοντελοποίησης
e_k_a = X_a(:,5) - k; % Σφάλμα Παραμετροποίησης του k
e_b_a = X_a(:,6) - b; % Σφάλμα Παραμετροποίησης του b

phi_ddot_a = -k/J * X_a(:,2) + b/J * 0.25 * sin(0.5 * pi * t_a);
phi_ddot_hat_a = -X_a(:,5)/J .* X_a(:,4) + X_a(:,6)/J * 0.25 .* sin(0.5 * pi * t_a);

% Θέμα 1.β : Με Διαταραχή (d = 0.02*sin(2*t))
d_b = @(t) 0.02 * sin(2*t); 
[t_b, X_b] = ode45(@(t,X) equations(t, X, J, k, b, gamma1, gamma2, d_b), tspan, X0);

% Εξαγωγή αποτελεσμάτων 1.β
phi_b = X_b(:,1); phi_hat_b = X_b(:,3);
e_phi_b = phi_b - phi_hat_b;
e_k_b = X_b(:,5) - k;
e_b_b = X_b(:,6) - b;



% Figure 1: Ερώτημα 1.α (Χωρίς Διαταραχή)
figure('Name', 'Θέμα 1.α: Εκτίμηση χωρίς Διαταραχή');
subplot(4,1,1);
plot(t_a, phi_a, 'b', t_a, phi_hat_a, '--'); title('Γωνία \phi(t) vs \phi_{hat}(t)'); grid on;
subplot(4,1,2);
plot(t_a, e_phi_a, 'g'); title('Σφάλμα Γωνίας e_\phi(t)'); grid on;
subplot(4,1,3);
plot(t_a, e_k_a, 'k'); title('Σφάλμα Παραμέτρου e_k(t)'); grid on;
subplot(4,1,4);
plot(t_a, e_b_a, 'm'); title('Σφάλμα Παραμέτρου e_b(t)'); xlabel('Time (s)'); grid on;

% Figure 2: Ερώτημα 1.β (Με Διαταραχή)
figure('Name', 'Θέμα 1.β: Εκτίμηση με Διαταραχή');
subplot(4,1,1);
plot(t_b, phi_b, 'b', t_b, phi_hat_b, 'r--'); title('Γωνία \phi(t) vs \phi_{hat}(t) (με Διαταραχή)'); grid on;
subplot(4,1,2);
plot(t_b, e_phi_b, 'g'); title('Σφάλμα Γωνίας e_\phi(t)'); grid on;
subplot(4,1,3);
plot(t_b, e_k_b, 'k'); title('Σφάλμα Παραμέτρου e_k(t)'); grid on;
subplot(4,1,4);
plot(t_b, e_b_b, 'm'); title('Σφάλμα Παραμέτρου e_b(t)'); xlabel('Time (s)'); grid on;

% Figure 3: Σύγκριση Σφαλμάτων Παραμέτρων (Προαιρετικό για το σχολιασμό) 
figure('Name', 'Σύγκριση Σφαλμάτων k και b');
subplot(2,1,1);
plot(t_a, e_k_a, 'b', t_b, e_k_b, 'r'); title('Σφάλμα e_k: Χωρίς (μπλε) vs Με Διαταραχή (κόκκινο)'); grid on;
subplot(2,1,2);
plot(t_a, e_b_a, 'b', t_b, e_b_b, 'r'); title('Σφάλμα e_b: Χωρίς (μπλε) vs Με Διαταραχή (κόκκινο)'); grid on;


% --- Υπολογισμός Επιταχύνσεων (Για το Θέμα 1.α) ---
% Βοηθητικός υπολογισμός της εισόδου για τον χρόνο t_a
u_a = 0.25 * sin(0.5 * pi * t_a);

% Πραγματική επιτάχυνση (Φυσικό σύστημα)
% Χρησιμοποιούμε τα k, b και το πραγματικό phi_dot (X_a(:,2))
phi_ddot_a = (1/J) * (-k * X_a(:,2) + b * u_a);

% Εκτιμώμενη επιτάχυνση (Μοντέλο εκτίμησης)
% Χρησιμοποιούμε τα k_hat (X_a(:,5)), b_hat (X_a(:,6)) και το πραγματικό phi_dot (X_a(:,2))!
phi_ddot_hat_a = (1/J) * (-X_a(:,5) .* X_a(:,2) + X_a(:,6) .* u_a);

% Σφάλμα επιτάχυνσης (προαιρετικό αλλά πολύ χρήσιμο για το γράφημα)
e_phi_ddot_a = phi_ddot_a - phi_ddot_hat_a;

% --- Εκτύπωση Figure 20 ---
figure(20);
set(gcf, 'Name', 'Σύγκριση Επιταχύνσεων (1.α)'); % Ονομάζουμε το παράθυρο

% Πάνω γράφημα: Σύγκριση των δύο σημάτων
subplot(2,1,1);
plot(t_a, phi_ddot_a, 'b', 'LineWidth', 1.5); hold on;
plot(t_a, phi_ddot_hat_a, 'r--', 'LineWidth', 1.5);
title('Πραγματική Επιτάχυνση \phi_{ddot} vs Εκτιμώμενη \phi_{ddot\_hat}');
legend('\phi_{ddot}(t)', '\phi_{ddot\_hat}(t)');
grid on;

% Κάτω γράφημα: Το σφάλμα της επιτάχυνσης
subplot(2,1,2);
plot(t_a, e_phi_ddot_a, 'k', 'LineWidth', 1.5);
title('Σφάλμα Επιτάχυνσης (πρέπει να μηδενίζεται)');
xlabel('Χρόνος (sec)');
grid on;