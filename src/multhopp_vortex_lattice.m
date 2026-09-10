%% PRÁCTICA ALAS - AERODINÁMICA

clear; clc; close all;

%%  Datos Grupo 8
tau1  = 0.95;
y1    = 7;
tau2  = 0.45;
B     = 20;
alpha0 = 0.2;
epsilon = @(y) 0.025*abs(y);
alpha_eff = @(y) alpha0 + epsilon(y);   % torsión SUMA
b = B/2;
c_raiz = 1;
a0 = 2*pi;

chord = @(y) arrayfun(@(yi) ...
    (abs(yi) <= y1) .* (c_raiz + (tau1*c_raiz - c_raiz)*(abs(yi)/y1)) + ...
    (abs(yi) >  y1) .* (tau1*c_raiz + (tau2*c_raiz - tau1*c_raiz)*((abs(yi)-y1)/(b-y1))), y);

%% Parte 1 – Método de Multhopp
fprintf('Parte 1 - Método de Multhopp\n');
m = 50;
i_vec = (1:m)';
theta  = i_vec * pi / (m+1);
y_col  = -b * cos(theta);
k_col  = chord(y_col) / (2*b);
alp_col = alpha_eff(y_col);

% Matriz C  [m x m] 
C = zeros(m, m);
for ii = 1:m
    for jj = 1:m
        C(ii,jj) = (1 + pi*k_col(ii)/(2*sin(theta(ii))) * jj) * sin(jj*theta(ii));
    end
end

% RHS
rhs = pi * k_col .* alp_col;    

A_coef = C \ rhs;

N_plot = 500;
theta_plot = linspace(pi/(N_plot+1), pi - pi/(N_plot+1), N_plot)';
y_plot     = -b * cos(theta_plot);

Gamma_mult = zeros(N_plot, 1);
for jj = 1:m
    Gamma_mult = Gamma_mult + A_coef(jj) * sin(jj * theta_plot);
end

fprintf('Coeficiente A1 (Multhopp) = %.6f\n', A_coef(1));

%% Parte 2 – Método Vortex Lattice
fprintf('Parte 2 - Método Vortex Lattice\n');

Ny = 40;
Nx = 10;

y_edges = -b + (b*2) * (1 - cos(linspace(0, pi, Ny+1))) / 2;
y_panel_center = 0.5*(y_edges(1:end-1) + y_edges(2:end));
dy_panel       = diff(y_edges);

N = Ny * Nx;

xj  = zeros(N,1);   yj  = zeros(N,1);
xjp = zeros(N,1);   yjp = zeros(N,1);
xi  = zeros(N,1);   yi  = zeros(N,1);

for iy = 1:Ny
    y_L = y_edges(iy);
    y_R = y_edges(iy+1);
    c_L = chord(y_L);
    c_R = chord(y_R);
    c_mid = chord(y_panel_center(iy));
    dx_L = c_L / Nx;
    dx_R = c_R / Nx;
    dx_mid = c_mid / Nx;
    for ix = 1:Nx
        p = (iy-1)*Nx + ix;
        xj(p)  = (ix-1)*dx_L + dx_L/4;    yj(p)  = y_L;
        xjp(p) = (ix-1)*dx_R + dx_R/4;    yjp(p) = y_R;
        xi(p)  = (ix-1)*dx_mid + 3*dx_mid/4;
        yi(p)  = y_panel_center(iy);
    end
end

fprintf('Calculando matriz de influencia (%dx%d)...\n', N, N);
Omega = zeros(N, N);
for ii = 1:N
    for jj = 1:N
        xi_  = xi(ii);   yi_  = yi(ii);
        xj_  = xj(jj);   yj_  = yj(jj);
        xjp_ = xjp(jj);  yjp_ = yjp(jj);

        % Hilo libre izquierdo
        if abs(yi_-yj_) > 1e-10
            w1 = -1/(4*pi*(yi_-yj_)) * (1 - (xj_-xi_)/sqrt((xi_-xj_)^2+(yi_-yj_)^2));
        else; w1 = 0; end

        % Hilo libre derecho
        if abs(yi_-yjp_) > 1e-10
            w3 = 1/(4*pi*(yi_-yjp_)) * ((xi_-xjp_)/sqrt((xi_-xjp_)^2+(yi_-yjp_)^2)+1);
        else; w3 = 0; end

        % Cabeza del torbellino
        r1x = xi_-xj_;    r1y = yi_-yj_;
        r2x = xi_-xjp_;   r2y = yi_-yjp_;
        rjx = xjp_-xj_;   rjy = yjp_-yj_;
        mod_r1 = sqrt(r1x^2+r1y^2);
        mod_r2 = sqrt(r2x^2+r2y^2);
        mod_rj = sqrt(rjx^2+rjy^2);
        if mod_r1>1e-10 && mod_r2>1e-10 && mod_rj>1e-10
            cos_t1 = (r1x*rjx+r1y*rjy)/(mod_r1*mod_rj);
            cos_t2 = (r2x*rjx+r2y*rjy)/(mod_r2*mod_rj);
            cross_z = rjx*r1y - rjy*r1x;
            sin_t1_k = cross_z/(mod_rj*mod_r1);
            h2 = mod_r1*abs(sin_t1_k);
            if abs(h2)>1e-10
                w2 = sign(sin_t1_k)/(4*pi*h2)*(cos_t1-cos_t2);
            else; w2=0; end
        else; w2=0; end

        Omega(ii,jj) = w1+w2+w3;
    end
end

fprintf('Matriz calculada. Resolviendo sistema...\n');

% RHS VLM
rhs_vl = zeros(N,1);
for ii = 1:N
    rhs_vl(ii) = (1/B) * (-alpha0 - epsilon(yi(ii)));
end

gamma_hat = Omega \ rhs_vl;

Gamma_vl_col = zeros(Ny, 1);
for iy = 1:Ny
    for ix = 1:Nx
        p = (iy-1)*Nx + ix;
        Gamma_vl_col(iy) = Gamma_vl_col(iy) + gamma_hat(p);
    end
end
y_vl = y_panel_center';

%% CL a lo largo de la envergadura (Vortex Lattice)
c_vl = chord(y_vl);
CL_vl = 2 * B * Gamma_vl_col ./ c_vl;

%% Gráfcas
figure(1);
set(gcf,'Position',[100 100 700 450]);
plot(y_plot, Gamma_mult, 'b-', 'LineWidth', 2);
xlabel('y [m]','FontSize',12);
ylabel('G(y)','FontSize',12);
title('Método de Multhopp – Circulación adimensional Grupo 8','FontSize',13);
grid on; xlim([-b b]);
legend(sprintf('\\tau_1=%.2f, \\tau_2=%.2f, B=%d, \\alpha=%.2f rad',tau1,tau2,B,alpha0),...
       'Location','north');

figure(2);
set(gcf,'Position',[200 100 700 450]);
plot(y_vl, Gamma_vl_col, 'r-', 'LineWidth', 2, 'MarkerSize', 4);
xlabel('y [m]','FontSize',12);
ylabel('G(y)','FontSize',12);
title('Método Vortex Lattice – Circulación adimensional Grupo 8','FontSize',13);
grid on; xlim([-b b]);
legend(sprintf('\\tau_1=%.2f, \\tau_2=%.2f, B=%d, \\alpha=%.2f rad',tau1,tau2,B,alpha0),...
       'Location','north');

figure(3);
set(gcf,'Position',[300 100 700 450]);
plot(y_plot, Gamma_mult, 'b-', 'LineWidth', 2); hold on;
plot(y_vl,   Gamma_vl_col, 'r-', 'LineWidth', 2);
xlabel('y [m]','FontSize',12);
ylabel('G(y)','FontSize',12);
title('Comparación Multhopp vs Vortex Lattice – Grupo 8','FontSize',13);
legend('Multhopp','VLattice','Location','north');
grid on; xlim([-b b]);

figure(4);
set(gcf,'Position',[400 100 700 450]);
plot(y_vl, CL_vl, 'b-', 'LineWidth', 2);
xlabel('y [m]','FontSize',12);
ylabel('C_L','FontSize',12);
title('Coeficiente de sustentación C_L a lo largo de la envergadura – Grupo 8','FontSize',13);
grid on; xlim([-b b]);

fprintf('\n Resultados resumen \n');
fprintf('Gamma_max Multhopp  (adim): %.6f\n', max(Gamma_mult));
fprintf('Gamma_max VLattice  (adim): %.6f\n', max(Gamma_vl_col));
fprintf('CL_max (VL):                %.6f\n', max(CL_vl));