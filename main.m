%% =========================================================================
%  Code developed by Ali RAZA in scope of Torrent et al. (Phys. Rev. B 97, 014105, 2018)
%  Loss compensation in time-dependent elastic metamaterials
% =========================================================================

clear; close all; clc;

%% 1. SYSTEM PARAMETERS
% Defined in the text description of Figure 3.
% Material A
rhoA = 1.0; CA = 1.0; etaA = 0.1;
% Material B
rhoB = 1.0; CB = 3.0; etaB = 0.2;

% Derived quantities: Speed c, Impedance z (Eq. A6)
cA = sqrt(CA / rhoA); cB = sqrt(CB / rhoB);
zA = sqrt(CA * rhoA); zB = sqrt(CB * rhoB);

% Time modulation period
tau = 1.0;

%% 2. FIGURE 3: DISPERSION AND STABILITY


k0_c0 = linspace(0.001, 6.2, 1000);
tau_A_cases = [0.25, 0.5, 0.8];

figure('Name', 'Figure 3', 'Color', 'w', 'Position', [50 100 1200 600]);

k0_peak_gain = 2.2; % Initial guess

for i = 1:length(tau_A_cases)
    tauA_frac = tau_A_cases(i);
    tA = tauA_frac * tau;
    tB = (1.0 - tauA_frac) * tau;
    
    % Average viscosity parameter <eta/rho> from Eq. (12)
    eta_rho_avg = (tA*(etaA/rhoA) + tB*(etaB/rhoB))/tau;
    
    phiA = k0_c0 * cA * tA;
    phiB = k0_c0 * cB * tB;
    
    % Dispersion relation for two-layer system: Trace of Transfer Matrix
    % See Eq. (A27) and (A26): cos(Omega*T) = 0.5 * tr(P(T))
    cos_Om = cos(phiA).*cos(phiB) - 0.5*(zA/zB + zB/zA)*sin(phiA).*sin(phiB);
    
    % Bloch-Floquet frequency Omega (Eq. A25)
    Om_tau = acos(cos_Om); 
    
    % Stability Parameter Delta(k0) defined in Eq. (14)
    % Delta = 1/2 <eta/rho> k0^2 - Im(Omega)
    Delta = 0.5 * eta_rho_avg * (k0_c0.^2) - (abs(imag(Om_tau)) / tau);
    
    if tauA_frac == 0.5
        [~, idx_min] = min(Delta);
        k0_peak_gain = k0_c0(idx_min);
    end
    
    subplot(2, 3, i);
    plot(k0_c0, real(Om_tau), 'LineWidth', 2.5);
    title(['\tau_A = ' num2str(tauA_frac) '\tau']);
    xlabel('k_0 c_0 \tau'); ylabel('Re(\Omega) \tau');
    axis([0 6.2 0 3.2]); grid on;
    
    subplot(2, 3, i+3);
    plot(k0_c0, Delta, 'r', 'LineWidth', 2);
    hold on; yline(0, 'k--');
    xlabel('k_0 c_0 \tau'); ylabel('\Delta(k_0)');
    axis([0 6.2 -0.5 3.5]); grid on;
end

%% 3. FIGURE 4: ENERGY GAIN 

k_E = linspace(0.001, 3.1, 400);
N_list = [2, 4, 6, 8, 10];

% Interface Matrices based on Impedance Mismatch (Appendix A4, Eq. A22 structure)
rAB = zA/zB; T_AB = 0.5*[1+rAB, 1-rAB; 1-rAB, 1+rAB];
rBA = zB/zA; T_BA = 0.5*[1+rBA, 1-rBA; 1-rBA, 1+rBA];

figure('Name', 'Figure 4', 'Color', 'w', 'Position', [100 150 600 800]);

for i = 1:length(tau_A_cases)
    tauA_frac = tau_A_cases(i);
    tA = tauA_frac * tau;
    tB = (1.0 - tauA_frac) * tau;
    eta_rho_avg = (tA*(etaA/rhoA) + tB*(etaB/rhoB))/tau;
    
    subplot(3, 1, i); hold on;
    
    for j = 1:length(N_list)
        N = N_list(j);
        E_gain = zeros(size(k_E));
        
        for k_idx = 1:length(k_E)
            k = k_E(k_idx);
            
            % Propagation Phase Matrices (Eq. A21, diagonal form)
            phiA = k * cA * tA;
            phiB = k * cB * tB;
            PA = diag([exp(-1i*phiA), exp(1i*phiA)]);
            PB = diag([exp(-1i*phiB), exp(1i*phiB)]);
            
            % Unit Cell Matrix M_cell (product of interface and propagation matrices)
            M_cell = T_BA * PB * T_AB * PA;
            
            % Total Matrix for N periods (Eq. 7 implies M_N = M^N)
            M_total = M_cell^N;
            
            % Damping factor Gamma from Eq. (10) and (12)
            % Damping term = exp(-Gamma * k0^2) per Eq. (9)
            kappa = N * tau * 0.5 * eta_rho_avg;
            damp = exp(-kappa * k^2);
            
            % Reflection (rT) and Transmission (tT) Coefficients
            % Defined in Eq. (6a) and (6b): rT = M_21, tT = M_11
            tT = M_total(1,1) * damp;
            rT = M_total(2,1) * damp;
            
            % Total Energy Gain E_T/E_0 = |rT|^2 + |tT|^2 (Caption Fig. 4)
            E_gain(k_idx) = abs(tT)^2 + abs(rT)^2;
        end
        plot(k_E, E_gain, 'LineWidth', 1.5);
    end
    title(['\tau_A = ' num2str(tauA_frac) '\tau']);
    ylabel('E_T / E_0'); grid on; xlim([0 3.1]);
    if i == 1; legend('N=2','N=4','N=6','N=8','N=10'); end
    if i == 3; xlabel('k_0 c_0 \tau'); end
end


%% 4. FIGURE 5: PULSE EVOLUTION

% we Simulation Parameters 
L = 200;            % Domain length
Nx = 4096;          % Points for FFT
x = linspace(-L/2, L/2, Nx);
dx = x(2)-x(1);
k_axis = (2*pi/L) * [0:Nx/2-1, -Nx/2:-1]; % Standard FFT k-ordering

% we Pulse Parameters (Matches Fig. 5 text) 
N_pulse = 25;       % Modulation lasts 25 periods
A0 = 2.5;           % Initial Amplitude
sigma_x = 2.0;      % Pulse width
k_carrier = k0_peak_gain; 

% we Initial Pulse u(x,0) 
u0_x = A0 * exp(-x.^2/(2*sigma_x^2)) .* cos(k_carrier * x);
U0_k = fft(u0_x);   

% we Gain Parameters (Matches Eq. 12) 
tA = 0.5 * tau; tB = 0.5 * tau;
eta_rho_avg = (tA*(etaA/rhoA) + tB*(etaB/rhoB))/tau;

% we Pre-computation: Eigen-Decomposition for Speed

threshold = 1e-6 * max(abs(U0_k));
active_indices = find(abs(U0_k) > threshold);
num_active = length(active_indices);

% Storage for decomposition
V_store = zeros(2, 2, num_active);       % Eigenvectors
V_inv_store = zeros(2, 2, num_active);   % Inverse Eigenvectors
D_diag_store = zeros(2, num_active);     % Eigenvalues

for j = 1:num_active
    idx = active_indices(j);
    k = k_axis(idx);
    
    % Modulation Phase (Eq. 7 / Eq. 11 logic)
    phiA = k * cA * tA; phiB = k * cB * tB;
    PA = diag([exp(-1i*phiA), exp(1i*phiA)]);
    PB = diag([exp(-1i*phiB), exp(1i*phiB)]);
    
    % Single Unit Cell Matrix
    M_cell = T_BA * PB * T_AB * PA;
    
    % Eigen-decomposition
    [V, D] = eig(M_cell);
    
    V_store(:,:,j) = V;
    V_inv_store(:,:,j) = inv(V);
    D_diag_store(:,j) = diag(D);
end

%  Animation Setup 
figure('Name', 'Figure 5 Animation', 'Color', 'w', 'Position', [150 200 1000 400]);
hLine = plot(x, u0_x, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
grid on; 
xlabel('x / (c_0 \tau)'); ylabel('u(x)/u_0');
xlim([-60 60]); 
% Set initial Y-limits based on A0
ylim([-A0*1.5, A0*1.5]); 
titleHandle = title('t/\tau = 0');

% Animation Loop 
% Matches timeline: 0 -> 25 (Growth) -> 60 (Drift/Split)
times_anim = 0:0.5:50; 

for t_val = times_anim
    
    % 1. Determine Physics State (Inside vs Outside Time-Crystal)
    if t_val <= N_pulse * tau
        % Phase 1: Inside Crystal (Eq. 7 / Eq. 13)
        % The pulse stays localized (group velocity ~ 0 in gap) and grows.
        current_N = t_val / tau;
        dt_prop = 0; 
    else
        % Phase 2: Exited Crystal
        % Modulation stops. The accumulated field splits and drifts.
        current_N = N_pulse;
        dt_prop = t_val - N_pulse * tau;
    end
    
    % 2. Fourier Synthesis
    U_k_out = zeros(size(U0_k));
    
    % Viscous damping factor (Eq. 9): exp(-Gamma * k^2)
    % Kappa scales with the duration of interaction (current_N)
    kappa_eff = current_N * tau * 0.5 * eta_rho_avg;
    
    % Pre-calculate drift phase (only active after t > 25)
    if dt_prop > 0
        % Forward drift phase
        phase_drift_vec = exp(-1i * k_axis * cA * dt_prop);
    else
        phase_drift_vec = ones(size(k_axis));
    end

    for j = 1:num_active
        idx = active_indices(j);
        k = k_axis(idx);
        
        % Compute M^N using eigenvalues (Fast)
        D_pow = diag(D_diag_store(:,j).^current_N);
        M_total = V_store(:,:,j) * D_pow * V_inv_store(:,:,j);
        
        % Apply Damping (Eq. 9)
        damp = exp(-kappa_eff * k^2);
        
        % Extract Coefficients (Eq. 6a, 6b)
        % tT = M11 (Forward amplitude)
        % rT = M21 (Backward amplitude)
        tT = M_total(1,1) * damp;
        rT = M_total(2,1) * damp;
        
        % Apply Drift Phase
        p_prop = phase_drift_vec(idx);
        
        % Forward Wave Construction
        U_k_out(idx) = U_k_out(idx) + U0_k(idx) * tT * p_prop;
        
        % Backward Wave Construction
        %  "reflected" wave travels backwards in x.
        % We map it to the negative k bin. It inherits the source phase
        % which naturally creates the left-propagating behavior (x+ct).
        if idx == 1; idx_neg = 1; else; idx_neg = Nx - idx + 2; end
        U_k_out(idx_neg) = U_k_out(idx_neg) + U0_k(idx) * rT * p_prop;
    end
    
    % 3. Update Visuals
    u_real = real(ifft(U_k_out));
    
    set(hLine, 'YData', u_real);
    set(titleHandle, 'String', ['t/\tau = ' num2str(t_val, '%.1f')]);
    
    % Auto-scale Y-axis (Crucial as pulse grows significantly)
    current_max = max(abs(u_real));
    if current_max > 1e-6
        ylim([-current_max*1.2, current_max*1.2]);
    end
    
    drawnow limitrate;
end

