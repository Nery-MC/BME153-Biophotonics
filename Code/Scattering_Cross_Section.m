% Name : Nery Matias Calmo 
% BME 151 : Biophotonics 

% Title : Van de Haulst Approximation using the Cross Section Function 
% Purpose : to calculate the scattering cross section as a function of 
% wavelength for both normal and pre-cancerous cells. Assuming spherical 
% scattering

% Known Variables : 
% n(nucleus) = 1.42; n(cytoplasm) = 1.36; Wavelength Range = 300 - 800 nm;
% d(normal) = 8.5 um; d(pre-cancerous) = 10.2 um 

% Define Constants and Ranges 
lamda_range = linspace(300,800,100); % Wave Length Range (x-axis)
d_normal = 8.5e-6; % Diameter of Normal Nuclei (meters) 
d_pcancer = 10.2e-6; % Diameter of Pre-Cancerous Nuclei (meters) 
nn = 1.42; % Refractive Index of Nuclei 
nc = 1.36; % Refractive Index of Cytoplasm 

% Initialize Array to Store Results 
Cross_Section_Normal = zeros(size(lamda_range)); 
Cross_Section_PreCancerous = zeros(size(lamda_range)); 


% Calculate Scattering Cross Section for ea Wavelength 

for i = 1:length(lamda_range)
    lamda = lamda_range(i); 

    % Use the function to calculate the scattering coefficent for normal
    % cells 
    Cross_Section_Normal(i) = calculate_sigma(lamda, d_normal, nc, nn); 
    Cross_Section_PreCancerous(i) = calculate_sigma(lamda, d_pcancer, nc, nn); 

end 

% Plot Results 
figure 
plot(lamda_range, Cross_Section_Normal, 'b-', 'LineWidth', 2); 

hold on 
plot(lamda_range, Cross_Section_PreCancerous, 'r-', 'LineWidth', 2); 
xlabel('Wavelength (nm)');
ylabel('Scattering Cross-Section (um^{2})'); 
title('Scattering Cross-Section vs Wavelength'); 
legend('Normal', 'Pre-Cancerous'); 

grid on; 

% Scattering Cross-Section Function 
function sigma_s = calculate_sigma(lamda, d_particle, n1, n2) 

    % Convert Lamda nm to m
    lamda = lamda * 10e-9; 
    % Calculate the Wave Vector 
    k = (2 * pi) / (lamda); 
    % Calculate Relative Refractive Index 
    m = n2 / n1; 
    % Calculate the Phase Lag 
    g = k * d_particle * abs(m-1); 

    % Scattering Cross-Section using Equation 
    a = (pi * d_particle^2) / 4; 
    b = (4 * sin(g)) / g; 
    c = (4 / g) * (1 - cos(g)); 
    x_m = a * (2 - b + c); % Cross-section in Meters
    x = x_m * 10e6; % Converting to um 

    % Return calculated value for the scattering coefficent 
    sigma_s = x;

end 
