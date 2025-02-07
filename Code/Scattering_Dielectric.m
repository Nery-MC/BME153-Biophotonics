
% Name : Nery Matias Calmo 
% BME 151 : Biophotonics 

% Title : Scattering Coefficent using Mie Theory Model
% Purpose : to calculate the scattering coefficent using the Mie model for 
% scattering by diaelectric spheres as a function of lamda [range: visible 
% spectrum (350 nm - 750 nm)] 

% Known Variables : 
% r(cell) = 8 um; 3000 organelles ea cells w/ avg r(s) = 300 nm; 
% n(s) = 1.4; n(m) = 1.36 

% Define Constants and Ranges 
lamda_range = linspace(350,700,100); % Wave Length Range (x-axis)
r_cell = 8e-6; % Radius of Cell (meters) 
num_org = 3000; % Number of Organelles per cell 
r_org = 300e-9; % Radius of ea Organelle (meters) 
ns = 1.4; % Refractive Index of Cell (scatter) 
nm = 1.36; % Refractive Index of the Medium 

% Initialize Array to Store Results 
scattering_coefficent = zeros(size(lamda_range)); 

% Calculate scattering Coefficient for ea Wavelength 

for i = 1:length(lamda_range)
    lamda = lamda_range(i); 

    % Use the function to calculate the scattering coefficent 
    scattering_coefficent(i) = calculate_mu(lamda, r_cell, r_org, num_org, ns, nm); 
end 

f = calculate_mu(630, r_cell, r_org, num_org, ns, nm);
disp('The scattering coefficent @630 nm : '); 
disp(f);

% Plot Results 
figure; 
plot(lamda_range, scattering_coefficent, 'b-', 'LineWidth', 2); 
xlabel('Wavelength (nm)');
ylabel('Scattering Coefficent (m^{-1})'); 
title('Scattering Coefficent vs Wavelength'); 
grid on; 

% Scattering Function 
function mu_s = calculate_mu(lamda, radius_cell, radius_organelle, num_organelles, n_s, n_m)
    % Calculate the Concentration of the Scatter (assuming sphere) 
    volume_cell = (4/3) * pi * (radius_cell)^3; % Volume of Cell 
    concentration_organelle = num_organelles / volume_cell; % Concentration of the Organelles 

    % Scattering Coefficent using Equation 
    a = 3.28 * pi * ((radius_organelle)^2) * concentration_organelle; 
    b = ((2 * pi * radius_organelle)/(lamda/(1e9)))^0.37; 
    c = ((n_s/n_m) - 1)^2.09; 
    x = a * b * c; 

    % Return calculated value for the scattering coefficent 
    mu_s = x;
end 
        