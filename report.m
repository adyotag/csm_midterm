% Question 1
%% Part A
% Read Data
clear, close all, clc;

Assembly("Beam_Bending_Q4_16x4_Al.txt").run(); % 16x4 4-node bilinear mesh
ro = Read_output("Beam_Bending_Q4_16x4_Al.txt");
nsad = ro.nsad(); ips = ro.ips();
X = nsad(:,2); Y = nsad(:,3); XD = nsad(:,4); YD = nsad(:,5);
SXX = nsad(:, 6); SYY = nsad(:, 7); SXY = nsad(:, 8);
contour_num = 20;

% FEM solution
nPtsX = length(unique(X)); nPtsY = length(unique(Y));
Xr = reshape(X, [nPtsY, nPtsX]); XDr = reshape(XD, [nPtsY, nPtsX]);
Yr = reshape(Y, [nPtsY, nPtsX]); YDr = reshape(YD, [nPtsY, nPtsX]);
SXXr = reshape(SXX, [nPtsY, nPtsX]);
SYYr = reshape(SYY, [nPtsY, nPtsX]);
SXYr = reshape(SXY, [nPtsY, nPtsX]);

figure();
contourf(Xr + XDr, Yr + YDr, SXXr, contour_num);
title('$$\sigma_{xx}$$ Nodal Stress (MPa)', 'interpreter', 'latex'); colorbar(); axis equal; xlabel('X Position'); ylabel('Y Position');
figure();
contourf(Xr + XDr, Yr + YDr, SYYr, contour_num);
title('$$\sigma_{yy}$$ Nodal Stress (MPa)', 'interpreter', 'latex'); colorbar(); axis equal; xlabel('X Position'); ylabel('Y Position');
figure();
contourf(Xr + XDr, Yr + YDr, SXYr, contour_num);
title('$$\sigma_{xy}$$ Nodal Stress (MPa)', 'interpreter', 'latex'); colorbar(); axis equal; xlabel('X Position'); ylabel('Y Position');

%% Part B
clear, close all, clc;

% EB Solution
q = 2.71E-9 * 9.82E3;
I = 1./12.; E = 70E3; L = 12;
alpha = q/(E*I);
deflection = @(x) (1/6912 - 5*alpha/4.)*(x.^3) + (9*alpha - 1/192)*(x.^2) + (alpha/24.)*(x.^4);
moment = @(x) E*I*( (1/1152 - 7.5*alpha)*x + (18*alpha - 1/96.) + 0.5*alpha*(x.^2));
bending_stress_an = @(x, yn) -(moment(x).*yn)/I; bending_stress_LHS = @(yn) bending_stress_an(0, yn);
shear_stress_an = @(x) -E*I*(alpha*x + 1/1152 - 7.5*alpha);
reaction_force_an = [integral(bending_stress_LHS, -0.5, 0.5), -shear_stress_an(0)];

% 4-point integration
% Different element sizes
A441 = Assembly("Beam_Bending_Q4_4x1_Al.txt").run(); % 4x1, 4-node
ro = Read_output("Beam_Bending_Q4_4x1_Al.txt"); nsad = ro.nsad(); ips = ro.ips();
X441 = nsad(:,2); Y441 = nsad(:,3); XD441 = nsad(:,4); YD441 = nsad(:,5);
SXX441 = nsad(:, 6); SYY441 = nsad(:, 7); SXY441 = nsad(:, 8);

A482 = Assembly("Beam_Bending_Q4_8x2_Al.txt").run(); % 8x2, 4-node
ro = Read_output("Beam_Bending_Q4_8x2_Al.txt"); nsad = ro.nsad(); ips = ro.ips();
X482 = nsad(:,2); Y482 = nsad(:,3); XD441 = nsad(:,4); YD482 = nsad(:,5);
SXX482 = nsad(:, 6); SYY482 = nsad(:, 7); SXY482 = nsad(:, 8);

A4164 = Assembly("Beam_Bending_Q4_16x4_Al.txt").run(); % 16x4, 4-node
ro = Read_output("Beam_Bending_Q4_16x4_Al.txt"); nsad = ro.nsad(); ips = ro.ips();
X4164 = nsad(:,2); Y4164 = nsad(:,3); XD4164 = nsad(:,4); YD4164 = nsad(:,5);
SXX4164 = nsad(:, 6); SYY4164 = nsad(:, 7); SXY4164 = nsad(:, 8);

A4168 = Assembly("Beam_Bending_Q4_16x8_PU.txt").run(); % 16x8, 4-node
ro = Read_output("Beam_Bending_Q4_16x8_PU.txt"); nsad = ro.nsad(); ips = ro.ips();
X4168 = nsad(:,2); Y4168 = nsad(:,3); XD4168 = nsad(:,4); YD4168 = nsad(:,5);
SXX4168 = nsad(:, 6); SYY4168 = nsad(:, 7); SXY4168 = nsad(:, 8);
A4168.trf

figure();
scatter(X4168+XD4168, Y4168+YD4168);


%% Different nodes

% [ *** "Beam_Bending_Q4_4x1_Al.txt" already processed ! *** ]

A841 = Assembly("Beam_Bending_Q8_4x1_Al.txt").run(); % 4x1, 4-node
ro = Read_output("Beam_Bending_Q4_4x1_Al.txt");
nsad = ro.nsad(); ips = ro.ips();
X841 = nsad(:,2); Y841 = nsad(:,3); XD841 = nsad(:,4); YD841 = nsad(:,5);
SXX841 = nsad(:, 6); SYY841 = nsad(:, 7); SXY841 = nsad(:, 8);

A941 = Assembly("Beam_Bending_Q9_4x1_Al.txt").run(); % 4x1, 9-node
ro = Read_output("Beam_Bending_Q4_4x1_Al.txt"); nsad = ro.nsad(); ips = ro.ips();
X941 = nsad(:,2); Y941 = nsad(:,3); XD941 = nsad(:,4); YD941 = nsad(:,5);
SXX941 = nsad(:, 6); SYY941 = nsad(:, 7); SXY941 = nsad(:, 8);

% Compute errors in reaction force
rar = [A441,A482,A4164,A4168, A841, A941]';
rar.trf


