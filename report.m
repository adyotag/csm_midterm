%% Question 1 %%
%% Part A
% Read Data
clear, close all, clc; clean;

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
clear, close all, clc; clean;

% First run using 4-node integration

% EB Solution
q = 2.71E-9 * 9.82E3;
I = 1./12.; E = 70E3; L = 12;
alpha = q/(E*I);
deflection = @(x) (1/6912 - 5*alpha/4.)*(x.^3) + (9*alpha - 1/192)*(x.^2) + (alpha/24.)*(x.^4);
moment = @(x) E*I*( (1/1152 - 7.5*alpha)*x + (18*alpha - 1/96.) + 0.5*alpha*(x.^2));
bending_stress_an = @(x, yn) -(moment(x).*yn)/I; bending_stress_LHS = @(yn) bending_stress_an(0, yn);
shear_stress_an = @(x) -E*I*(alpha*x + 1/1152 - 7.5*alpha);
reaction_force_an = [integral(bending_stress_LHS, -0.5, 0.5), -shear_stress_an(0)];

% First load different assemblies
A441 = Assembly("Beam_Bending_Q4_4x1_Al.txt"); % 4x1, 4-node
A482 = Assembly("Beam_Bending_Q4_8x2_Al.txt"); % 8x2, 4-node
A4164 = Assembly("Beam_Bending_Q4_16x4_Al.txt"); % 16x4, 4-node
A841 = Assembly("Beam_Bending_Q8_4x1_Al.txt"); % 4x1, 4-node
A941 = Assembly("Beam_Bending_Q9_4x1_Al.txt"); % 4x1, 9-node

assem_list = [A441,A482,A4164, A841, A941]';
[assem_list.nIntPts] = deal(4);
for iter = 1:length(assem_list)
   assem_list(iter).run();
end

% !!!   The warning for the badly conditioned matrix isn't from my code:            !!!
% !!!   try running the benchmark, which goes through each input file provided.     !!!

% Different element sizes
ro = A441.readout(); nsad = ro.nsad();
X441 = nsad(:,2); Y441 = nsad(:,3); XD441 = nsad(:,4); YD441 = nsad(:,5);
SXX441 = nsad(:, 6); SYY441 = nsad(:, 7); SXY441 = nsad(:, 8);

ro = A482.readout(); nsad = ro.nsad();
X482 = nsad(:,2); Y482 = nsad(:,3); XD482 = nsad(:,4); YD482 = nsad(:,5);
SXX482 = nsad(:, 6); SYY482 = nsad(:, 7); SXY482 = nsad(:, 8);

ro = A4164.readout(); nsad = ro.nsad();
X4164 = nsad(:,2); Y4164 = nsad(:,3); XD4164 = nsad(:,4); YD4164 = nsad(:,5);
SXX4164 = nsad(:, 6); SYY4164 = nsad(:, 7); SXY4164 = nsad(:, 8);

% Different nodes

% [ *** "Beam_Bending_Q4_4x1_Al.txt" already processed ! *** ]

ro = A841.readout(); nsad = ro.nsad();
X841 = nsad(:,2); Y841 = nsad(:,3); XD841 = nsad(:,4); YD841 = nsad(:,5);
SXX841 = nsad(:, 6); SYY841 = nsad(:, 7); SXY841 = nsad(:, 8);

ro = A941.readout(); nsad = ro.nsad();
X941 = nsad(:,2); Y941 = nsad(:,3); XD941 = nsad(:,4); YD941 = nsad(:,5);
SXX941 = nsad(:, 6); SYY941 = nsad(:, 7); SXY941 = nsad(:, 8);

% Compute errors in reaction force
numer = abs([assem_list.trf] - repmat(reaction_force_an, 1, 5));
denom = repmat(reaction_force_an, 1, 5);
log_error = log(numer(2:2:end)./denom(2:2:end));
log_char_size = log([3., 1.5, 0.75, 3., 3.]);

figure();
scatter(log_char_size(1), log_error(1)); hold on
scatter(log_char_size(2), log_error(2));
scatter(log_char_size(3), log_error(3));
plot(log_char_size(1:3), log_error(1:3));
legend('4x1', '8x2', '16x4')
title('$$\log(e)$$ vs. $$\log(l)$$ for 4x1, 8x2, and 16x4 Assemblies with 4-Point Integration for Al Q4 Element', 'interpreter', 'latex');
xlabel('l, characteristic length'); ylabel('e, Error in y-reaction force');

figure();
scatter(log_char_size(1), log_error(1)); hold on
scatter(log_char_size(4), log_error(4));
scatter(log_char_size(5), log_error(5));
legend('Q4', 'Q8', 'Q9');
title('$$\log(e)$$ vs. $$\log(l)$$ for Q4, Q8, and Q9 Elements with 4-Point Integration for Al 4x1 Assembly', 'interpreter', 'latex');
xlabel('l, characteristic length'); ylabel('e, Error in y-reaction force');


% Now run using 9-node integration
clean; clear; clc;

% EB Solution
q = 2.71E-9 * 9.82E3;
I = 1./12.; E = 70E3; L = 12;
alpha = q/(E*I);
deflection = @(x) (1/6912 - 5*alpha/4.)*(x.^3) + (9*alpha - 1/192)*(x.^2) + (alpha/24.)*(x.^4);
moment = @(x) E*I*( (1/1152 - 7.5*alpha)*x + (18*alpha - 1/96.) + 0.5*alpha*(x.^2));
bending_stress_an = @(x, yn) -(moment(x).*yn)/I; bending_stress_LHS = @(yn) bending_stress_an(0, yn);
shear_stress_an = @(x) -E*I*(alpha*x + 1/1152 - 7.5*alpha);
reaction_force_an = [integral(bending_stress_LHS, -0.5, 0.5), -shear_stress_an(0)];


% First load different assemblies
A441 = Assembly("Beam_Bending_Q4_4x1_Al.txt"); % 4x1, 4-node
A482 = Assembly("Beam_Bending_Q4_8x2_Al.txt"); % 8x2, 4-node
A4164 = Assembly("Beam_Bending_Q4_16x4_Al.txt"); % 16x4, 4-node
A841 = Assembly("Beam_Bending_Q8_4x1_Al.txt"); % 4x1, 4-node
A941 = Assembly("Beam_Bending_Q9_4x1_Al.txt"); % 4x1, 9-node

assem_list = [A441,A482,A4164, A841, A941]';
[assem_list.nIntPts] = deal(9);
for iter = 1:length(assem_list)
   assem_list(iter).run();
end

% Different element sizes
ro = A441.readout(); nsad = ro.nsad();
X441 = nsad(:,2); Y441 = nsad(:,3); XD441 = nsad(:,4); YD441 = nsad(:,5);
SXX441 = nsad(:, 6); SYY441 = nsad(:, 7); SXY441 = nsad(:, 8);

ro = A482.readout(); nsad = ro.nsad();
X482 = nsad(:,2); Y482 = nsad(:,3); XD482 = nsad(:,4); YD482 = nsad(:,5);
SXX482 = nsad(:, 6); SYY482 = nsad(:, 7); SXY482 = nsad(:, 8);

ro = A4164.readout(); nsad = ro.nsad();
X4164 = nsad(:,2); Y4164 = nsad(:,3); XD4164 = nsad(:,4); YD4164 = nsad(:,5);
SXX4164 = nsad(:, 6); SYY4164 = nsad(:, 7); SXY4164 = nsad(:, 8);

% Different nodes

% [ *** "Beam_Bending_Q4_4x1_Al.txt" already processed ! *** ]

ro = A841.readout(); nsad = ro.nsad();
X841 = nsad(:,2); Y841 = nsad(:,3); XD841 = nsad(:,4); YD841 = nsad(:,5);
SXX841 = nsad(:, 6); SYY841 = nsad(:, 7); SXY841 = nsad(:, 8);

ro = A941.readout(); nsad = ro.nsad();
X941 = nsad(:,2); Y941 = nsad(:,3); XD941 = nsad(:,4); YD941 = nsad(:,5);
SXX941 = nsad(:, 6); SYY941 = nsad(:, 7); SXY941 = nsad(:, 8);

% Compute errors in reaction force
numer = abs([assem_list.trf] - repmat(reaction_force_an, 1, 5));
denom = repmat(reaction_force_an, 1, 5);
log_error = log(numer(2:2:end)./denom(2:2:end));
log_char_size = log([3., 1.5, 0.75, 3., 3.]);

figure();
scatter(log_char_size(1), log_error(1)); hold on
scatter(log_char_size(2), log_error(2));
scatter(log_char_size(3), log_error(3));
plot(log_char_size(1:3), log_error(1:3));
legend('4x1', '8x2', '16x4')
title('$$\log(e)$$ vs. $$\log(l)$$ for 4x1, 8x2, and 16x4 Assemblies with 9-Point Integration for Al Q4 Element', 'interpreter', 'latex');
xlabel('l, characteristic length'); ylabel('e, Error in y-reaction force');

figure();
scatter(log_char_size(1), log_error(1)); hold on
scatter(log_char_size(4), log_error(4));
scatter(log_char_size(5), log_error(5));
legend('Q4', 'Q8', 'Q9');
title('$$\log(e)$$ vs. $$\log(l)$$ for Q4, Q8, and Q9 Elements with 9-Point Integration for Al 4x1 Assembly', 'interpreter', 'latex');
xlabel('l, characteristic length'); ylabel('e, Error in y-reaction force');


%% Part C


























%%