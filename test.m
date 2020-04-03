%% Benchmark for time

a = Assembly("Biaxial_Q4_2x2.txt").run();
b = Assembly("Beam_Bending_Q4_4x1_Al.txt").run();
c = Assembly("Beam_Bending_Q8_4x1_Al.txt").run();
d = Assembly("Beam_Bending_Q9_4x1_Al.txt").run();
e = Assembly("Beam_Bending_Q4_8x2_Al.txt").run();
f = Assembly("Beam_Bending_Q4_16x4_Al.txt").run();
g = Assembly("Beam_Bending_Q4_16x4_PU.txt").run();
h = Assembly("Beam_Bending_Q9_16x4_PU.txt").run();


%% Sample
clear; clean; close all; clc;
% Since this is structured using OOP, one can use the workspace to
% explore the different attributes easily

g = Assembly("Beam_Bending_Q4_16x4_PU.txt");    % load text file

% At this stage, we can still make changes to certain parameters. For example,
g.nIntPts;
g.nIntPts = 4;
g.nIntPts;

g.run();    % Run calculations

g.elemList;    % list of elements in the Assembly g

% Each assembly has useful attributes. These attributes are listed under
% "property" in the variable explorer/workspace. For example,

g.k_global_final;    % Global Stiffness Matrix
g.f_global_final;   % global loading vector
g.d_global_final;   % global displacements
g.m_global_final;   % Global mass matrix
g.p_global_final;   % Global projection vector
g.n_stress_final;   % Nodal stress vector
g.rf_global_final;  % Reaction forces experienced at specified nodes
g.trf;              % total reaction force

% If needed, you can browse even further and look at each individul element. 
% Suppose we wanted to look at element 1, then

g.elemList(1);

% One can use the workspace to explore each element's attributes accordingly
% Finally, to plot stresses, one can load that information as follows.

% Alternatively, to read the integration point stresses file,
info = g.readout().ips();

% To load the nodal stresses and displacements file,
info = g.readout().nsad();

% One can refer to the relevant text files, to understand which columns refer
% to what. One can use this information to graph displacements for example.

X = info(:,2); Y = info(:,3);   % Read in X and Y
XD = info(:,4); YD = info(:,5); % Read in X and Y Displacements
SXX = info(:, 6);
scatter(X+XD, Y+YD, 10, SXX, 'filled')       % plot of deformed configuration with bending stress
colorbar()









%