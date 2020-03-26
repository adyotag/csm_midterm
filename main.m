% Read in filetype
[node, element, elemType, nel, nen, nIntPts,nnd,ps, nu, E, Force_Node, bforce, disp_BC] = ...
                Read_input("Beam_Bending_Q9_16x4_PU.txt");

% Construct each element
elemList = [];
for i = 1:nel
    local_num = element(i,:);
    local_pos = reshape(node(local_num,:)',1,[]);
    elemList = [elemList, Element(i, local_num, local_pos, nen, nIntPts, ps, E, nu, bforce)];
end

% Perform assembly of Global Stiffness Matrix
k_global = zeros(nnd*2);
f_global = zeros(nnd*2,1); 
for i = 1:nel
    k_global = k_global + elemList(i).getGlobalStiffness(nnd*2);
    f_global = f_global + elemList(i).getGlobalLoading(nnd*2);    
end

% Apply Displacement BC
globaleqns = 2*(disp_BC(:,1) - 1) + disp_BC(:,2);
dispamt = disp_BC(:,3);

for i = 1:length(dispamt)
    idx = globaleqns(i);
    f_global = f_global - k_global(:,idx)*dispamt(i);
    k_global(:,idx) = zeros(nnd*2,1); k_global(idx,:) = zeros(nnd*2,1);
    k_global(idx,idx) = 1.;
end

% Reduce condition of matrix
mask = ~logical(sum(globaleqns == (1:nnd*2),1));
k_global = k_global(:, mask); k_global = k_global(mask, :);
f_global = f_global(mask);

% Reassemble final solution
d_global = zeros(2*nnd,1);
d_global(mask) = inv(k_global)*f_global; d_global(~mask) = dispamt;




