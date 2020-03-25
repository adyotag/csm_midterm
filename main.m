% Read in filetype
[node, element, elemType, nel, nen, nIntPts,nnd,ps, nu, E, Force_Node, bforce, disp_BC] = Read_input("Beam_Bending_Q9_16x4_PU.txt");

% Construct each element
elemList = [];
for i = 1:nel
    local_num = element(i,:);
    local_pos = reshape(node(local_num,:)',1,[]);
    elemList = [elemList, Element(i, local_num, local_pos, nen, nIntPts, ps, E, nu)];
    break
end

% Perform assembly of Global Stiffness Matrix


% element
% disp_BC

% nIntPts
% nnd % = number of nodes


% bforce
