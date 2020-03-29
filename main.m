% Read in filetype
[node, element, elemType, nel, nen, nIntPts,nnd,ps, nu, E, Force_Node, bforce, disp_BC] = ...
                Read_input("Biaxial_Q4_2x2.txt");

% Construct each element. Preallocating space for objects seems like a rabbit hole ;-)
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
d_global(mask) = k_global\f_global;
d_global(globaleqns) = dispamt;

% Update local displacements for each element
for i = 1:nel
   elemList(i).updateLocalDisplacements(d_global); 
end

% Compute Nodal Stresses
m_global = zeros(3*nnd);
p_global = zeros(3*nnd,1);
for i = 1:nel
    m_global = m_global + elemList(i).getGlobalMassMatrix(3*nnd);
    p_global = p_global + elemList(i).getGlobalProjectionVector(3*nnd);
end
n_stress = m_global\p_global;

% Compute Reaction Forces
rf_global = zeros(2*nnd,1);
num_rf = 2*length(Force_Node);
for i = 1:nel
    rf_global = rf_global + elemList(i).getGlobalRFVector(2*nnd);
end

rfgidx = 2*(reshape(repmat(Force_Node',2,1),[num_rf,1]) - 1) + 1 + cast(mod((1:num_rf) - 1, 2)', 'uint32');
rf_global = rf_global(rfgidx);

% Compute Total Reaction Force
trf = sum(reshape(rf_global, [num_rf/2, 2]),1);

% Create output file for nodal information
farr = zeros(nnd, 8);
farr(:,1) = 1:nnd; farr(:,2:3) = node;
farr(:, 4:5) = reshape(d_global',[2,nnd])';
farr(:, 6:end) = reshape(n_stress',[3,nnd])';

fileID = fopen('NodalStressAndDisp.txt','w');
fprintf(fileID, '|\tNode_Num\t|\tX\t|\tY\t|\tX_Disp\t|\tY_Disp\t|\tNode_Str_XX\t|\tNode_Str_YY\t|\tNode_Str_XY\t|\n\n');
fprintf(fileID, '%9d\t\t%7.2f\t%7.2f\t%10.5f\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\n', farr');
fclose(fileID);

% Create output file for IP Stresses
farr = zeros(nel*nIntPts, 6);
for i = 1:nel
    row_start = (i-1)*nIntPts+1; row_end = (i-1)*nIntPts+nIntPts;
    farr(row_start:row_end,1) = ones(row_end-row_start+1,1)*elemList(i).elem_num;
    farr(row_start:row_end, 2:end) = elemList(i).getStressIP();
end

fileID = fopen('IPStresses.txt','w');
fprintf(fileID, '|\tElem_Num\t|\tIP_X\t|\tIP_Y\t|\tNode_Str_XX\t|\tNode_Str_YY\t|\tNode_Str_XY\t|\n\n');
fprintf(fileID, '%9d\t\t%11.5f\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\n', farr');
fclose(fileID);














%

