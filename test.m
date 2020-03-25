gamma = Element(1, [1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8; 1,2,3,4,5,6,7,8]');
% gamma.getStiffnessMatrix(1)
% gamma.getStiffnessMatrix(2)
% gamma.getStiffnessMatrix(3)
% q = gamma.getShapeFunctions();
% q(-1,-1)
% q(1,-1)
% q(1,1)
% q(-1,1)
% q(0,-1)
% q(1,0)
% q(0,1)
% q(-1,0)
% q(0,0)

% gamma.getStiffnessMatrix(2);
% [node, element, elemType, nel, nen, nIntPts,nnd,ps, nu, E, Force_Node, bforce, disp_BC] = Read_input("Biaxial_Q4_2x2.txt");
% y = gamma.getShapeFunctions()
% z = gamma.getShapeFunctionsD()
% y(1,1)
% z(1,1)

% parfor i = 1:10
%  i  
% end

[node, element, elemType, nel, nen, nIntPts,nnd,ps, nu, E, Force_Node, bforce, disp_BC] = Read_input("Biaxial_Q4_2x2.txt");
