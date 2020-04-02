function [node, element, elemType, nel, nen, nIntPts,nnd,ps, nu, E, Force_Node, bforce, disp_BC]=Read_input(filename)

% you can remove the two cd ../ line if your files are in the same
% directories
%cd ../file
filename
fid=fopen(filename, 'r');
%cd ../Solution

C=textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
E=C{1,1};

C=textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nu=C{1,1};

C=textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
ps=C{1,1};

C=textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nnd=C{1,1};
C=textscan(fid,'%n %n','MultipleDelimsAsOne', 1);
X=C{1};
Y=C{2};
node=[X Y];

C=textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nel=C{1,1};
C=textscan(fid,'%n',1);
nen=C{1,1};
if (nen==4)
    elemType='Q4';
    C=textscan(fid,'%u %u %u %u','MultipleDelimsAsOne', 1);
    element=cell2mat(C(1:4));
elseif (nen==9)
    elemType='Q9';
    C=textscan(fid,'%u %u %u %u %u %u %u %u %u','MultipleDelimsAsOne', 1);
    element=cell2mat(C(1:9));
elseif (nen==8)
    elemType='Q8';
    C=textscan(fid,'%u %u %u %u %u %u %u %u','MultipleDelimsAsOne', 1);
    element=cell2mat(C(1:8));
end

C=textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nIntPts=C{1,1};

C=textscan(fid,'%s',1);
C=textscan(fid,'%s',6);
C=textscan(fid,'%n %n %n','MultipleDelimsAsOne', 1);
clear disp_BC
z=C{3};
disp_BC=[C{1} C{2} z];

C=textscan(fid,'%s',4);
C=textscan(fid,'%u','MultipleDelimsAsOne', 1);
Force_Node=C{1};
C=textscan(fid,'%s',2);
C=textscan(fid,'%n',1);
bforce=C{1,1};
fclose(fid);

