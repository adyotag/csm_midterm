% Main function that can be called by the report.m
classdef Assembly < handle
    properties
        % file attributes
        filetype = NaN;
        elemList = [];
        node = NaN; element= NaN; elemType = NaN;
        nel = 0; nen = 0; nIntPts = 0; nnd = 0;
        ps = 0; nu = 0; E = 0; Force_Node = NaN;
        bforce = 0; disp_BC = NaN;

        % assembly attributes
        trf = zeros(2);

    end
    methods
        function self = Assembly(f)
            % Read in filetype
            self.filetype = f;
            [self.node, self.element, self.elemType, self.nel, self.nen, ...
                self.nIntPts, self.nnd, self.ps, self.nu, self.E, self.Force_Node, ...
                self.bforce, self.disp_BC] = Read_input(self.filetype);
        end
        
        function r = run(self)
            % Construct each element. Preallocating space for objects seems like a rabbit hole ;-)
            for i = 1:self.nel
                local_num = self.element(i,:);
                local_pos = reshape(self.node(local_num,:)',1,[]);
                self.elemList = [self.elemList, Element(self.elemType, i, local_num, ...
                    local_pos, self.nen, self.nIntPts, self.ps, self.E, self.nu, self.bforce)];
            end

            % Perform assembly of Global Stiffness Matrix
            k_global = zeros(self.nnd*2);
            f_global = zeros(self.nnd*2,1);
            for i = 1:self.nel
                k_global = k_global + self.elemList(i).getGlobalStiffness(self.nnd*2);
                f_global = f_global + self.elemList(i).getGlobalLoading(self.nnd*2);
            end

            % Apply Displacement BC
            globaleqns = 2*(self.disp_BC(:,1) - 1) + self.disp_BC(:,2);
            dispamt = self.disp_BC(:,3);
            for i = 1:length(dispamt)
                idx = globaleqns(i);
                f_global = f_global - k_global(:,idx)*dispamt(i);
                k_global(:,idx) = zeros(self.nnd*2,1); k_global(idx,:) = zeros(self.nnd*2,1);
                k_global(idx,idx) = 1.;
            end
            
            % Reduce condition of matrix
            mask = ~logical(sum(globaleqns == (1:self.nnd*2),1));
            k_global = k_global(:, mask); k_global = k_global(mask, :);
            f_global = f_global(mask);

            % Reassemble final solution
            d_global = zeros(2*self.nnd,1);
            d_global(mask) = k_global\f_global;
            d_global(globaleqns) = dispamt;
            
            % Update local displacements for each element
            for i = 1:self.nel
               self.elemList(i).updateLocalDisplacements(d_global); 
            end

            % Compute Nodal Stresses
            m_global = zeros(3*self.nnd);
            p_global = zeros(3*self.nnd,1);
            for i = 1:self.nel
                m_global = m_global + self.elemList(i).getGlobalMassMatrix(3*self.nnd);
                p_global = p_global + self.elemList(i).getGlobalProjectionVector(3*self.nnd);
            end
            n_stress = m_global\p_global;

            % Compute Reaction Forces
            rf_global = zeros(2*self.nnd,1);
            num_rf = 2*length(self.Force_Node);
            for i = 1:self.nel
                rf_global = rf_global + self.elemList(i).getGlobalRFVector(2*self.nnd);
            end

            rfgidx = 2*(reshape(repmat(self.Force_Node',2,1),[num_rf,1]) - 1) + 1 + cast(mod((1:num_rf) - 1, 2)', 'uint32');
            rf_global = rf_global(rfgidx);

            % Compute Total Reaction Force
            self.trf = sum(reshape(rf_global, [2, num_rf/2]),2)';

            % Create output file for nodal information
            farr = zeros(self.nnd, 8);
            farr(:,1) = 1:self.nnd; farr(:,2:3) = self.node;
            farr(:, 4:5) = reshape(d_global',[2,self.nnd])';
            farr(:, 6:end) = reshape(n_stress',[3,self.nnd])';
            % Why is importing textfiles so damn difficult in MATLAB?
            save("NodalStressAndDisp_"+extractBefore(self.filetype,".")+".mat", 'farr');
            
            fileID = fopen('NodalStressAndDisp_'+self.filetype,'w');
            fprintf(fileID, '|\tNode_Num\t|\tX\t|\tY\t|\tX_Disp\t|\tY_Disp\t|\tNode_Str_XX\t|\tNode_Str_YY\t|\tNode_Str_XY\t|\n\n');
            fprintf(fileID, '%9f\t\t%7.2f\t%7.2f\t%10.5f\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\n', farr');
            fclose(fileID);

            % Create output file for IP Stresses
            farr = zeros(self.nel*self.nIntPts, 6);
            for i = 1:self.nel
                row_start = (i-1)*self.nIntPts+1; row_end = (i-1)*self.nIntPts+self.nIntPts;
                farr(row_start:row_end,1) = ones(row_end-row_start+1,1)*self.elemList(i).elem_num;
                farr(row_start:row_end, 2:end) = self.elemList(i).getStressIP();
            end
            % Why is importing textfiles so damn difficult in MATLAB?
            save("IPStresses_" + extractBefore(self.filetype,".") + ".mat", 'farr'); 
            fileID =fopen("IPStresses_"+self.filetype,'w');
            fprintf(fileID, '|\tElem_Num\t|\tIP_X\t|\tIP_Y\t|\tNode_Str_XX\t|\tNode_Str_YY\t|\tNode_Str_XY\t|\n\n');
            fprintf(fileID, '%9f\t\t%11.5f\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\n', farr');
            fclose(fileID);
            
            r = self;            
        end
        
        % Function to help read output files. Useful keeping report.m clean
        function r = readout(self)
            r = Read_output(self.filetype);
        end
        
    end
end