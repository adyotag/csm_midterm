classdef Element < handle
    % Each element has these attributes
    properties
       elem_num = 0;
       elem_name = NaN;
       elem_type = 0;
       lambda = 0;
       mu = 0;
       elem_D = zeros(3);
       local_num = [];
       positions = [];
       int_pos = NaN;
       int_weights = NaN;
       local_stiffness_matrix = NaN;
       local_loading_vector = NaN;
       local_mass_matrix = NaN;
       local_displacement_vector = NaN;
       local_stress_IP = NaN;
       bforce = zeros(2,1);
       nIntPts = NaN;
       shiftamt = 0;
    end
    methods
        % Constructor for Element object
        function self = Element(elem_name, elem_num, local_num, pos, nen, nIntPts, ps, E, nu, bforce)
            assert(nen==4 | nen==8 | nen==9 , "Improper Element type !");
            assert(length(pos)==2*length(local_num), "Global node number and node position arrays not same length !");
            assert(ps == 1 | ps == 2, "Must choose (1) Plane strain OR (2) Plane stress !")
            assert(nIntPts == 4 | nIntPts == 9, "Only implemented Gaussian quadrature rules for 4 and 9 point integration !")
            self.elem_num = elem_num;
            self.elem_name = elem_name;
            self.nIntPts = nIntPts;
            self.elem_type = nen;            
            self.local_num = local_num;
            self.positions = pos;
            self.rollSF();
            self.bforce(2) = bforce;    % y-component only specified in body force
            self.lambda = E*nu/((1.+nu)*(1.-2.*nu));
            self.mu = E/(2*(1+nu));
            if ps == 1  % Plane Strain case
                self.elem_D(1,1) = self.lambda + 2*self.mu; self.elem_D(1,2) = self.lambda;
                self.elem_D(2,1) = self.lambda; self.elem_D(2,2) = self.lambda + 2*self.mu;
                self.elem_D(3,3) = self.mu;
            elseif ps == 2 % Plane Stress case
                lambda_bar = 2*self.lambda*self.mu / (self.lambda + 2*self.mu);
                self.elem_D(1,1) = lambda_bar + 2*self.mu; self.elem_D(1,2) = lambda_bar;
                self.elem_D(2,1) = lambda_bar; self.elem_D(2,2) = lambda_bar + 2*self.mu;
                self.elem_D(3,3) = self.mu;
            end
            if nIntPts == 4 % 4 point Gauss rule
               self.int_pos = [-(3.)^-0.5, (3.)^-0.5]; self.int_weights = [1., 1.];
            elseif nIntPts == 9 % 9 point Gauss rule
               self.int_pos = [-(3.*0.2)^0.5, 0., (3.*0.2)^0.5]; self.int_weights = [5./9., 8./9.,5./9.];
            end
            self.local_stiffness_matrix = self.generateStiffnessMatrix();
            self.local_loading_vector = self.generateLoadingVector();
            self.local_mass_matrix = self.generateMassMatrix();
            self.local_displacement_vector = zeros(2*self.elem_type);
            self.local_stress_IP = zeros(nIntPts, 3);
        end
        
        % Rolls position and connectivity matrices accordingly to match implementation
        function rollSF(self)
            pos2d = reshape(self.positions, [2, self.elem_type]);
            shiftamt = 0:3; 
            if self.elem_type == 4
               mask = all(pos2d == [min(pos2d(1,:)); min(pos2d(2,:))], 1);
               self.shiftamt = shiftamt(mask);
           elseif self.elem_type == 8
               corners = pos2d(:, 1:4); mids = pos2d(:, 5:8);
               mask = all(corners == [min(corners(1,:)); min(corners(2,:))], 1);
               self.shiftamt = shiftamt(mask);

            elseif self.elem_type == 9
               corners = pos2d(:, 1:4); mids = pos2d(:, 5:8); center = pos2d(:, 9);
               mask = all(corners == [min(corners(1,:)); min(corners(2,:))], 1);
               self.shiftamt = shiftamt(mask);
           end
        end
           
        % Returns the Shape Functions depending on the type of element at hand
        function r = getShapeFunctions(self)
            if self.elem_type == 4 % Correct
                r = @(k, n) circshift( ...
                    0.25*[(1-k)*(1-n), (1+k)*(1-n), (1+k)*(1+n), (1-k)*(1+n)], self.shiftamt);
            elseif self.elem_type == 8 % Correct
                r = @(k, n) circshift( ...
                    [0.25*(k-1)*(1-n)*(1+k+n), 0.25*(1+k)*(n-1)*(1-k+n), 0.25*(1+k)*(1+n)*(-1+k+n), 0.25*(k-1)*(1+n)*(1+k-n),  ...
                     0.5*(1-k^2)*(1-n), 0.5*(1+k)*(1-n^2) , 0.5*(1-k^2)*(1+n), 0.5*(1-k)*(1-n^2) ], ...
                     self.shiftamt);
            elseif self.elem_type == 9 % Correct
                r = @(k, n) circshift( ...
                    [0.25*k*n*(k-1)*(n-1), 0.25*k*n*(k+1)*(n-1), 0.25*k*n*(k+1)*(n+1), 0.25*k*n*(k-1)*(n+1),  ...
                     0.5*n*(1-n)*(k^2 - 1), 0.5*k*(-k-1)*(n^2 - 1) , 0.5*n*(-n-1)*(k^2 - 1), 0.5*k*(1-k)*(n^2 - 1),  ...
                     (1-k^2)*(1-n^2) ], self.shiftamt);
            end
            
        end

        % Returns the derivatives of the Shape Functions depending on the type of element at hand
        function r = getShapeFunctionsD(self)
            if self.elem_type == 4 % Correct
                r = @(k,n) 0.25 * circshift( ...
                    [n-1, 1-n, 1+n, -1-n; k-1, -1-k, 1+k, 1-k], self.shiftamt, 2);
            elseif self.elem_type == 8 % Correct
                r = @(k,n) circshift(  ...
                    [0.25*(1-n)*(2*k+n), 0.25*(1-n)*(2*k-n), 0.25*(1+n)*(2*k+n), 0.25*(1+n)*(2*k-n), k*(n-1), 0.5*(1-n*n), k*(-n-1), 0.5*(-1+n*n); ...
                     0.25*(1-k)*(2*n+k), 0.25*(1+k)*(2*n-k), 0.25*(1+k)*(2*n+k), 0.25*(1-k)*(2*n-k), 0.5*(-1+k*k), n*(-k-1), 0.5*(1-k*k), n*(k-1)], ...
                     self.shiftamt, 2);
            elseif self.elem_type == 9
                r = @(k,n) circshift(  ...
                    [0.25*n*(2*k-1)*(n-1), 0.25*n*(2*k+1)*(n-1), 0.25*n*(2*k+1)*(n+1), 0.25*n*(2*k-1)*(n+1), k*n*(1-n), 0.5*(2*k+1)*(1-n^2), -k*n*(n+1), 0.5*(1-2*k)*(n^2-1), 2*k*(n^2-1); ...
                     0.25*k*(2*n-1)*(k-1), 0.25*k*(k+1)*(2*n-1), 0.25*k*(2*n+1)*(k+1), 0.25*k*(2*n+1)*(k-1), 0.5*(k^2-1)*(1-2*n), -n*k*(k+1), 0.5*(2*n+1)*(1-k^2), n*k*(1-k), 2*n*(k^2-1) ], ...
                     self.shiftamt, 2);
            end
        end
        
        % Returns Jacobian matrix
        function r = getJ(self,k,n)
            % Computation of the jacobian
            SFDH = self.getShapeFunctionsD(); % shape function derivative handle
            SFDM = SFDH(k,n); % shape function derivative matrix  
            temp1 = sum(reshape(SFDM(:)' .* self.positions, [2, self.elem_type]),2);            
            temp2 = sum(reshape(reshape(flip(SFDM, 1),[1, self.elem_type*2]) .* self.positions, [2, self.elem_type]),2);
            r = [temp1(1), temp2(1); temp2(2), temp1(2)];
        end
        
        % Returns B matrix
        function r = getB(self,k,n)
            J = self.getJ(k,n);
            SFDH = self.getShapeFunctionsD(); % shape function derivative handle
            SFDM = SFDH(k,n); % shape function derivative matrix 
            
            % Computation of B
            A = (det(J))^(-1.) * [J(2,2),-J(1,2),0,0; ...
                                  0,0,-J(2,1),J(1,1); ...
                                  -J(2,1),J(1,1),J(2,2),-J(1,2)];

            SFDM_interleaved = reshape([SFDM; zeros(size(SFDM))], [2,2*self.elem_type]);
            G = [SFDM_interleaved; circshift(SFDM_interleaved,1,2)];
            r = A * G;          
        end
        
        % Returns local stiffness matrix of the element
        function r = generateStiffnessMatrix(self)
            r = zeros(2*self.elem_type);
            [kk, nn] = meshgrid(self.int_pos, self.int_pos);
            [wk, wn] = meshgrid(self.int_weights, self.int_weights); 
            ipos = [kk(:),nn(:)];
            iweights = prod([wk(:), wn(:)] ,2);
            
            for iter = 1:length(iweights)
                k = ipos(iter,1); n = ipos(iter,2); w = iweights(iter);
                B = self.getB(k,n); J = self.getJ(k,n);
                r = r + B'*self.elem_D*B*det(J)*w;                
            end
        end
        
        % Returns local loading vector of the element
        function r = generateLoadingVector(self)
            r = zeros(2*self.elem_type,1);
            [kk, nn] = meshgrid(self.int_pos, self.int_pos);
            [wk, wn] = meshgrid(self.int_weights, self.int_weights); 
            ipos = [kk(:),nn(:)];
            iweights = prod([wk(:), wn(:)] ,2);
            SFH = self.getShapeFunctions(); % shape function handle
            
            for iter = 1:length(iweights)
                k = ipos(iter,1); n = ipos(iter,2); w = iweights(iter);
                J = self.getJ(k,n);
                SFM = SFH(k,n); % shape function matrix
                r = r + reshape(repmat(SFM,2,1),[2*self.elem_type,1]) * det(J) * w;

            end
            r = repmat(self.bforce,self.elem_type,1) .* r;
            
        end
        
        % Returns local stiffness matrix set in the global matrix
        function r = getGlobalStiffness(self, m)
                r = zeros(m);
                [L_r, L_c] = meshgrid(1:self.elem_type*2, 1:self.elem_type*2);
                localidx = [L_r(:), L_c(:)];
                
                globalarrpos = 2*(reshape(repmat(self.local_num,2,1),[self.elem_type*2,1]) - 1) + 1 + cast(mod((1:self.elem_type*2) - 1, 2)', 'uint32');
                [G_r, G_c] = meshgrid(globalarrpos, globalarrpos);
                globalidx = [G_r(:), G_c(:)];
                
                localidx1D = sub2ind(size(self.local_stiffness_matrix), localidx(:,1), localidx(:,2));
                globalidx1D = sub2ind(size(r), globalidx(:,1), globalidx(:,2));

                r(globalidx1D) = self.local_stiffness_matrix(localidx1D);
                
        end
        
        % Returns local loading vector in the global vector 
        function r = getGlobalLoading(self, m)
            r = zeros(m,1);  
            globalidx =  2*(reshape(repmat(self.local_num,2,1),[self.elem_type*2,1]) - 1) + 1 + cast(mod((1:self.elem_type*2) - 1, 2)', 'uint32');     
            r(globalidx) = self.local_loading_vector;

        end

        % Updates the local displacements of each element after computation for post processing
        function updateLocalDisplacements(self, d)
            global_idx = 2*(reshape(repmat(self.local_num,2,1),[self.elem_type*2,1]) - 1) + 1 + cast(mod((1:self.elem_type*2) - 1, 2)', 'uint32');
            self.local_displacement_vector = d(global_idx);
            self.local_displacement_vector;
        end
        
        % Generates local mass matrix
        function r = generateMassMatrix(self)
            r = zeros(3*self.elem_type);
            [kk, nn] = meshgrid(self.int_pos, self.int_pos);
            [wk, wn] = meshgrid(self.int_weights, self.int_weights); 
            ipos = [kk(:),nn(:)];
            iweights = prod([wk(:), wn(:)] ,2);
            SFH = self.getShapeFunctions(); % shape function handle
            
            for iter = 1:length(iweights)
                k = ipos(iter,1); n = ipos(iter,2); w = iweights(iter);
                J = self.getJ(k,n);
            
                SFM = SFH(k,n); % shape function matrix
                SFV = reshape(repmat(SFM,3,1),[3*self.elem_type,1]); % shape function vector
                
                r = r + repmat(eye(3), [self.elem_type, self.elem_type]) .* (SFV * SFV') * det(J) * w;
            end
  
        end
        
        % Obtain local mass matrix in the global form
        function r = getGlobalMassMatrix(self, m)
            r = zeros(m);
            
            [L_r, L_c] = meshgrid(1:self.elem_type*3, 1:self.elem_type*3);
            localidx = [L_r(:), L_c(:)];
            
            yyidx = cast(mod((1:self.elem_type*3) - 1, 3)' == 1, 'uint32');
            xyidx = 2*cast(mod((1:self.elem_type*3) - 1, 3)' == 2, 'uint32');
            globalarrpos = 3*(reshape(repmat(self.local_num,3,1),[self.elem_type*3,1]) - 1) + 1 + yyidx + xyidx;
            [G_r, G_c] = meshgrid(globalarrpos, globalarrpos);
            globalidx = [G_r(:), G_c(:)];
            
            localidx1D = sub2ind(size(self.local_mass_matrix), localidx(:,1), localidx(:,2));
            globalidx1D = sub2ind(size(r), globalidx(:,1), globalidx(:,2));

            r(globalidx1D) = self.local_mass_matrix(localidx1D);   
        end
        
        % Obtain local projection vector for each element
        function r = getLocalProjectionVector(self)       
            r = zeros(3*self.elem_type,1);
            [kk, nn] = meshgrid(self.int_pos, self.int_pos);
            [wk, wn] = meshgrid(self.int_weights, self.int_weights); 
            ipos = [kk(:),nn(:)];
            iweights = prod([wk(:), wn(:)] ,2);
            SFH = self.getShapeFunctions(); % shape function handle
            
            for iter = 1:length(iweights)
                k = ipos(iter,1); n = ipos(iter,2); w = iweights(iter);
                J = self.getJ(k,n);
                SFM = SFH(k,n); % shape function matrix
                
                self.local_stress_IP(iter,:) = self.elem_D * self.getB(k,n) * ... 
                    self.local_displacement_vector; % Save stress computation along the way
                
                r = r + reshape(repmat(SFM,3,1),[3*self.elem_type,1]) .* ...
                    reshape(repmat(self.local_stress_IP(iter,:)',self.elem_type,1), [self.elem_type*3,1]) * ...
                    det(J) * w;
            end
        
        end
        
        % Returns local projection vector in the global vector
        function r = getGlobalProjectionVector(self, m)
            r = zeros(m,1);
            yyidx = cast(mod((1:self.elem_type*3) - 1, 3)' == 1, 'uint32');
            xyidx = 2*cast(mod((1:self.elem_type*3) - 1, 3)' == 2, 'uint32');
            globalidx = 3*(reshape(repmat(self.local_num,3,1),[self.elem_type*3,1]) - 1) + 1 + yyidx + xyidx;
            r(globalidx) = self.getLocalProjectionVector();

        end
        
        % Returns local reaction force vector for each element
        function r = generateLocalRFVector(self)
            r = self.local_stiffness_matrix * self.local_displacement_vector ...
                    - self.local_loading_vector;  
        end
        
        % Returns local reaction force vector in global form
        function r = getGlobalRFVector(self, m)
            r = zeros(m,1);
            globalidx =  2*(reshape(repmat(self.local_num,2,1),[self.elem_type*2,1]) - 1) + 1 + cast(mod((1:self.elem_type*2) - 1, 2)', 'uint32');
            r(globalidx) = self.generateLocalRFVector();
        end
        
        % Returns stresses at integration points
        function r = getStressIP(self)
            r = zeros(self.nIntPts, 5);
            [kk, nn] = meshgrid(self.int_pos, self.int_pos);
            r(:,1) = kk(:); r(:,2) = nn(:);
            r(:, 3:5) = self.local_stress_IP;  
        end
       
        % Returns displacement at center of element
        function r = getDispInCtr(self)
            SFH = self.getShapeFunctions(); SFM = SFH(0,0);
            product = repmat(SFM, 2,1) .* reshape(self.local_displacement_vector, [2, self.elem_type]);
            r = sum(product, 2);
        end        
        
        % Returns actual IP points in terms of x and y
        function r = getRealIPs(self)
            r = zeros(self.nIntPts, 2);
            [kk, nn] = meshgrid(self.int_pos, self.int_pos);
            kk = kk(:); nn = nn(:); SFH = self.getShapeFunctions();
            for i = 1:self.nIntPts
                k = kk(i); n = nn(i);
                SFM = repmat(SFH(k,n),2,1);
                r(i, :) = sum(SFM .* reshape(self.positions,[2,self.elem_type]), 2);
            end
            elem_arr = self.elem_num * ones(self.nIntPts,1);
            r = [r, kk, nn, elem_arr];
        end
        
    end
    
end
