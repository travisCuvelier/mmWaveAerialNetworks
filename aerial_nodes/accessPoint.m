%maybe this is refactored.
classdef accessPoint < handle
    properties
        nAntennasY;
        nAntennasZ;
        nAntennas;
        frequency;
        interElementSpacing;
        elementFactorSq %(phi,theta)
        elementTypeStr %Be sure and be specific.
        bfStrategy
        
        %derived fields 
        gainDenominatorMatrix
        gainDenominatorDiag
        gainDenominatorEigenvectors
        wavelength;
        ueList = [];
        numUEs
        steeredAtIdx;   
        phiSteer
        thetaSteer
        arrayFilter =[];
        rxGains
    end
    
    
    methods
        
        function obj = accessPoint(frequency, nAntennasY, nAntennasZ,interElementSpacing,elementFactorSq,elementTypeStr)
           
            obj.bfStrategy = 'optimal';
            
            obj.nAntennas = nAntennasY*nAntennasZ;
            obj.nAntennasY = nAntennasY;
            obj.nAntennasZ = nAntennasZ;
            obj.frequency = frequency;
            obj.interElementSpacing = interElementSpacing;
            obj.elementTypeStr = elementTypeStr;
            
            obj.wavelength = 3e8/frequency;
            
            obj.elementFactorSq= elementFactorSq;
            
            obj.numUEs = 0;
            
            dimString = sprintf('%dX%d_iis_%f_meters_', obj.nAntennasY, obj.nAntennasZ, obj.interElementSpacing);
            fstring = strcat('aerial_nodes/antennaData/gainNormalizationMatrix_',dimString);
            fstring = strcat(fstring,obj.elementTypeStr);
               
            try
                
               load(strcat(fstring,'.mat'));
               obj.gainDenominatorMatrix = gainData.W;
    
               obj.gainDenominatorDiag = gainData.diag;
               obj.gainDenominatorEigenvectors = gainData.eigs;
               %warning('Using precomputed gain pattern, ensure element name has all relevant element information.');
            catch
                
                W = zeros(obj.nAntennas,obj.nAntennas);
                k = 2*pi/obj.wavelength;

                
                for row = 1:obj.nAntennas
                   for column = 1:row
                       
                      integrand = @(phi,theta)  exp(1i*k*obj.interElementSpacing*(floor((row-1)/obj.nAntennasZ)*sin(phi).*sin(theta) + mod(row-1,obj.nAntennasZ)*cos(theta)) )...
                            .*exp(-1i*k*obj.interElementSpacing*(floor((column-1)/obj.nAntennasZ)*sin(phi).*sin(theta) + mod(column-1,obj.nAntennasZ)*cos(theta)) ).*obj.elementFactorSq(phi,theta).*sin(theta);
                               
                       W(row,column) =integral2(integrand,0,pi*2,0,pi,'Method','iterated'); %highly oscilliatory integral 
                      
                       
                       if(row~=column)
                           W(column,row) = conj(W(row,column));
                       end
                   end
                end
                
                [U,S,~] = svd(W); 
                
                gainData.W = W;
                gainData.diag = S;
                gainData.eigs = U;
                save(strcat(fstring,'.mat'),'gainData');
                
                obj.gainDenominatorMatrix = W;
                obj.gainDenominatorDiag = S;
                obj.gainDenominatorEigenvectors = U;
    
           end
 
            
        end
        

        function associate(obj, ue)
           
           obj.ueList = [obj.ueList ue]; 
           obj.numUEs = obj.numUEs+1;

        end
        
        function dissociate(obj, ue)
            for idx = 1:obj.numUEs
               if(obj.ueList(idx) == ue)
                   obj.ueList = [obj.ueList(1:(idx-1)) obj.ueList((idx+1):end)];
                   obj.numUEs = obj.numUEs-1;
                   ue.associatedAP = [];
                   break;
               end
            end
        end
        
        function steerAtUE(obj,ueTarget)
            
            %put these displacements back into actual spherical coordinates
            obj.phiSteer = atan2(ueTarget.cartesianDisplacement(2),ueTarget.cartesianDisplacement(1));
            obj.thetaSteer =  acos(ueTarget.cartesianDisplacement(3)/sqrt(ueTarget.cartesianDisplacement'*ueTarget.cartesianDisplacement));
            deluey = (0:1:(obj.nAntennasY-1))*obj.interElementSpacing;
            deluey = deluey.';
            deluez = (0:1:(obj.nAntennasZ-1))*obj.interElementSpacing;
            deluez = deluez.';
            
            arraySteeringVectorY = exp(1i*2*pi*deluey*sin(obj.phiSteer)*sin(obj.thetaSteer)/obj.wavelength);
            arraySteeringVectorZ = exp(1i*2*pi*deluez*cos(obj.thetaSteer)/obj.wavelength);
            arraySteeringVector = kron(arraySteeringVectorY,arraySteeringVectorZ);
            A = arraySteeringVector*arraySteeringVector';
            
            SsqrtInv = diag(1./(diag(obj.gainDenominatorDiag).^(1/2)));
            
            M = SsqrtInv*obj.gainDenominatorEigenvectors'*A*obj.gainDenominatorEigenvectors*SsqrtInv;
            [u,~] = svd(M);
            
            optimalWeights = obj.gainDenominatorEigenvectors*SsqrtInv*u(:,1);
            
            if(strcmp('optimal',obj.bfStrategy))
                obj.arrayFilter = optimalWeights/sqrt(optimalWeights'*optimalWeights);
            else
                obj.arrayFilter = arraySteeringVector;
            end
                  
            obj.rxGains = zeros(obj.numUEs,1);
            
            aidx = 0;
                
            for idx = 1:obj.numUEs
                
                ue = obj.ueList(idx);
                phi = atan2(ue.cartesianDisplacement(2),ue.cartesianDisplacement(1));
                theta =  acos(ue.cartesianDisplacement(3)/sqrt(ue.cartesianDisplacement'*ue.cartesianDisplacement));
                arraySteeringVectorY = exp(1i*2*pi*deluey*sin(phi)*sin(theta)/obj.wavelength);
                arraySteeringVectorZ = exp(1i*2*pi*deluez*cos(theta)/obj.wavelength);
                arraySteeringVector = kron(arraySteeringVectorY,arraySteeringVectorZ);
                A = arraySteeringVector*arraySteeringVector';
                gainNumerator = 4*pi*obj.arrayFilter'*A*obj.arrayFilter*obj.elementFactorSq(phi,theta);
                gainDenominator = obj.arrayFilter'*obj.gainDenominatorMatrix*obj.arrayFilter;
                
                obj.rxGains(idx) = real(gainNumerator/gainDenominator);
                ue.rxGain = obj.rxGains(idx);

                if(ue == ueTarget)
                   aidx = idx;
                end
               
            end
        
            if(aidx == 0)
                error('This ue isn''t associated with this access point. If this is desirable you may handle this exception.');
            end
            
            obj.steeredAtIdx = aidx;

        end
        
        function steerAtUEidx(obj,ueIdx)
            
            if(ueIdx > obj.numUEs)
                error('There''s not that many ue''s associated with this base station.');
            end
            
            ueTarget= obj.ueList(ueIdx);
            obj.steeredAtIdx = ueIdx;
            
            obj.phiSteer = atan2(ueTarget.cartesianDisplacement(2),ueTarget.cartesianDisplacement(1));
            obj.thetaSteer =  acos(ueTarget.cartesianDisplacement(3)/sqrt(ueTarget.cartesianDisplacement'*ueTarget.cartesianDisplacement));
            deluey = (0:1:(obj.nAntennasY-1))*obj.interElementSpacing;
            deluey = deluey.';
            deluez = (0:1:(obj.nAntennasZ-1))*obj.interElementSpacing;
            deluez = deluez.';
            
            arraySteeringVectorY = exp(1i*2*pi*deluey*sin(obj.phiSteer)*sin(obj.thetaSteer)/obj.wavelength);
            arraySteeringVectorZ = exp(1i*2*pi*deluez*cos(obj.thetaSteer)/obj.wavelength);
            arraySteeringVector = kron(arraySteeringVectorY,arraySteeringVectorZ);
            A = arraySteeringVector*arraySteeringVector';
            
            SsqrtInv = diag(1./(diag(obj.gainDenominatorDiag).^(1/2)));
            
            M = SsqrtInv*obj.gainDenominatorEigenvectors'*A*obj.gainDenominatorEigenvectors*SsqrtInv;
            [u,~] = svd(M);
            
            optimalWeights = obj.gainDenominatorEigenvectors*SsqrtInv*u(:,1);
            
            if(strcmp('optimal',obj.bfStrategy))
                obj.arrayFilter = optimalWeights/sqrt(optimalWeights'*optimalWeights);
            else
                obj.arrayFilter = arraySteeringVector;
            end
            
            obj.rxGains = zeros(obj.numUEs,1);
            
                
            for idx = 1:obj.numUEs
       
                ue = obj.ueList(idx);
                phi = atan2(ue.cartesianDisplacement(2),ue.cartesianDisplacement(1));
                theta =  acos(ue.cartesianDisplacement(3)/sqrt(ue.cartesianDisplacement'*ue.cartesianDisplacement));
                arraySteeringVectorY = exp(1i*2*pi*deluey*sin(phi)*sin(theta)/obj.wavelength);
                arraySteeringVectorZ = exp(1i*2*pi*deluez*cos(theta)/obj.wavelength);
                arraySteeringVector = kron(arraySteeringVectorY,arraySteeringVectorZ);
                A = arraySteeringVector*arraySteeringVector';
                gainNumerator = 4*pi*obj.arrayFilter'*A*obj.arrayFilter*obj.elementFactorSq(phi,theta);
                gainDenominator = obj.arrayFilter'*obj.gainDenominatorMatrix*obj.arrayFilter;
                obj.rxGains(idx) = real(gainNumerator/gainDenominator);
                ue.rxGain = obj.rxGains(idx);
            end
            
        end
            
           
    end
end
