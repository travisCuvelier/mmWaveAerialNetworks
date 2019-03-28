classdef uplinkUser < handle
    properties
        
        %set
        
        associatedAP %it assumed that the antennas are being steering towards this ap.
        
        nAntennasY
        
        nAntennasZ
        
        interElementSpacing
        
        cartesianDisplacement %displacement vector from access point to drone wrt standard Cartesian basis. 
        
        orientationBasis %three column unit vectors that specify the orientaiton of the array. The array normal is the second column, by convention.
                         % This matrix defines a coordinate transform. 
                         %When it premultiplies a 
                         %vector defined in terms of the local basis, the
                         %resulting coordinates are cartesian. The affine
                         %offset is not included (obviously).
        
        elementFactorSq %function(phi,theta) handle specifying element gain pattern
        elementTypeStr
        %derived
        
        bfStrategy
        
        range %distance from access point.        
        cartesianToLocalTransform %transform matrix from coorinates wrt. Cartesian basis to local orientation basis
        nAntennas
        frequency %determined by access point
        wavelength
        localDisplacement %displacement from drone to access point in local coordinate system.
        localPhiSteer
        localThetaSteer
        slantAngle
        arrayFilter
        steeredGain %txGain (static, usually);
        rxGain %gainAtAccessPoint (dynamic)
        
        gainDenominatorMatrix
        gainDenominatorDiag
        gainDenominatorEigenvectors 
    end
    
    
    methods
        
        function obj = uplinkUser(associatedAP, displacementFromAP, nAntennasY, nAntennasZ, interElementSpacing, orientationBasis, elementFactorSq,elementTypeStr)
            
            obj.bfStrategy = 'optimal';
            obj.associatedAP = associatedAP;
            associatedAP.associate(obj);

            obj.frequency = obj.associatedAP.frequency;
            obj.wavelength = 3e8/obj.frequency;
            
            obj.nAntennasY = nAntennasY;
            obj.nAntennasZ = nAntennasZ;
            obj.nAntennas = nAntennasY*nAntennasZ;

            obj.interElementSpacing = interElementSpacing;  
            
            obj.elementFactorSq =elementFactorSq;
            
            obj.elementTypeStr = elementTypeStr;
            
            dimString = sprintf('%dX%d_iis_%f_meters_', obj.nAntennasY, obj.nAntennasZ, obj.interElementSpacing);
            fstring = strcat('aerial_nodes/antennaData/gainNormalizationMatrix_',dimString);
            fstring = strcat(fstring,obj.elementTypeStr);
               
            try
                
               load(strcat(fstring,'.mat'));
               obj.gainDenominatorMatrix = gainData.W;
               obj.gainDenominatorDiag = gainData.diag;
               obj.gainDenominatorEigenvectors = gainData.eigs;
               %fprintf('WARNING: USING PRECOMPUTED GAIN PATTERN. ENSURE ELEMENT NAME INCLUDES ALL RELEVANT ELEMENT INFORMATION\n');
               %fprintf(strcat(fstring,'.mat\n\n'));
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
            
            obj.cartesianDisplacement = displacementFromAP;
            obj.range = sqrt(displacementFromAP'*displacementFromAP);
    
            obj.orientationBasis = orientationBasis;
            obj.cartesianToLocalTransform = orientationBasis';

            obj.localDisplacement = -obj.cartesianToLocalTransform*displacementFromAP;
            obj.localPhiSteer = atan2(obj.localDisplacement(2),obj.localDisplacement(1));
            obj.localThetaSteer =  acos(obj.localDisplacement(3)/sqrt(obj.localDisplacement'*obj.localDisplacement));
           
            phiAP = atan2(obj.cartesianDisplacement(2),obj.cartesianDisplacement(1));
            thetaAP =  acos(obj.cartesianDisplacement(3)/sqrt(obj.cartesianDisplacement'*obj.cartesianDisplacement));
            
            obj.slantAngle= pi/2-acos(cos(phiAP)*sin(thetaAP)); %This assumes ap is facing straight down. 
            
            deluey = (0:1:(obj.nAntennasY-1))*obj.interElementSpacing;
            deluey = deluey.';
            deluez = (0:1:(obj.nAntennasZ-1))*obj.interElementSpacing;
            deluez = deluez.';
            arraySteeringVectorY = exp(1i*2*pi*deluey*sin(obj.localPhiSteer)*sin(obj.localThetaSteer)/obj.wavelength);
            arraySteeringVectorZ = exp(1i*2*pi*deluez*cos(obj.localThetaSteer)/obj.wavelength);
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
            
            
            gainNumerator = 4*pi*obj.arrayFilter'*A*obj.arrayFilter*obj.elementFactorSq(obj.localPhiSteer,obj.localThetaSteer);
            gainDenominator = obj.arrayFilter'*obj.gainDenominatorMatrix*obj.arrayFilter;
            obj.steeredGain = real(gainNumerator/gainDenominator);
            
           
        end
           
    end
end
