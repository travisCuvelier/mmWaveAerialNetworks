%This object defines a half wave dipole aligned along the z direction. 
%There are no parameters in the construction.
classdef isotropicRadiator< handle
    
    
    methods
        
        function obj = isotropicRadiator()
            
           
        end
        

        %This function evaluates the radient intensity pattern of the
        %constructed isotropic radiating element at a given azimuth and 
        %elevation. Of course, since it's isotropic, it radiates in all
        %directions. 
        
        %Inputs: 
        %phi = azimuth (radians)
        %theta = elevation (radians)
        
        %Outputs:
        %unNormalizedGain = radiant intensity pattern of isotropic radiator
        %evaluated at given angles
        
        function unNormalizedGain =  patternSq(~,~)
             unNormalizedGain = 1;
        end
           
    end
end
