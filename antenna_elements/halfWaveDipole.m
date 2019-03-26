%This object defines a half wave dipole aligned along the z direction. 
%There are no parameters in the construction.
classdef halfWaveDipole< handle
    properties
        
        electricalLength
        
    end
    
    
    methods
        
        function obj = halfWaveDipole()
            
            obj.electricalLength = .5;
           
        end
        

        %This function evaluates the radient intensity pattern of the
        %constructed half wave dipole element at a given azimuth and 
        %elevation
        
        %Inputs: 
        %phi = azimuth (radians)
        %theta = elevation (radians)
        
        %Outputs:
        %unNormalizedGain = radiant intensity pattern of dipole
        %evaluated at given angles
        
        function unNormalizedGain =  patternSq(~,theta)
             unNormalizedGain = (cos(pi*cos(theta)/2)./sin(theta)).^2; %z oriented
        end
           
    end
end
