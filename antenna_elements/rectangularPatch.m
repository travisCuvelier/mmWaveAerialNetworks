%This object defines a model for a rectangular patch radiator that 
%radiates in x direction, short dimension (width) in z direction, long
%dimension along y direction. See Balanis (antenna theory) 14.2

%See constructor for initialization guidance.

classdef rectangularPatch < handle
    properties
        
        %set
        frequency
        substrateEpsR
        substrateThickness 
        
        %derived
        width
        lambdaFs
        kFs
        physicalLength
        effectiveLength
        deltaL
        effectiveSubstrateEpsR
        
        
    end
    
    
    methods
        
        %Inputs:
        %frequency = operating frequency (Hz)
        %substrateEpsR = relative dielectric permittivity (unitless)
        %subtrateThickness = thickness of substrate (meters)
        
        %Some reasonable parameters would be to set:
        %substrateEpsR = 9.8, which would imply an alumina substrate
        %subtrateThickness = .025*lambda, where lambda is the operating wavelength
  
        function obj = rectangularPatch(frequency,substrateEpsR,subtrateThickness)
            
            c = 3e8;
            obj.frequency = frequency;
            obj.lambdaFs = c/frequency;
            obj.kFs = 2*pi/obj.lambdaFs;
            obj.substrateEpsR = substrateEpsR;
            obj.substrateThickness = subtrateThickness;
            

            obj.width = (c/(2*obj.frequency))*sqrt(2/(obj.substrateEpsR+1));
            obj.effectiveSubstrateEpsR = (obj.substrateEpsR+1)/2+(obj.substrateEpsR-1)/(sqrt(1+12*obj.substrateThickness/obj.width)*2);

            deltaLnum = (obj.effectiveSubstrateEpsR+.3)*(obj.width/obj.substrateThickness+.264);
            deltaLden = (obj.effectiveSubstrateEpsR-.258)*(obj.width/obj.substrateThickness+.8);
            obj.deltaL = obj.substrateThickness*.412*(deltaLnum/deltaLden);
            obj.effectiveLength = c/(sqrt(obj.effectiveSubstrateEpsR)*2*obj.frequency);
            obj.physicalLength = obj.effectiveLength-2*obj.deltaL;

           
        end
        

        %This function evaluates the radient intensity pattern of the
        %constructed patch element at a given azimuth and elevation
        
        %Inputs: 
        %phi = azimuth (radians)
        %theta = elevation (radians)
        
        %Outputs:
        %unNormalizedGain = radiant intensity pattern of patch radiator
        %evaluated at given angles
        
        function unNormalizedGain =  patternSq(obj,phi,theta)
                              
            cosPhi = cos(phi);
            X = (obj.kFs)*(obj.substrateThickness)*sin(theta).*cosPhi/2;
            Z = (obj.kFs)*(obj.width)*cos(theta)/2;
            Y = (obj.kFs)*(obj.effectiveLength)*sin(theta).*sin(phi)/2;
            pat = sin(theta).*sinc(X/pi).*sinc(Z/pi).*cos(Y).*(cosPhi>0);
            unNormalizedGain = conj(pat).*pat;
            
        end
           
    end
end
