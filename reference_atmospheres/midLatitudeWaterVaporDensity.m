%This function implements the ITU Mid-Latitude reference atmosphere.
%It can be used to estimate water vapor density at various heights for 
%latitutes beween 22° and 45° in either the winter or summer.

%Inputs:
%seasonFlag = 'w' for winter, 's' for summer
%hkm = height (kilometers)

%Output:
%waterVaporDensitygpm3 = water vapor density (grams per cubic meter)

%Reference:
% "Reference Standard Atmospheres", ITU Std. ITU-R P.835-5, 2012.
function waterVaporDensitygpm3 = midLatitudeWaterVaporDensity( seasonFlag, hkm )

    waterVaporDensitygpm3 = zeros(size(hkm));

    if(seasonFlag == 's')
        
        s1 = hkm<=15;
        
        waterVaporDensitygpm3(s1) = 14.3542*exp(-.4174*hkm(s1)...
            -.02290*hkm(s1).^2+.001007*hkm(s1).^3);
        
        
    elseif(seasonFlag == 'w')
        
        s1 = hkm<=10;
        waterVaporDensitygpm3(s1) = 3.4742*exp(-.2697*hkm(s1)...
            -.03604*hkm(s1).^2+.0004489*hkm(s1).^3);
        
    else
       
        error('Wrong season flag');
        
    end


end

