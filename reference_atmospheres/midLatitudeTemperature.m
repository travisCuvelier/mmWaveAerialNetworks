%This function implements the ITU Mid-Latitude reference atmosphere.
%It can be used to estimate temperature at various heights for latitutes
%beween 22° and 45° in either the winter or summer.

%Inputs:
%seasonFlag = 'w' for winter, 's' for summer
%hkm = height (kilometers)

%Output:
%temperatureK = temperature (Kelvin)

%Reference:
% "Reference Standard Atmospheres", ITU Std. ITU-R P.835-5, 2012.

function temperatureK = midLatitudeTemperature( seasonFlag, hkm )

    temperatureK = zeros(size(hkm));
    
    
    
    if(seasonFlag == 's')
        
        s1 = (0<=hkm)&(hkm<13);
        s2 = (13<=hkm)&(hkm<17);
        s3 = (17<=hkm)&(hkm<47);
        s4 = (47<=hkm)&(hkm<53);
        s5 = (53<=hkm)&(hkm<80);
        s6 = (80<=hkm)&(hkm<=100);
        
        temperatureK(s1) = 294.9838 - 5.2159*hkm(s1) -0.07109*hkm(s1).^2;
        temperatureK(s2) = 215.15;
        temperatureK(s3) = 215.15*exp((hkm(s3)-17)*0.008128);
        temperatureK(s4) = 275;
        temperatureK(s5) = 275+(1-exp((hkm(s5)-53)*.06))*20;
        temperatureK(s6) = 175;
        
    elseif(seasonFlag == 'w')
        
        s1 = (0<=hkm)&(hkm<10);
        s2 = (10<=hkm)&(hkm<33);
        s3 = (33<=hkm)&(hkm<47);
        s4 = (47<=hkm)&(hkm<53);
        s5 = (53<=hkm)&(hkm<80);
        s6 = (80<=hkm)&(hkm<=100);
        
        temperatureK(s1) = 272.7241 - 3.6217*hkm(s1) - 0.1759*hkm(s1).^2;
        temperatureK(s2) = 218;
        temperatureK(s3) = 218+(hkm(s3)-33)*3.3571;
        temperatureK(s4) = 265;
        temperatureK(s5) = 265-(hkm(s5)-53)*2.0370;
        temperatureK(s6) = 210;
        
    else
       
        error('Wrong season flag');
        
    end


end

