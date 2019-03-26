%This function implements the ITU Mid-Latitude reference atmosphere.
%It can be used to estimate barometric pressure at various heights for 
%latitutes beween 22° and 45° in either the winter or summer.

%It's unclear from the reccomendation as to whether this is the total
%atmospheric pressure or just the pressure from "dry air". Since the 
%pressures from dry air usually dominates it probably doesn't make much of
%a difference. A safe bet would be to consider it the "dry pressure", which
%would at worst overestimate the attenuation.

%Inputs:
%seasonFlag = 'w' for winter, 's' for summer
%hkm = height (kilometers)

%Output:
%pressurehPa = pressure (hPa)

%Reference:
% "Reference Standard Atmospheres", ITU Std. ITU-R P.835-5, 2012.

function pressurehPa = midLatitudePressure( seasonFlag, hkm )

    pressurehPa = zeros(size(hkm));
    s1 = (0<=hkm)&(hkm<=10);        
    s2 = (10<hkm)&(hkm<=72);
    s3 = (72<hkm)&(hkm<=100);

    if(seasonFlag == 's')

        
        pressurehPa(s1) = 1012.8186-111.5569*hkm(s1)+3.8646*hkm(s1).^2;
        P10 = 1012.8186-111.5569*10+3.8646*10^2;
        pressurehPa(s2) = P10*exp(-.147*(hkm(s2)-10));
        P72 = P10*exp(-.147*(72-10));
        pressurehPa(s3) = P72*exp(-.165*(hkm(s3)-72));
       
    elseif(seasonFlag == 'w')
        
        pressurehPa(s1) = 1018.8627-124.2954*hkm(s1)+4.8307*hkm(s1).^2;
        P10 = 1018.8627-124.2954*10+4.8307*10.^2;
        pressurehPa(s2) = P10*exp(-.147*(hkm(s2)-10));
        P72 = P10*exp(-.147*(72-10));
        pressurehPa(s3) = P72*exp(-.155*(hkm(s3)-72));
        
    else
       
        error('Wrong season flag');
        
    end


end
