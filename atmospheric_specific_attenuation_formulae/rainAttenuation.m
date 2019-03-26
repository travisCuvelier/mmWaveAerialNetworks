%From ITU-R P.838-3 
%attenuationdBperKm gives attenuation in dB per kilometers
%for horizontally (1 index) and vertically (2 index) waves.
%Approximately valid from 1 to 1000 GHz.
%Theta is elevation angle, tau is polarization wrt. horizontal (tau = pi/4
%for circular polarization). This has been verified.

%This script computes the specific attenuation (dB/km) caused by rain
%between 1 and 1000 GHz. The attenuation is negative, by convention. 

%Inputs:
%fHz = operating frequency (Hertz) 
%RainRatemmph = rate of rainfall in millimeters per hour
%theta = path elevation angle (radians)
%tau = polarization angle with respect to the horizontal (radians) 
%For circular polarization, tau should be set to pi/4.
%Note that there is no dependence on path elevation for circular
%polarization.

%Outputs:
%attenuationdBperKm = specific attenuation due to rainfall (dB/kilometers)
%Note that attenuationdBperKm will be negative.  

function attenuationdBperKm = rainAttenuation( fHz, RainRatemmph,theta,tau )

    if(sum(fHz<10^9))
        error('Not valid below 1 GHz');
    elseif(sum(fHz>1000*10^9))
        error('Not valid above 1000 GHz');
    end
    
    f(1,:) = fHz/(10^9);

    %index one is subscript used in sum.
    %Index 2 is polarization horizontal (1) or vertical (2). 
    ak = zeros(4,2);
    bk= zeros(4,2);
    ck = zeros(4,2);
    mk = zeros(1,2);
    scriptck = zeros(1,2);

    aa = zeros(5,2);
    ba= zeros(5,2);
    ca = zeros(5,2);
    ma = zeros(1,2);
    scriptca = zeros(1,2);

    ak(:,1) = [-5.33980 -0.35351 -0.23789 -0.94158].';
    bk(:,1) = [-0.10008 1.26970 0.86036 0.64552].';
    ck(:,1) = [1.13098 0.45400 0.15354 0.16817].';
    mk(:,1) = -0.18961;
    scriptck(:,1) = 0.71147;

    ak(:,2) = [-3.80595 -3.44965 -0.39902 0.50167].';
    bk(:,2) = [0.56934 -0.22911 0.73042 1.07319].';
    ck(:,2) = [0.81061 0.51059 0.11899 0.27195].';
    mk(:,2) = -0.16398;
    scriptck(:,2) = 0.63297;

    aa(:,1) = [-0.14318 0.29591 0.32177 -5.37610 16.1721].';
    ba(:,1) = [1.82442 0.77564 0.63773 -0.96230 -3.29980].';
    ca(:,1) = [-0.55187 0.19822 0.13164 1.47828 3.43990].';
    ma(:,1) = 0.67849;
    scriptca(:,1) = -1.95537;

    aa(:,2) = [-0.07771 0.56727 -0.20238 -48.2991 48.5833].';
    ba(:,2) = [2.33840 0.95545 1.14520 0.791669 0.791459].';
    ca(:,2) = [-0.76284 0.54039 0.26809 0.116226 0.116479].';
    ma(:,2) = -0.053739;
    scriptca(:,2) = 0.83433;



    lgks = zeros(2,length(f));
    as = zeros(2,length(f));


    for hvidx = 1:2

        insideExp = -((log10(f) - bk(:,hvidx))./ck(:,hvidx)).^2;
        lgks(hvidx,:) = sum(ak(:,hvidx).*exp(insideExp),1)+mk(:,hvidx)*log10(f)+scriptck(:,hvidx);

    end
    
    ks = 10.^lgks;
   
    for hvidx = 1:2

        insideExp = -((log10(f) - ba(:,hvidx))./ca(:,hvidx)).^2;
        as(hvidx,:) = sum(aa(:,hvidx).*exp(insideExp),1)+ma(:,hvidx)*log10(f)+scriptca(:,hvidx);

    end

    kh = ks(1,:);
    kv = ks(2,:);
    ah = as(1,:);
    av = as(2,:);

    
    k = (kh+kv+(kh-kv)*(cos(theta))^2*cos(2*tau))/2;
    a = (kh.*ah+kv.*av+(kh.*ah-kv.*av)*(cos(theta))^2*cos(2*tau))./(2*k);
 
    attenuationdBperKm = -k.*RainRatemmph.^(a);
   
end

