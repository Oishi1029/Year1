clear all;
%--------------------------Coordinates------------------------------------%
tx = 0; ty = 4; tz = 2; % transmitter [x,y,z] location
rx = 10; ry = 3; rz = 2; % recevier [x,y,z] location
rxv = rx:.01:22; % from p to q

%----------------------------Frequency------------------------------------%
f1 = 2.4*10^9; %operation frequency of 2.4 GHz
f2 = 24*10^9; % operation frequency of 24 GHz


%-------------------------------------------------------------------------%
%To determine E-fields for 2.4GHz operating frequency
[ET0, ET1, ET2, ET3, ET4] = determineET(tx, ty, tz, ry, rz, rxv, f1);

%To plot variation graph of E-field from p to q
plotVariationGraph(rxv, f1, ET0, ET1, ET2, ET3, ET4);

%To determine E-fields for 2.4GHz operating frequency
[ET0, ET1, ET2, ET3, ET4] = determineET(tx, ty, tz, ry, rz, rxv, f2);

%To plot variation graph of E-field from p to q
plotVariationGraph(rxv, f2, ET0, ET1, ET2, ET3, ET4);
%-------------------------------------------------------------------------%

function [ET0, ET1, ET2, ET3, ET4] = determineET(tx, ty, tz, ry, rz, rxv, f)

    tx_img = 0; % x-coordinate of all transmitter image
    bz = 2;     % z-coordinate of all boundary interception point

    C = 3*10^8;      % speed of light 
    beta = 2*pi*f/C; % wave number

    epsilon_air = 1;  % air
    epsilon_wall = 3; % wall

    % constant to use in snells formula
    snells_constant = sqrt(epsilon_wall/epsilon_air); 


    %---------------------------------LOS---------------------------------%

    % Vector distance between receiver and transmitter
    d = sqrt((rxv-tx).^2+(ry-ty)^2+(rz-tz)^2);

    % E-field of direct ray (line of sight)
    Ed = (1./d).*exp(-j*beta*d);
    
    %--------------------------First Order Upper--------------------------%

    ty_img = 8;  % y-coordinate of transmitter image
    by = 6;      % y-coordinate of the boundary interception point
    C = ty_img;  % y-intercept 
    
    m = (ty_img-ry)./(tx_img-rxv); % gradient 
    bxv = (by-C)./m; % vector x-coordinate of the boundary interception point
    
    d1 = sqrt((bxv-tx).^2+(by-ty)^2+(bz-tz)^2);  % distance between transmitter and boundary interception point
    d2 = sqrt((bxv-rxv).^2+(by-ry)^2+(bz-rz)^2); % distance between boundary interception point and receiver

    angle_inc = acosd((d1.^2+d2.^2-d.^2)./(2*d1.*d2))/2;  % incident angle
    angle_trans = asind(sind(angle_inc)/snells_constant); % transmission angle
    tau = (cosd(angle_inc)-snells_constant*cosd(angle_trans))...
        ./(cosd(angle_inc)+snells_constant*cosd(angle_trans)); %reflection coefficient 
    
    E1_upper = (1./(d1+d2)).*exp(-j*beta*(d1+d2)).*tau; % E-field 
    
    %----------------------------First Order Lower------------------------%
    ty_img = -4; % y-coordinate of transmitter image
    by = 0;      % y-coordinate of the boundary interception point
    C = ty_img;  % y-intercept 
    
    m = (ty_img-ry)./(tx_img-rxv); % gradient 
    bxv = (by-C)./m; % vector x-coordinate of the boundary interception point
    
    d1 = sqrt((bxv-tx).^2+(by-ty)^2+(bz-tz)^2);  % distance between boundary interception point and the transmitter
    d2 = sqrt((bxv-rxv).^2+(by-ry)^2+(bz-rz)^2); % distance between boundary interception point and the receiver
    angle_inc = acosd((d1.^2+d2.^2-d.^2)./(2*d1.*d2))/2;   % incident angle of first ray
    angle_trans = asind(sind(angle_inc)./snells_constant); % transmission angle of first ray

    tau = (cosd(angle_inc)-snells_constant*cosd(angle_trans))...
        ./(cosd(angle_inc)+snells_constant*cosd(angle_trans)); % reflection coefficient
    
    E1_lower = ((1./(d1+d2)).*exp(-j*beta*(d1+d2))).*tau; % E-field 
    
    %----------------------------Second Order Upper------------------------%
    
    ty1_img = 8; ty2_img = -8;  % y-coordinate of transmitter images
    by1 = 6; by2 = 0;           % y-coordinate of boundary interception points
    C1 = ty1_img; C2 = ty2_img; % y-intercepts 
    
    m2 = (ty2_img-ry)./(tx_img-rxv); % gradient 
    bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
    m1 = (ty1_img-by2)./(tx_img-bxv2); % gradient 
    bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point
    
    d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz-tz)^2);    % distance between boundary interception point and the transmitter
    d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz-bz)^2); % distance between first and second boundary interception point
    d3 = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz-rz)^2);   % distance between second boundary interception point and the receiver
    dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz-tz)^2);    % distance between second boundary interception point and the transmitter
    dB = sqrt((bxv1-rxv).^2+(by1-ry)^2+(bz-rz)^2);   % distance between first boundary interception point and the receiver
    
    angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*d1.*d2)))/2; % incident angle of first ray
    angle_trans1 = asind(sind(angle_inc1)./snells_constant); % transmission angle of first ray
    tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))...
        ./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1)); % reflection coefficient 
    
    angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*d2.*d3)))/2; % incident angle of second ray
    angle_trans2 = asind(sind(angle_inc2)./snells_constant); % transmission angle of second ray
    tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))...
        ./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2)); % reflection coefficient 
    
    
    E2_upper = (1./(d1+d2+d3)).*exp(-j*beta*(d1+d2+d3)).*tau1.*tau2; % E-field 
    
    
    %=========================Second Order Lower===============================%
    ty1_img = -4; ty2_img = 16;  % y-coordinate of transmitter images
    by1 = 0; by2 = 6; % y-coordinate of boundary interception points
    C1 = ty1_img; C2 = ty2_img;% y-intercept 
    
    m2 = (ty2_img-ry)./(tx_img-rxv); % gradient 
    bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
    m1 = (ty1_img-by2)./(tx_img-bxv2); % gradient 
    bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point
    
    d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz-tz)^2); % distance between boundary interception point and the transmitter
    d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz-bz)^2); % distance between first and second boundary interception point
    d3 = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz-rz)^2); % distance between second boundary interception point and the receiver
    dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz-tz)^2); %distance between second boundary interception point and the transmitter
    dB = sqrt((bxv1-rxv).^2+(by1-ry)^2+(bz-rz)^2); %distance between first boundary interception point and the receiver
    
    angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*(d1.*d2))))/2; % incident angle
    angle_trans1 = asind(sind(angle_inc1)./snells_constant); % transmittion angle
    tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))...
        ./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1)); %reflection coefficient
    
    angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % incident angle
    angle_trans2 = asind(sind(angle_inc2)./snells_constant); % transmittion angle
    tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))...
        ./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2)); % reflection coefficient
    
    
    E2_lower = (1./(d1+d2+d3)).*exp(-j*beta*(d1+d2+d3)).*tau1.*tau2; % E-field
    
    %=========================Third Order Upper===============================%
    ty1_img = 8; ty2_img = -8; ty3_img = 20; % y-coordinate of transmitter images
    C1 = ty1_img; C2 = ty2_img; C3 = ty3_img;% y-intercept 
    by1 = 6; by2 = 0; by3 = 6;% y-coordinate of the boundary interception point
    
    m3 = (ty3_img-ry)./(tx_img-rxv); % gradient 
    bxv3 = (by3-C3)./m3; % x-coordinate of the boundary interception point
    m2 = (ty2_img-by3)./(tx_img-bxv3); % gradient 
    bxv2 = (by2-C2)./m2; % x-coordinate of the boundary interception point
    m1 = (ty1_img-by2)./(tx_img-bxv2); % gradient 
    bxv1 = (by1-C1)./m1; % x-coordinate of the boundary interception point
    
    d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz-tz)^2); % distance between first boundary interception point and the transmitter
    d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz-bz)^2); % distance between first and second boundary interception point 
    d3 = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz-bz)^2); % distance between second and third boundary interception point 
    d4 = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz-rz)^2); % distance between third boundary interception point and receiver
    dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz-tz)^2); % distance between second boundary interception point and transmitter
    dB = sqrt((bxv3-bxv1).^2+(by3-by1)^2+(bz-bz)^2); % distance between third and first boundary interception point 
    dC = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz-rz)^2); % distance between second boundary interception point and receiver
    
    angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*d1.*d2)))/2; % incident angle
    angle_trans1 = asind(sind(angle_inc1)./snells_constant); % transmission angle
    tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))...
        ./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1)); %reflection coefficient
    
    angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % incident angle
    angle_trans2 = asind(sind(angle_inc2)./snells_constant); % transmission angle
    tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))...
        ./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2)); %reflection coefficient
    
    angle_inc3 = (acosd((d3.^2+d4.^2-dC.^2)./(2*(d3.*d4))))/2; % incident angle
    angle_trans3 = asind(sind(angle_inc3)./snells_constant); % transmission angle
    tau3 = (cosd(angle_inc3)-snells_constant*cosd(angle_trans3))...
        ./(cosd(angle_inc3)+snells_constant*cosd(angle_trans3));% reflection coefficie
    
    E3_upper = (1./(d1+d2+d3+d4)).*exp(-j*beta*(d1+d2+d3+d4)).*tau1.*tau2.*tau3; % E-field
    
    %=========================Third Order Lower===============================%
    ty1_img = -4; ty2_img = 16; ty3_img = -16; % y-coordinate of transmitter images
    C1 = ty1_img; C2 = ty2_img; C3 = ty3_img;  % y-intercept 
    by1 = 0; by2 = 6; by3 = 0;                 % y-coordinate of the boundary interception point
    
    m3 = (ty3_img-ry)./(tx_img-rxv); % gradient 
    bxv3 = (by3-C3)./m3; % vector x-coordinate of the boundary interception point
    m2 = (ty2_img-by3)./(tx_img-bxv3); % gradient 
    bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
    m1 = (ty1_img-by2)./(tx_img-bxv2); % gradient 
    bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point
    
    d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz-tz)^2); % distance between first boundary interception point and the transmitter
    d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz-bz)^2); % distance between first and second boundary interception point 
    d3 = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz-bz)^2); % distance between second and third boundary interception point 
    d4 = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz-rz)^2); % distance between third boundary interception point and receiver
    dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz-tz)^2); % distance between second boundary interception point and transmitter
    dB = sqrt((bxv3-bxv1).^2+(by3-by1)^2+(bz-bz)^2); % distance between first and third boundary interception point 
    dC = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz-rz)^2); % distance between second boundary interception point and receiver
    
    angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*(d1.*d2))))/2; % incident angle
    angle_trans1 = asind(sind(angle_inc1)./snells_constant); % transmission angle
    tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))...
        ./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1));% reflection coefficient
    
    angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % incident angle
    angle_trans2 = asind(sind(angle_inc2)./snells_constant); % transmission angle
    tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))...
        ./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2)); % reflection coefficient
    
    angle_inc3 = (acosd((d3.^2+d4.^2-dC.^2)./(2*(d3.*d4))))/2; % incident angle
    angle_trans3 = asind(sind(angle_inc3)./snells_constant); % transmission angle
    tau3 = (cosd(angle_inc3)-snells_constant*cosd(angle_trans3))...
        ./(cosd(angle_inc3)+snells_constant*cosd(angle_trans3));% reflection coefficient
    
    E3_lower = (1./(d1+d2+d3+d4)).*exp(-j*beta*(d1+d2+d3+d4)).*tau1.*tau2.*tau3; % E-field
    
    %=========================Fourth Order Upper===============================%
    ty1_img = 8; ty2_img = -8; ty3_img = 20; ty4_img = -20; % y-coordinate of transmitter images
    C1 = ty1_img; C2 = ty2_img; C3 = ty3_img; C4 = ty4_img;% y-intercept 
    by1 = 6; by2 = 0; by3 = 6; by4 = 0;% y-coordinate of the boundary interception point
    
    m4 = (ty4_img - ry)./(tx_img - rxv); % gradient 
    bxv4 = (by4-C4)./m4; % x-coordinate of the boundary interception point
    m3 = (ty3_img - by4)./(tx_img - bxv4); % gradient 
    bxv3 = (by3-C3)./m3; % x-coordinate of the boundary interception point
    m2 = (ty2_img-by3)./(tx_img-bxv3); % gradient 
    bxv2 = (by2-C2)./m2; % x-coordinate of the boundary interception point
    m1 = (ty1_img-by2)./(tx_img-bxv2); % gradient 
    bxv1 = (by1-C1)./m1; % x-coordinate of the boundary interception point
    
    d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz-tz)^2); % distance between first boundary interception point and the transmitter
    d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz-bz)^2); % distance between first and second boundary interception point
    d3 = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz-bz)^2); % distance between second and third boundary interception point
    d4 = sqrt((bxv3-bxv4).^2+(by3-by4)^2+(bz-bz)^2); % distance between third and fourth boundary interception point
    d5 = sqrt((bxv4-rxv).^2+(by4-ry)^2+(bz-rz)^2);   % distance between fourth boundary interception point and receiver
    dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz-tz)^2);    % distance between second boundary interception point and transmitter
    dB = sqrt((bxv3-bxv1).^2+(by3-by1)^2+(bz-bz)^2); % distance between third and first boundary interception point 
    dC = sqrt((bxv2-bxv4).^2+(by2-by4)^2+(bz-bz)^2); % distance between second and fourth boundary interception point 
    dD = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz-rz)^2);   % distance between third boundary interception point and receiver
    
    angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*(d1.*d2))))/2; % incident angle
    angle_trans1 = asind(sind(angle_inc1)./snells_constant); % transmission angle
    tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))...
        ./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1)); % reflection coefficient
    
    angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % incident angle
    angle_trans2 = asind(sind(angle_inc2)./snells_constant); % transmission angle
    tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))...
        ./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2));% reflection coefficient
    
    angle_inc3 = (acosd((d3.^2+d4.^2-dC.^2)./(2*(d3.*d4))))/2; % incident angle
    angle_trans3 = asind(sind(angle_inc3)./snells_constant); % transmission angle
    tau3 = (cosd(angle_inc3)-snells_constant*cosd(angle_trans3))...
        ./(cosd(angle_inc3)+snells_constant*cosd(angle_trans3));% reflection coefficient
    
    angle_inc4 = (acosd((d4.^2+d5.^2-dD.^2)./(2*(d4.*d5))))/2; % incident angle
    angle_trans4 = asind(sind(angle_inc4)./snells_constant); % transmission angle
    tau4 = (cosd(angle_inc4)-snells_constant*cosd(angle_trans4))...
        ./(cosd(angle_inc4)+snells_constant*cosd(angle_trans4));% reflection coefficient
    
    E4_upper = (1./(d1+d2+d3+d4+d5)).*exp(-j*beta*(d1+d2+d3+d4+d5)).*tau1.*tau2.*tau3.*tau4; % E-field
    
    %=========================Fourth Order Lower===============================%
    ty1_img = -4; ty2_img = 16; ty3_img = -16; ty4_img = 28; % y-coordinate of transmitter images
    C1 = ty1_img; C2 = ty2_img; C3 = ty3_img; C4 = ty4_img;% y-intercept 
    by1 = 0; by2 = 6; by3 = 0; by4 = 6;% y-coordinate of the boundary interception point
    
    m4 = (ty4_img - ry)./(tx_img - rxv); % gradient 
    bxv4 = (by4-C4)./m4; % x-coordinate of the boundary interception point
    m3 = (ty3_img - by4)./(tx_img - bxv4); % gradient 
    bxv3 = (by3-C3)./m3; % x-coordinate of the boundary interception point
    m2 = (ty2_img-by3)./(tx_img-bxv3); % gradient 
    bxv2 = (by2-C2)./m2; % x-coordinate of the boundary interception point
    m1 = (ty1_img-by2)./(tx_img-bxv2); % gradient 
    bxv1 = (by1-C1)./m1; % x-coordinate of the boundary interception point
    
    d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz-tz)^2); % distance between boundary interception point and the transmitter
    d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz-bz)^2); % distance between first and second boundary interception point
    d3 = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz-bz)^2); % distance between second and third boundary interception point 
    d4 = sqrt((bxv3-bxv4).^2+(by3-by4)^2+(bz-bz)^2); % distance between third and fourth boundary interception point 
    d5 = sqrt((bxv4-rxv).^2+(by4-ry)^2+(bz-rz)^2);   % distance between fourth boundary interception point and the receiver
    dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz-tz)^2);    % distance between second boundary interception point and the transmitter
    dB = sqrt((bxv3-bxv1).^2+(by3-by1)^2+(bz-bz)^2); % distance between third and first boundary interception point 
    dC = sqrt((bxv2-bxv4).^2+(by2-by4)^2+(bz-bz)^2); % distance between second and fourth boundary interception point 
    dD = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz-rz)^2);   % distance between third boundary interception point and receiver
    
    angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*(d1.*d2))))/2; % incident angle
    angle_trans1 = asind(sind(angle_inc1)./snells_constant); % transmission angle
    tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1));
    
    angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % incident angle
    angle_trans2 = asind(sind(angle_inc2)./snells_constant); % transmission angle
    tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))...
        ./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2)); % reflection coefficient
    
    angle_inc3 = (acosd((d3.^2+d4.^2-dC.^2)./(2*(d3.*d4))))/2; % incident angle
    angle_trans3 = asind(sind(angle_inc3)./snells_constant); % transmission angle
    tau3 = (cosd(angle_inc3)-snells_constant*cosd(angle_trans3))...
        ./(cosd(angle_inc3)+snells_constant*cosd(angle_trans3)); % reflection coefficient 
    
    angle_inc4 = (acosd((d4.^2+d5.^2-dD.^2)./(2*(d4.*d5))))/2; % incident angle
    angle_trans4 = asind(sind(angle_inc4)./snells_constant); % transmission angle
    tau4 = (cosd(angle_inc4)-snells_constant*cosd(angle_trans4))...
        ./(cosd(angle_inc4)+snells_constant*cosd(angle_trans4)); % reflection coefficient 
    
    E4_lower = (1./(d1+d2+d3+d4+d5)).*exp(-j*beta*(d1+d2+d3+d4+d5)).*tau1.*tau2.*tau3.*tau4; %E-field

    %--------------------Calculate E-field Summation----------------------%
    ET0 = 20*log10(abs(Ed));
    ET1 = 20*log10(abs(Ed + E1_upper + E1_lower));
    ET2 = 20*log10(abs(Ed + E1_upper + E1_lower + E2_upper));
    ET3 = 20*log10(abs(Ed + E1_upper + E1_lower + E2_upper + E2_lower ...
        + E3_upper + E3_lower));
    ET4 = 20*log10(abs(Ed + E1_upper + E1_lower + E2_upper + E2_lower ...
        + E3_upper + E3_lower + E4_upper + E4_lower));
end

function plotVariationGraph(rxv, f, ET0, ET1, ET2, ET3, ET4)
    
    % to create figure
    figure
    
    % to plot LOS 
    plot(rxv,ET0,'b--')

    %to retain figure
    hold on
    
    % to plot LOS + 1st ref
    plot(rxv,ET1,'m--')

    % to plot LOS + 1st ref + 2nd ref
    plot(rxv,ET2,'k--')

    % to plot LOS + 1st ref + 2nd ref + 3rd ref
    plot(rxv,ET3,'g--')

    % to plot LOS + 1st ref + 2nd ref + 3rd ref + 4th ref
    plot(rxv,ET4, 'r-')
    
    %to turn off retain figure
    hold off

    %to label each E-field summation
    legend("ET0","ET1","ET2","ET3","ET4")
    

    if f == 2.4*10^9
        title('Variation of signal strength of 2.4 GHz E-field')
    else
        title('Variation of signal strength of 24 GHz E-field')
    end

    xlabel('Line Segment "pq"/m')
    ylabel('E-field/dB')
    
end
