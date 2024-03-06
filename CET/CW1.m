clear all;

%--------------------------Constants------------------------------------%
f = 2.4*10^9; % first frequency (2.4 GHZ) 
C = 3*10^8;   % speed of light 
beta = 2*pi*f/C; % wave number
epsilon_air = 1; % air
epsilon_wall = 3; % wall
snells_constant = sqrt(epsilon_wall/epsilon_air);
bz = 2; bz1 = 2; bz2 = 2; bz3 = 2; bz4 = 2; bz5 =2; %boundary z-coordinate
tx = 0; ty = 4; tz = 2; % transmitter [x,y,z] location
rx = 10; ry = 3; rz = 2; % recevier [x,y,z] location
rxv = rx:.01:22; % from p to q
tx_img = 0; tx1_img = 0; tx2_img = 0; tx3_img = 0; tx4_img = 0; 
tz_img = 2; tz1_img = 2; tz2_img = 2; tz3_img = 2; tz4_img = 2;

%---------------------------------LOS ok-------------------------------------%

% to solve for the direct ray 
d = sqrt((rxv-tx).^2+(ry-ty)^2+(rz-tz)^2); % vector distance between receiver and transmitter
Ed = (1./d).*exp(-j*beta*d); % direct ray (line of sight)
Edf = 20*log10(abs(Ed));

%--------------------------First Order Upper------------------------------%
% to solve for first order reflection (upper boundary)
ty_img = 8; 
by = 6; % y-coordinate of the boundary interception point
C = ty_img;

m = (ty_img-ry)./(tx_img-rxv); % gradient 
bxv = (by-C)./m; % vector x-coordinate of the boundary interception point

d1 = sqrt((bxv-tx).^2+(by-ty)^2+(bz-tz)^2);
d2 = sqrt((bxv-rxv).^2+(by-ry)^2+(bz-rz)^2);
angle_inc = acosd((d1.^2+d2.^2-d.^2)./(2*d1.*d2))/2; 
angle_trans = asind(sind(angle_inc)/snells_constant); 
tau = (cosd(angle_inc)-snells_constant*cosd(angle_trans))./(cosd(angle_inc)+snells_constant*cosd(angle_trans));

E1_upper = (1./(d1+d2)).*exp(-j*beta*(d1+d2)).*tau;

%=========================First Order Lower===============================%
ty_img = -4;
by = 0; % y-coordinate of the boundary interception point
C = ty_img; % y-intercept 

m = (ty_img-ry)./(tx_img-rxv); % gradient 
bxv = (by-C)./m; % vector x-coordinate of the boundary interception point

d1 = sqrt((bxv-tx).^2+(by-ty)^2+(bz-tz)^2); % distance between boundary interception point and the transmitter
d2 = sqrt((bxv-rxv).^2+(by-ry)^2+(bz-rz)^2); % distance between boundary interception point and the receiver
angle_inc = acosd((d1.^2+d2.^2-d.^2)./(2*d1.*d2))/2; % vector incident angle
angle_trans = asind(sind(angle_inc)./snells_constant); % tau angle
tau = (cosd(angle_inc)-snells_constant*cosd(angle_trans))./(cosd(angle_inc)+snells_constant*cosd(angle_trans));

E1_lower = ((1./(d1+d2)).*exp(-j*beta*(d1+d2))).*tau;

%=========================Second Order Upper===============================%
% to solve for first order reflection (upper boundary)
ty1_img = 8; 
ty2_img = -8; 
by1 = 6; by2 = 0;% y-coordinate of the boundary interception point
C1 = ty1_img; C2 = ty2_img;% y-intercept 

m2 = (ty2_img-ry)./(tx2_img-rxv); % gradient 
bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point

d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); 
d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); 
d3 = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz2-rz)^2); 
dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz2-tz)^2);
dB = sqrt((bxv1-rxv).^2+(by1-ry)^2+(bz1-rz)^2);

angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*d1.*d2)))/2; % vector incident angle
angle_trans1 = asind(sind(angle_inc1)./snells_constant); % tau angle
tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1));

angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*d2.*d3)))/2; % vector incident angle
angle_trans2 = asind(sind(angle_inc2)./snells_constant); % tau angle
tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2));


E2_upper = (1./(d1+d2+d3)).*exp(-j*beta*(d1+d2+d3)).*tau1.*tau2;


%=========================Second Order Lower===============================%
ty1_img = -4; 
ty2_img = 16; 
by1 = 0; by2 = 6;% y-coordinate of the boundary interception point
C1 = ty1_img; C2 = ty2_img;% y-intercept 

m2 = (ty2_img-ry)./(tx2_img-rxv); % gradient 
bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point

d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3 = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz2-rz)^2); % distance between boundary interception point and the receiver
dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz2-tz)^2);
dB = sqrt((bxv1-rxv).^2+(by1-ry)^2+(bz1-rz)^2);

angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*(d1.*d2))))/2; % vector incident angle
angle_trans1 = asind(sind(angle_inc1)./snells_constant); % tau angle
tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1));

angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % vector incident angle
angle_trans2 = asind(sind(angle_inc2)./snells_constant); % tau angle
tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2));


E2_lower = (1./(d1+d2+d3)).*exp(-j*beta*(d1+d2+d3)).*tau1.*tau2;

%=========================Third Order Upper===============================%
ty1_img = 8; 
ty2_img = -8; 
ty3_img = 20; 
C1 = ty1_img; C2 = ty2_img; C3 = ty3_img;% y-intercept 
by1 = 6; by2 = 0; by3 = 6;% y-coordinate of the boundary interception point

m3 = (ty3_img-ry)./(tx3_img-rxv); % gradient 
bxv3 = (by3-C3)./m3; % vector x-coordinate of the boundary interception point
m2 = (ty2_img-by3)./(tx2_img-bxv3); % gradient 
bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point

d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3 = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz2-bz3)^2); % distance between boundary interception point and the receiver
d4 = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz3-rz)^2);
dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz2-tz)^2);
dB = sqrt((bxv3-bxv1).^2+(by3-by1)^2+(bz3-bz1)^2);
dC = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz2-rz)^2);

angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*d1.*d2)))/2; % vector incident angle
angle_trans1 = asind(sind(angle_inc1)./snells_constant); % tau angle
tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1));

angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % vector incident angle
angle_trans2 = asind(sind(angle_inc2)./snells_constant); % tau angle
tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2));

angle_inc3 = (acosd((d3.^2+d4.^2-dC.^2)./(2*(d3.*d4))))/2; % angle between the line d1 and d2
angle_trans3 = asind(sind(angle_inc3)./snells_constant); % tau angle
tau3 = (cosd(angle_inc3)-snells_constant*cosd(angle_trans3))./(cosd(angle_inc3)+snells_constant*cosd(angle_trans3));

E3_upper = (1./(d1+d2+d3+d4)).*exp(-j*beta*(d1+d2+d3+d4)).*tau1.*tau2.*tau3;

%=========================Third Order Lower===============================%
ty1_img = -4; 
ty2_img = 16; 
ty3_img = -16; 
C1 = ty1_img; C2 = ty2_img; C3 = ty3_img;% y-intercept 
by1 = 0; by2 = 6; by3 = 0;% y-coordinate of the boundary interception point

m3 = (ty3_img-ry)./(tx3_img-rxv); % gradient 
bxv3 = (by3-C3)./m3; % vector x-coordinate of the boundary interception point
m2 = (ty2_img-by3)./(tx2_img-bxv3); % gradient 
bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point

d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3 = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz2-bz3)^2); % distance between boundary interception point and the receiver
d4 = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz3-rz)^2);
dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz2-tz)^2);
dB = sqrt((bxv3-bxv1).^2+(by3-by1)^2+(bz3-bz1)^2);
dC = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz2-rz)^2);

angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*(d1.*d2))))/2; % angle between the line d1 and d2
angle_trans1 = asind(sind(angle_inc1)./snells_constant); % tau angle
tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1));

angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % angle between the line d1 and d2
angle_trans2 = asind(sind(angle_inc2)./snells_constant); % tau angle
tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2));

angle_inc3 = (acosd((d3.^2+d4.^2-dC.^2)./(2*(d3.*d4))))/2; % angle between the line d1 and d2
angle_trans3 = asind(sind(angle_inc3)./snells_constant); % tau angle
tau3 = (cosd(angle_inc3)-snells_constant*cosd(angle_trans3))./(cosd(angle_inc3)+snells_constant*cosd(angle_trans3));

E3_lower = (1./(d1+d2+d3+d4)).*exp(-j*beta*(d1+d2+d3+d4)).*tau1.*tau2.*tau3;

%=========================Fourth Order Upper===============================%
ty1_img = 8; 
ty2_img = -8; 
ty3_img = 20; 
ty4_img = -20;
C1 = ty1_img; C2 = ty2_img; C3 = ty3_img; C4 = ty4_img;% y-intercept 
by1 = 6; by2 = 0; by3 = 6; by4 = 0;% y-coordinate of the boundary interception point

m4 = (ty4_img - ry)./(tx4_img - rxv); % gradient 
bxv4 = (by4-C4)./m4; % vector x-coordinate of the boundary interception point
m3 = (ty3_img - by4)./(tx3_img - bxv4); % gradient 
bxv3 = (by3-C3)./m3; % vector x-coordinate of the boundary interception point
m2 = (ty2_img-by3)./(tx2_img-bxv3); % gradient 
bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point

d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3 = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz2-bz3)^2); % distance between boundary interception point and the receiver
d4 = sqrt((bxv3-bxv4).^2+(by3-by4)^2+(bz3-bz4)^2);
d5 = sqrt((bxv4-rxv).^2+(by4-ry)^2+(bz4-rz)^2);
dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz2-tz)^2);
dB = sqrt((bxv3-bxv1).^2+(by3-by1)^2+(bz3-bz1)^2);
dC = sqrt((bxv2-bxv4).^2+(by2-by4)^2+(bz2-bz4)^2);
dD = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz3-rz)^2);

angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*(d1.*d2))))/2; % angle between the line d1 and d2
angle_trans1 = asind(sind(angle_inc1)./snells_constant); % tau angle
tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1));

angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % angle between the line d1 and d2
angle_trans2 = asind(sind(angle_inc2)./snells_constant); % tau angle
tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2));

angle_inc3 = (acosd((d3.^2+d4.^2-dC.^2)./(2*(d3.*d4))))/2; % angle between the line d1 and d2
angle_trans3 = asind(sind(angle_inc3)./snells_constant); % tau angle
tau3 = (cosd(angle_inc3)-snells_constant*cosd(angle_trans3))./(cosd(angle_inc3)+snells_constant*cosd(angle_trans3));

angle_inc4 = (acosd((d4.^2+d5.^2-dD.^2)./(2*(d4.*d5))))/2; % angle between the line d1 and d2
angle_trans4 = asind(sind(angle_inc4)./snells_constant); % tau angle
tau4 = (cosd(angle_inc4)-snells_constant*cosd(angle_trans4))./(cosd(angle_inc4)+snells_constant*cosd(angle_trans4));

E4_upper = (1./(d1+d2+d3+d4+d5)).*exp(-j*beta*(d1+d2+d3+d4+d5)).*tau1.*tau2.*tau3.*tau4;

%=========================Fourth Order Lower===============================%
ty1_img = -4; 
ty2_img = 16; 
ty3_img = -16; 
ty4_img = 28;
C1 = ty1_img; C2 = ty2_img; C3 = ty3_img; C4 = ty4_img;% y-intercept 
by1 = 0; by2 = 6; by3 = 0; by4 = 6;% y-coordinate of the boundary interception point

m4 = (ty4_img - ry)./(tx4_img - rxv); % gradient 
bxv4 = (by4-C4)./m4; % vector x-coordinate of the boundary interception point
m3 = (ty3_img - by4)./(tx3_img - bxv4); % gradient 
bxv3 = (by3-C3)./m3; % vector x-coordinate of the boundary interception point
m2 = (ty2_img-by3)./(tx2_img-bxv3); % gradient 
bxv2 = (by2-C2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-C1)./m1; % vector x-coordinate of the boundary interception point

d1 = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2 = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3 = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz2-bz3)^2); % distance between boundary interception point and the receiver
d4 = sqrt((bxv3-bxv4).^2+(by3-by4)^2+(bz3-bz4)^2);
d5 = sqrt((bxv4-rxv).^2+(by4-ry)^2+(bz4-rz)^2);
dA = sqrt((bxv2-tx).^2+(by2-ty)^2+(bz2-tz)^2);
dB = sqrt((bxv3-bxv1).^2+(by3-by1)^2+(bz3-bz1)^2);
dC = sqrt((bxv2-bxv4).^2+(by2-by4)^2+(bz2-bz4)^2);
dD = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz3-rz)^2);

angle_inc1 = (acosd((d1.^2+d2.^2-dA.^2)./(2*(d1.*d2))))/2; % angle between the line d1 and d2
angle_trans1 = asind(sind(angle_inc1)./snells_constant); % tau angle
tau1 = (cosd(angle_inc1)-snells_constant*cosd(angle_trans1))./(cosd(angle_inc1)+snells_constant*cosd(angle_trans1));

angle_inc2 = (acosd((d2.^2+d3.^2-dB.^2)./(2*(d2.*d3))))/2; % angle between the line d1 and d2
angle_trans2 = asind(sind(angle_inc2)./snells_constant); % tau angle
tau2 = (cosd(angle_inc2)-snells_constant*cosd(angle_trans2))./(cosd(angle_inc2)+snells_constant*cosd(angle_trans2));

angle_inc3 = (acosd((d3.^2+d4.^2-dC.^2)./(2*(d3.*d4))))/2; % angle between the line d1 and d2
angle_trans3 = asind(sind(angle_inc3)./snells_constant); % tau angle
tau3 = (cosd(angle_inc3)-snells_constant*cosd(angle_trans3))./(cosd(angle_inc3)+snells_constant*cosd(angle_trans3));

angle_inc4 = (acosd((d4.^2+d5.^2-dD.^2)./(2*(d4.*d5))))/2; % angle between the line d1 and d2
angle_trans4 = asind(sind(angle_inc4)./snells_constant); % tau angle
tau4 = (cosd(angle_inc4)-snells_constant*cosd(angle_trans4))./(cosd(angle_inc4)+snells_constant*cosd(angle_trans4));

E4_lower = (1./(d1+d2+d3+d4+d5)).*exp(-j*beta*(d1+d2+d3+d4+d5)).*tau1.*tau2.*tau3.*tau4;





%plot(rxv,Etest0)
%hold on
%Etest = 20*log10(abs(E3_lower));
%Etest0 = 20*log10(abs(Ed));
%plot(rxv,Etest0)
%hold on
%Etest1 = 20*log10(abs(Ed + E1_upper+E1_lower));
%plot(rxv,Etest1)
%hold on
%Etest2 = 20*log10(abs(Ed+E1_upper+E1_lower+E2_upper+E2_lower));
%plot(rxv,Etest2)
%hold on
Etest3 = 20*log10(abs(Ed+E1_upper+E1_lower+E2_upper+E2_lower+E3_upper+E3_lower));
plot(rxv,Etest3,'r')
hold on
Etest4 = 20*log10(abs(Ed+E1_upper+E1_lower+E2_upper+E2_lower+E3_upper+E3_lower+E4_lower+E4_upper));
plot(rxv,Etest4,'m')
hold off



%xlabel('Line Segment "pq"/m'), ylabel('E-field/dB'), title("LOS");
%legend("LOS");
