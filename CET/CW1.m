clear all;

f = 2.4*10^9; % first frequency (2.4 GHZ) 
c = 3*10^8;   % speed of light 
beta = 2*pi*f/c; % wave number

%==============================LOS==================================%
% to solve for the direct ray 
tx = 0; ty = 4; tz = 2; % transmitter [x,y,z] location
rx = 10; ry = 3; rz = 2; % recevier [x,y,z] location
rxv = rx:.01:22; % from p to q
dv = sqrt((rxv-tx).^2+(ry-ty)^2+(rz-tz)^2); % vector distance between receiver and transmitter
Ed = (1./dv).*exp(-j*beta*dv); % direct ray (line of sight)

Edf = 20*log10(abs(Ed));
plot(rxv,Edf,"r--")
hold on
%xlabel('Line Segment "pq"/m'), ylabel('E-field/dB'), title("LOS");
%legend("LOS");


%=========================First Order Upper===============================%
% to solve for first order reflection (upper boundary)
tx_img = 0; ty_img = 8; tz_img = 2; 
bz=2; %boundary z-coordinate
c = ty_img; % y-intercept 
by = 6; % y-coordinate of the boundary interception point
epsilon_air = 1; % air
epsilon_wall = 3; % wall
snells_constant = sqrt(epsilon_wall/epsilon_air);
m = (ty_img-ry)./(tx_img-rxv); % gradient 
bxv = (by-c)./m; % vector x-coordinate of the boundary interception point
d1v = sqrt((bxv-tx).^2+(by-ty)^2+(bz-tz)^2); % distance between boundary interception point and the transmitter
d2v = sqrt((bxv-rxv).^2+(by-ry)^2+(bz-rz)^2); % distance between boundary interception point and the receiver
angle1v = acosd((d1v.^2+d2v.^2-dv.^2)./(2*(d1v.*d2v))); % angle between the line d1 and d2
angle_incv = angle1v/2; % vector incident angle

angle_tauv = asind(sind(angle_incv)./snells_constant); % tau angle
tauv = (cosd(angle_incv)-snells_constant*cos(angle_tauv))./(cosd(angle_incv)+snells_constant*cos(angle_tauv));
E1_upper = (1./(d1v+d2v)).*exp(-j*beta*(d1v+d2v));

Etest = 20*log10(abs(Ed+E1_upper));
plot(rxv,Etest,'b')
hold on


%=========================First Order Lower===============================%
tx_img = 0; ty_img = -4; tz_img = 2; 
bz=2; %boundary z-coordinate
c = ty_img; % y-intercept 
by = 0; % y-coordinate of the boundary interception point
epsilon_air = 1; % air
epsilon_wall = 3; % wall
snells_constant = sqrt(epsilon_wall/epsilon_air);
m = (ty_img-ry)./(tx_img-rxv); % gradient 
bxv = (by-c)./m; % vector x-coordinate of the boundary interception point
d1v = sqrt((bxv-tx).^2+(by-ty)^2+(bz-tz)^2); % distance between boundary interception point and the transmitter
d2v = sqrt((bxv-rxv).^2+(by-ry)^2+(bz-rz)^2); % distance between boundary interception point and the receiver
angle1v = acosd((d1v.^2+d2v.^2-dv.^2)./(2*(d1v.*d2v))); % angle between the line d1 and d2
angle_incv = angle1v/2; % vector incident angle

angle_tauv = asind(sind(angle_incv)./snells_constant); % tau angle
tauv = (cosd(angle_incv)-snells_constant*cos(angle_tauv))./(cosd(angle_incv)+snells_constant*cos(angle_tauv));
E1_lower = (1./(d1v+d2v)).*exp(-j*beta*(d1v+d2v));
Etest = 20*log10(abs(Ed+E1_upper+E1_lower));
plot(rxv,Etest,'m')
hold on

%=========================Second Order Upper===============================%
% to solve for first order reflection (upper boundary)
tx1_img = 0; ty1_img = 8; tz1_img = 2; 
tx2_img = 0; ty2_img = -8; tz2_img = 2; 
bz1 = 2; bz2 = 2;%boundary z-coordinate
c1 = ty1_img; c2 = ty2_img;% y-intercept 
by1 = 6; by2 = 0;% y-coordinate of the boundary interception point
epsilon_air = 1; % air
epsilon_wall = 3; % wall
snells_constant = sqrt(epsilon_wall/epsilon_air);
m2 = (ty2_img-ry)./(tx2_img-rxv); % gradient 
bxv2 = (by2-c2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-c1)./m1; % vector x-coordinate of the boundary interception point
d1v = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2v = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3v = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz2-rz)^2); % distance between boundary interception point and the receiver
d4v = sqrt((bxv1-rxv).^2+(by1-ry)^2+(bz1-rz)^2);
angle1v = acosd((d2v.^2+d3v.^2-d4v.^2)./(2*(d2v.*d3v))); % angle between the line d1 and d2
angle_incv = angle1v/2; % vector incident angle
angle_tauv = asind(sind(angle_incv)./snells_constant); % tau angle
tauv = (cosd(angle_incv)-snells_constant*cos(angle_tauv))./(cosd(angle_incv)+snells_constant*cos(angle_tauv));
E2_upper = (1./(d1v+d2v+d3v)).*exp(-j*beta*(d1v+d2v+d3v));
%Ef1_upper = 20*log10(abs(E1_upper));
Etest = 20*log10(abs(Ed+E1_upper+E1_lower+E2_upper));
plot(rxv,Etest,'r')
hold on


%=========================Second Order Lower===============================%
tx1_img = 0; ty1_img = -4; tz1_img = 2; 
tx2_img = 0; ty2_img = 16; tz2_img = 2; 
bz1 = 2; bz2 = 2;%boundary z-coordinate
epsilon_air = 1; % air
epsilon_wall = 3; % wall
snells_constant = sqrt(epsilon_wall/epsilon_air);
c1 = ty1_img; c2 = ty2_img;% y-intercept 
by1 = 0; by2 = 6;% y-coordinate of the boundary interception point
m2 = (ty2_img-ry)./(tx2_img-rxv); % gradient 
bxv2 = (by2-c2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-c1)./m1; % vector x-coordinate of the boundary interception point
d1v = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2v = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3v = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz2-rz)^2); % distance between boundary interception point and the receiver
d4v = sqrt((bxv1-rxv).^2+(by1-ry)^2+(bz1-rz)^2);
angle1v = acosd((d2v.^2+d3v.^2-d4v.^2)./(2*(d2v.*d3v))); % angle between the line d1 and d2
angle_incv = angle1v/2; % vector incident angle
angle_tauv = asind(sind(angle_incv)./snells_constant); % tau angle
tauv = (cosd(angle_incv)-snells_constant*cos(angle_tauv))./(cosd(angle_incv)+snells_constant*cos(angle_tauv));
E2_lower = (1./(d1v+d2v+d3v)).*exp(-j*beta*(d1v+d2v+d3v));
%Ef1_upper = 20*log10(abs(E1_upper));
Etest = 20*log10(abs(Ed+E1_upper+E1_lower+E2_upper+E2_lower));
plot(rxv,Etest,'g')
hold on

%=========================Third Order Upper===============================%
tx1_img = 0; ty1_img = 8; tz1_img = 2; 
tx2_img = 0; ty2_img = -8; tz2_img = 2; 
tx3_img = 0; ty3_img = 20; tz3_img = 2; 

bz1 = 2; bz2 = 2; bz3 = 2;%boundary z-coordinate
epsilon_air = 1; % air
epsilon_wall = 3; % wall
snells_constant = sqrt(epsilon_wall/epsilon_air);
c1 = ty1_img; c2 = ty2_img; c3 = ty3_img;% y-intercept 
by1 = 6; by2 = 0; by3 = 6;% y-coordinate of the boundary interception point
m3 = (ty3_img-ry)./(tx3_img-rxv); % gradient 
bxv3 = (by3-c3)./m3; % vector x-coordinate of the boundary interception point
m2 = (ty2_img-by3)./(tx2_img-bxv3); % gradient 
bxv2 = (by2-c2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-c1)./m1; % vector x-coordinate of the boundary interception point
d1v = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2v = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3v = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz2-bz3)^2); % distance between boundary interception point and the receiver
d4v = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz3-rz)^2);
d5v = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz2-rz)^2);
angle1v = acosd((d3v.^2+d4v.^2-d5v.^2)./(2*(d3v.*d4v))); % angle between the line d1 and d2
angle_incv = angle1v/2; % vector incident angle
angle_tauv = asind(sind(angle_incv)./snells_constant); % tau angle
tauv = (cosd(angle_incv)-snells_constant*cos(angle_tauv))./(cosd(angle_incv)+snells_constant*cos(angle_tauv));
E3_upper = (1./(d1v+d2v+d3v+d4v)).*exp(-j*beta*(d1v+d2v+d3v+d4v));
%Ef1_upper = 20*log10(abs(E1_upper));
Etest = 20*log10(abs(Ed+E1_upper+E1_lower+E2_upper+E2_lower+E3_upper));
plot(rxv,Etest,'k')
hold on

%=========================Third Order Lower===============================%
tx1_img = 0; ty1_img = -4; tz1_img = 2; 
tx2_img = 0; ty2_img = 16; tz2_img = 2; 
tx3_img = 0; ty3_img = -16; tz3_img = 2; 

bz1 = 2; bz2 = 2; bz3 = 2;%boundary z-coordinate
epsilon_air = 1; % air
epsilon_wall = 3; % wall
snells_constant = sqrt(epsilon_wall/epsilon_air);
c1 = ty1_img; c2 = ty2_img; c3 = ty3_img;% y-intercept 
by1 = 6; by2 = 0; by3 = 6;% y-coordinate of the boundary interception point
m3 = (ty3_img-ry)./(tx3_img-rxv); % gradient 
bxv3 = (by3-c3)./m3; % vector x-coordinate of the boundary interception point
m2 = (ty2_img-by3)./(tx2_img-bxv3); % gradient 
bxv2 = (by2-c2)./m2; % vector x-coordinate of the boundary interception point
m1 = (ty1_img-by2)./(tx1_img-bxv2); % gradient 
bxv1 = (by1-c1)./m1; % vector x-coordinate of the boundary interception point
d1v = sqrt((bxv1-tx).^2+(by1-ty)^2+(bz1-tz)^2); % distance between boundary interception point and the transmitter
d2v = sqrt((bxv1-bxv2).^2+(by1-by2)^2+(bz1-bz2)^2); % distance between boundary interception point and the receiver
d3v = sqrt((bxv2-bxv3).^2+(by2-by3)^2+(bz2-bz3)^2); % distance between boundary interception point and the receiver
d4v = sqrt((bxv3-rxv).^2+(by3-ry)^2+(bz3-rz)^2);
d5v = sqrt((bxv2-rxv).^2+(by2-ry)^2+(bz2-rz)^2);
angle1v = acosd((d3v.^2+d4v.^2-d5v.^2)./(2*(d3v.*d4v))); % angle between the line d1 and d2
angle_incv = angle1v/2; % vector incident angle
angle_tauv = asind(sind(angle_incv)./snells_constant); % tau angle
tauv = (cosd(angle_incv)-snells_constant*cos(angle_tauv))./(cosd(angle_incv)+snells_constant*cos(angle_tauv));
E3_upper = (1./(d1v+d2v+d3v+d4v)).*exp(-j*beta*(d1v+d2v+d3v+d4v));
%Ef1_upper = 20*log10(abs(E1_upper));
Etest = 20*log10(abs(Ed+E1_upper+E1_lower+E2_upper+E2_lower+E3_upper));
plot(rxv,Etest,'k')
hold on


