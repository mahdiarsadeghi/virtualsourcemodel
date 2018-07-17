%% Mahdiar Sadeghi
%% Northeastern University
%% December 6th 2017
%
% Inputs : Desired location (x0, y0), metal panel distance to the panel H
% dielectric loss, dielectric permittivity eps, and dielectric thickness d
%
% notes:
% -(x0, y0)=(0,0) is at the center of reflector panel
% -all distance values should be in cm
% -to have lossless dielectric, use loss=0
% -for bare metal case, use eps=1 and d=0
%
% Outputs : E response singal as a function z (distance from the detector)
%
%
% ease case:
% try this for bare metal: [E, z] = VSModel3D(0, 0, 0, 1, 0, 80);
% 

function [E, z] = VSModel3D(x0, y0, loss, eps, d, H) 

%% Global Variables
c = 3e10; f=24.16e9;
k = 2*pi*f/c;
n = sqrt(eps);
gam = (1-n)/(1+n); 
n1 = -(1/.96)*n;

% The exact (x,y) location of each antenna element on the eqo panel
xtile=0:6.248:125; xoff = 145.5 * (0:6);
x= [xtile+xoff(1), xtile+xoff(2), xtile+xoff(3), xtile+xoff(4), ...
    xtile+xoff(5), xtile+xoff(6), xtile+xoff(7)];
x=round(x/10, 2); 

ytile=0:6.248:294; yoff = 952+ 314.2 * (0:2);
y= [ytile+yoff(1), ytile+yoff(2), ytile+yoff(3)];
y=round(y/10, 2);

% Shifting x, y to have the origin at the center of the panel
xav = 49.9;  % (max(x)+min(x)) /2 = 49.9
yav = 141.3; % (max(y)+min(y)) /2 = 141.3
x = x - xav;
y = y - yav;  

%% Removing non-effective elements
x(x<(min(x)-2*x0)) = [];
x(x>(max(x)-2*x0)) = [];
y(y<(min(y)-2*y0)) = [];
y(y>(max(y)-2*y0)) = [];

%% Building 2D (x, y) location of each antenna element
[~,xsize] = size(x);    
[~,ysize] = size(y);    

X = ones(ysize,1)*x;
Y = y'*ones(1,xsize);

% Distance of effective antenna elements from the focal point in xy plane
R = (X-x0).^2+(Y-y0).^2;

%% depth axis resolution per cm
f = -20:.01:20;
[~, fsize] = size(f);

%% Initialize final matrix
E = zeros(size(fsize));

% Loss factors
loss2 =  exp(-k*2*d*n*loss);
loss3 =  exp(-k*4*d*n*loss);
loss4 =  exp(-k*6*d*n*loss);
loss5 =  exp(-k*8*d*n*loss);

for counter=1:fsize
    
    %% Receiver/transmitter adjusted phase(as function of path length)
    l0 = sqrt(R + (H-f(counter))^2 * ones(ysize,xsize));
%     l0 = round(k.*sqrt(R + (H-f(counter))^2 * ones(ysize,xsize))./pi)*pi/k;
    
    %% Path length derivations
    l1 = sqrt(R + (-2*d+H+f(counter))^2 * ones(ysize,xsize));
    l2 = sqrt(R + (-2*d*(1+1/n1)+H+f(counter))^2 * ones(1,xsize));
    l3 = sqrt(R + (-2*d*(1+2/n1)+H+f(counter))^2 * ones(1,xsize));
    l4 = sqrt(R + (-2*d*(1+3/n1)+H+f(counter))^2 * ones(1,xsize));
    l5 = sqrt(R + (-2*d*(1+4/n1)+H+f(counter))^2 * ones(1,xsize));
      
    %% Path length effect
    L1 = (1./l1);
    L2 = (1./l2);  
    L3 = (1./l3);
    L4 = (1./l4);
    L5 = (1./l5);
    
    %% Final result
    comp1 = gam*sum(sum(L1 .* exp(-1i*k*(-l0+l1 ))));
    comp2 = loss2* (-1+gam^2)*sum( sum((L2 .* exp(-1i*k*(2*d*(n+1/n1)*ones(1,xsize)  -l0 + l2 )))));
    comp3 = loss3* gam*(-1+gam^2)*sum( sum((L3 .* exp(-1i*k*(4*d*(n+(1/n1))*ones(1,xsize)  -l0 + l3 )))));
    comp4 = loss4* gam^2 * (-1+gam^2)*sum( sum(L4 .* exp(-1i*k*(6*d*(n+(1/n1))*ones(1,xsize)  -l0 + l4 ))));
    comp5 = loss5* gam^3 * (-1+gam^2)*sum( sum(L5 .* exp(-1i*k*(8*d*(n+(1/n1))*ones(1,xsize)  -l0 + l5 ))));
    
 
    E (counter) = comp1+comp2+comp3+comp4+comp5;
end

%% normalize the results
z = H-f;
E = (100/238.16)*E;
% E = E*100/max(abs(E));
end