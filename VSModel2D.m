%% Mahdiar Sadeghi
%% Northeastern University
%% December 6th 2017
%
% Inputs : dielectric permittivity (eps), dielectric thickness(d) and
% array distance from metal plate (H)
%
% notes:
% -all distance values should be in cm
% -to have lossless dielectric, use loss=0
% -for bare metal case, use eps=1 and d=0
%
% Outputs : E response singal as a function z (distance from the detector)

function [E, z] = VSModel(loss, eps, d, H)

%% Global Variables
c = 3e10; f=24.16e9;
n = sqrt(eps);
k = 2*pi*f/c;
gam = (1-n)/(1+n); 
n1 = -1.1*n;

loss1 =  exp(-k*2*d*n*loss);
loss2 =  exp(-k*4*d*n*loss);
loss3 =  exp(-k*6*d*n*loss);
loss4 =  exp(-k*8*d*n*loss);

%% half width of the panel
hwidth = 50;

x = -hwidth:.25:hwidth;
[~,xsize] = size(x);

%% depth axis resolution per cm
res = 20;
Lowerbound = res*(-20);
Upperbound = res*(+20);

z = (Lowerbound:Upperbound)/res;
E = zeros(size(z));  

for count=Lowerbound:Upperbound
    del = count/20;
    
    comp1 = gam* sum(exp(-1i*k*(-sqrt(x.^2 + (H-del)^2 * ones(1,xsize))+ sqrt(x.^2 + (-2*d+H+del)^2 * ones(1,xsize)) )));
    comp2 = loss1*(-1+gam^2)*sum(exp(-1i*k*(2*d*(n+1/n1)*ones(1,xsize)  -sqrt(x.^2 + (H-del)^2 * ones(1,xsize))+ sqrt(x.^2 + (-2*d*(1+1/n1)+H+del)^2 * ones(1,xsize)) )));
    comp3 = loss2* gam * (-1+gam^2)*sum(exp(-1i*k*(4*d*(n+(1/n1))*ones(1,xsize)  -sqrt(x.^2 + (H-del)^2 * ones(1,xsize))+ sqrt(x.^2 + (-2*d*(1+2/n1)+H+del)^2 * ones(1,xsize)) )));
    comp4 = loss3* gam^2 * (-1+gam^2)*sum(exp(-1i*k*(6*d*(n+(1/n1))*ones(1,xsize)  -sqrt(x.^2 + (H-del)^2 * ones(1,xsize))+ sqrt(x.^2 + (-2*d*(1+3/n1)+H+del)^2 * ones(1,xsize)) )));
    comp5 = loss4* gam^3 * (-1+gam^2)*sum(exp(-1i*k*(8*d*(n+(1/n1))*ones(1,xsize)  -sqrt(x.^2 + (H-del)^2 * ones(1,xsize))+ sqrt(x.^2 + (-2*d*(1+4/n1)+H+del)^2 * ones(1,xsize)) )));
    
    E (count-Lowerbound+1) = comp1+comp2+comp3+comp4+comp5;
end

%% normalize the results
z = H-z;
E = (100/401)*E;
end