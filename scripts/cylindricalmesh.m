%%
clear all;

Nr = 200;
Ntheta = 160;

firstoff = 0.001;
R = 0.5;

growthrate = 1.055;
Nlayers = 50;

for i = 1:Nlayers
   if(i == 1)
      r(i) = R; 
   elseif(i == 2)
      r(i) = r(i-1) + firstoff; 
   else
      r(i) = r(i-1) + growthrate*(r(i-1)-r(i-2));
   end
end

outer_growth = 1.025;
for i = (Nlayers+1):Nr
    r(i) = r(i-1) + outer_growth*(r(i-1)-r(i-2));
end

theta = linspace(0, 2*pi, Ntheta+1);
theta = theta(1:(end-1));

for i = 1:Nr
    for j = 1:Ntheta
        ii = (i-1)*Ntheta + j;
        x(ii) = r(i)*cos(theta(j));
        y(ii) = r(i)*sin(theta(j));
    end
end