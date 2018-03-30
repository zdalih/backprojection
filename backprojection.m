%MATH 406 - Project 2
%Author: Mohamed-Ali Hached (49356132)
%Date: December 10th 2017
%
%Utilize voltage measurements at predetermined points around a circle to image the 
%center of the circle. Voltage measurements are in file Z.dat
%Using the Electrical Impedance Tomography Method of Back-projection
%
clc
clear



%Import Z.dat data
Z = importdata("Z.dat");

%parameters
a = 10; %radius of the circle boundary
N = 800;

%make the 'pixels' using meshgrid
drho = a/N;
dtheta = 2*pi/N;
rho = 0 : drho : a;
theta = 0 : dtheta : 2*pi;
[meshRho, meshTheta] = meshgrid(rho,theta); 

%position of the terminals and the axes
alpha = (0:1:7)*(2*pi)./16;
k = horzcat(alpha, alpha+pi, alpha(1)); %here alpha(1) since it's a circle

%now we're gonna rotate the axes for each terminal and find
%the point (X,Y) where the equipotential lines for
%the uniform field intersects the boundary of the circle
%using this X,Y we want to find between which terminals it lies
sigma = 0*meshRho + 1; %the conductivity 
for i = 0 : 1 : 7
    alphaShift = circshift(-alpha,i);
    ZShift = circshift(Z, 16*i);
    angle = alphaShift(1);
    %coordinates of the pixels in cartesian
    %rotated so that the x axis lies with the 
    %excitation axis
    x_nr = meshRho.*cos(meshTheta);
    y_nr = meshRho.*sin(meshTheta);
    %do the rotation counterclockwise by the angle
    x = x_nr.*cos(angle) - y_nr.*sin(angle);
    y = x_nr.*sin(angle) + y_nr.*cos(angle);
    
    %intersections with the boundary
    %this is again with the rotated axes
    X = 2*a^2*x./(x.^2+y.^2+a^2);
    Y = sqrt(a^2-X.^2);
    %polar coordinates for these intersections
    [THETA, RHO] = cart2pol(X,Y); 
    
    
    %get the uniform field vector, for 
    for j = 1 : 16
        %upper uniform voltage
        uu = log((2*a^2*(1+cos(k(j+1))))/(2*a^2*(1 -cos(k(j+1)))))/(2*pi);
        %lower uniform volage
        ul = log((2*a^2*(1+cos(k(j))))/(2*a^2*(1 -cos(k(j)))))/(2*pi);
        %voltage difference
        Zh(j,1) = uu-ul;
        %if Infinity make it 1E-09
        if abs(uu-ul) == inf
            Zh(j,1) = 1E-09;
        end
    end
    
    %the lowest integer to n should be the number of the lower terminal
    %that traps the pixel
    n1 = floor(16*(THETA)./(2*pi) + 1);
    n2 = (18-n1) - 1;
    
    %now find the change in conductivity for the pixels
    sigma = sigma + 2 - ZShift(n2)./Zh(n2) - ZShift(n1)./Zh(n1) ;
end


h = surf(x_nr,y_nr,sigma);
set(h, 'edgecolor', 'none');
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');
title('Electrical Impedance Tomography Method of Back-projection', 'interpreter', 'Latex');
