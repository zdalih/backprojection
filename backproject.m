%MATH 406 - Project 2
%Author: Mohamed-Ali Hached (49356132)
%Date: December 10th 2017

clc
clear

%Import Z.dat data
Z = importdata("Z.dat");

%parameters
a = 10; %radius of the circle boundary
alpha = (0:1:7)*(2*pi)./16;
k = horzcat(alpha, alpha+pi, alpha(1)); %here alpha(1) is since it's a circle
N = 100; %NxN pixel division
pixelDim = 2*a/N; %width of each pixel

%the circle
U = zeros(N,N); %0 means no circle, 1 means yes circle
XY = cell(1,8); % a cell aray where the (X,+-Y) element of each pixel
                % in the circle is stored
m = [a a]; %vector from origin of suare to mid point of the circle

%get the nominal field measurements for a uniform conductivity
for i = 1 : 8
    for j = 1 : 16
        %uper uniform voltage
        uu = log((2*a^2*(1-cos(k(j+1)-alpha(i))))/(2*a^2*(1 -cos(k(j+1)-alpha(i)-pi))))/(2*pi);
        %lower uniform volage
        ul = log((2*a^2*(1-cos(k(j)-alpha(i))))/(2*a^2*(1 -cos(k(j)-alpha(i)-pi))))/(2*pi);
        %voltage difference 
        Zh(16*(i-1)+j,1) = uu-ul;
        %if Infinity make it 1E-09
        if abs(uu-ul) == inf
            Zh(16*(i-1)+j,1) = 1E-09;
        end
    end
end

%figure out if a pixel is in the circle or not
for r = 1 : N
    for c = 1 : N
        q = [(r-0.5)*pixelDim (c-0.5)*pixelDim]; %distance to mid of pixel in question
        %if the mid point of the pixel is in the circle make it a 0
        if norm(q - m) < a
            U(r,c) = 1;
            %For this point figure out the equipotential point on the
            %boundary
            point = q-m;
            x = point(1);
            y = point(2);
            X = (2*x*a^2)/(x^2+y^2+a^2);
            Y = sqrt(a^2-X^2);
            [theta,rho] = cart2pol(x,y); 
            %now for each pair of axis find the equipotential line
            for i = 1 : 8
                %equipotential boundary points for each axial node.
                XY{1,i}{r,c} = [theta-pi/2+(i-1)*(2*pi)/16,rho];
            end
        end
    end
end

%for each pixel find out which terminals it is stuck between 
%using the equipotential lines
sigma = U; %a grid of the conductivities

for i = 1 : 8
    for r = 1 : N
        for c = 1 : N
            if U(r,c) == 1
                
                TR = XY{1,i}{r,c};
                [x,y] = pol2cart(TR(1),TR(2)); 
                [theta,rho] = cart2pol(x,y); 
                n = abs(16*theta/(2*pi));
                
                if n - round(n,0) > 0
                    nl = round(n,0)+1;
                else
                    nl = round(n,0);
                end
                
                sigma(r,c) = sigma(r,c) - Z(16*(i-1)+nl)/Zh(16*(i-1)+nl) + 1;
                
            end
        end
    end
end

%display the region
image(sigma, 'CDataMapping','scaled');
colorbar
