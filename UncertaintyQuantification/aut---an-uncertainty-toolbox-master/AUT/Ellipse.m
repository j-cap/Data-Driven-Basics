function [ellipse_x1,ellipse_x2] = Ellipse(Cov,Mean_X1,Mean_X2,Prob)
% Ellipse of the equal probability for 2D Gaussian Distribution    
%   This function generates the parameters for plotting the constant
%   probability ellipse considering a two dimensional Gaussian distribution. 
%   Cov --- Covariance matrix
%   Mean_X1 --- Mean of the variable X1
%   Mean_X2 --- Mean of the variable X2
%   Prob --- Constant probability at which we want to obtain the ellipse
%
%   ellipse_x1 --- x1 component for plotting the ellipse
%   ellipse_x2 --- x2 component for plotting the ellipse

% Calculate the eigenvectors and eigenvalues of the Covariance matrix 
[eigenvec, eigenval ] = eig(Cov);

% Get the index of the largest eigenvector. It gives the direction of the
% major axis of the ellipse 
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue. Scaling factor for the major axis
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue. Direction and scaling factor
% for the minor axis
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Defining the theta_grid for plotting the ellipse and the lenght of the major
% and minor axis. 
theta_grid = linspace(0,2*pi);
phi = angle;
a=sqrt(chi2inv(Prob,2))*sqrt(largest_eigenval); % major axis
b=sqrt(chi2inv(Prob,2))*sqrt(smallest_eigenval); % minor axis

% the ellipse in x and y coordinates 
ellipse_x1_r  = a*cos( theta_grid );
ellipse_x2_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [cos(phi) sin(phi); -sin(phi) cos(phi) ];
% R = eigenvec;

%let's rotate the ellipse to angle phi
r_ellipse = [ellipse_x1_r;ellipse_x2_r]' * R;

ellipse_x1= r_ellipse(:,1)+ Mean_X1;
ellipse_x2 = r_ellipse(:,2)+ Mean_X2;
end

