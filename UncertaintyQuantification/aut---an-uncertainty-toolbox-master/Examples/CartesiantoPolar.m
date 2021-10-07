% This example consider the simple case of converting a point from Cartesian to Polar Coordinates.
% The resulting uncertainty are computed using Linearisation, Monte-Carlo(MC) and Unscented
% transformation(UT) methods.

%   unc=@unc_t; % uncomment this line for use MC method
   unc=@unc_ut; % uncomment this line for use UT method
 % coordinates in the cartesian plane(x,y)
x = unc(0.4,0.03,'x');
y = unc(0.3,0.01,'y');

% coordinates in the polar plane (r, ph)
r = sqrt(x^2 + y^2);
th = atan2(y,x);

% plot the mean and Covariance ellipse for all the methods in the same figure 
figure(1);
title('Cartesian coordinates(x,y)','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
if isa(x,'unc')
    
        plot(gmv(x),gmv(y),'or','LineWidth',2); 
        hold on;
        % plot covariance ellipse 2sigma
        CovI = get_cov_mat([x y]);
        [ellipse_x1Lin,ellipse_x2Lin] = Ellipse(CovI,gmv(x),gmv(y),0.95);
        lin = plot(ellipse_x1Lin,ellipse_x2Lin,'r','LineWidth',2);
        
elseif isa(x,'unc_t')
    
        plot(gmv(x),gmv(y),'ob','LineWidth',2); 
        hold on;
        %scatter(x.values,y.values,'b.'); % plot MC Samples
        %hold on;
        % plot covariance ellipse 2sigma
        CovI = get_cov_mat([x y]);
        [ellipse_x1MC,ellipse_x2MC] = Ellipse(CovI,gmv(x),gmv(y),0.95);
        MC = plot(ellipse_x1MC,ellipse_x2MC,'b','LineWidth',2);
        
elseif isa(x,'unc_ut')
    
        plot(gmv(x),gmv(y),'om','LineWidth',2); 
        hold on;
        % plot covariance ellipse 2sigma
        CovI = get_cov_mat([x y]);
        [ellipse_x1UT,ellipse_x2UT] = Ellipse(CovI,gmv(x),gmv(y),0.95);
        UT = plot(ellipse_x1UT,ellipse_x2UT,'m','LineWidth',2);
end 

figure(2);
title('Polar Coordinates(r,$phi$)','interpreter','latex');
xlabel('$phi$','interpreter','latex');
ylabel('r','interpreter','latex');
if isa(r,'unc')
    
        plot(gmv(th),gmv(r),'or','LineWidth',2); 
        hold on;
        % plot covariance ellipse 2sigma
        CovO = get_cov_mat([ th r]);
        [ellipse_x1Lin1,ellipse_x2Lin1] = Ellipse(CovO,gmv(th),gmv(r),0.95);
        lin1 = plot(ellipse_x1Lin1,ellipse_x2Lin1,'r','LineWidth',2);
        
elseif isa(r,'unc_t')
    
        plot(gmv(th),gmv(r),'ob','LineWidth',2); 
        hold on;
%         scatter(ph.values,r.values,'b.'); % plot MC Samples
%         hold on;
        % plot covariance ellipse 2sigma
        CovO = get_cov_mat([th r ]);
        [ellipse_x1MC1,ellipse_x2MC1] = Ellipse(CovO,gmv(th),gmv(r),0.95);
        MC1 = plot(ellipse_x1MC1,ellipse_x2MC1,'b','LineWidth',2);
        
elseif isa(r,'unc_ut')
    
        plot(gmv(th), gmv(r),'om','LineWidth',2); 
        hold on;
        % plot covariance ellipse 2sigma
        CovO = get_cov_mat([th r]);
        [ellipse_x1UT1,ellipse_x2UT1] = Ellipse(CovO,gmv(th),gmv(r),0.95);
        UT1 = plot(ellipse_x1UT1,ellipse_x2UT1,'m','LineWidth',2);
end 
% uncomment the next lines for plotting the legend after plot the data with the 3 Methods  
%   legend([lin, MC, UT],'Linearisation','MC ','UT');
%   legend([lin1, MC1, UT1],'Linearisation','MC ','UT');

