classdef unc_ut < handle
  % AUT - An Uncertainty Toolbox V 0.91
  % Date: 27/11/2019
  % Type "help disclaimer" for more information

   
  
  %% Properties of the unc_utertainty class
   
       properties  
       % Define the properties of the unc object to be created by the Class Constructor 
          nom_value ;  % Nominal value of the unc object
          dep;     % Array containing all the random variables (unc objects) 
                   % on which the unc object depends
          sigma;   % Array containing the sigma points of the unc object
       
       end
       
                   
  %% Methods of the uncertainty class    
       
       methods
      
  %% ------ Class constructor ------- %
  
       function obj = unc_ut(arg1,arg2,arg3)
                           
          global AUTTOOL;
          if isempty(AUTTOOL)
              AUTTOOL=1;
              help disclaimer;
              disp(['Press any key to continue...',13])
              pause;
          end
                           
                           
          global Number_of_Digits_to_Display;
          if isempty(Number_of_Digits_to_Display)
               Number_of_Digits_to_Display=2;
               fprintf(['Using %d digits to report uncertainty (global Number_of_Digits_to_Display)\n',13],Number_of_Digits_to_Display);
          end
          
          % Defining parameters for the Scaled Unscented Transform 
          global kappa;
          global beta;
          if isempty(kappa)&& isempty(beta) 
               kappa = 0; % kappa >=0
               beta = 2; % 
               fprintf(['Using scaling parameters \x03BA = %d and \x03B2 = %d.(Scaled Unscented Transform) \n',13],kappa, beta);
          end
          
          % Case(1): Number of input arguments is equal to 3:
          % Case(2): Number of input arguments is equal to 2:
               
          if (nargin>=2) 
              % Syntax (2):  unc_ut(Value, Uncertainty)
                 
                if isa(arg1,'double')&& isa(arg2, 'double') && ((size(arg1,1)==size(arg2,1) && size(arg1,2)==size(arg2,2))) %&& all(all(isreal(arg1))) 
                      obj(size(arg1,1),size(arg1,2))= unc_ut;
                      for i=1:size(arg1,1)
                          for j=1:size(arg1,2)
                              obj(i,j).nom_value = arg1(i,j);
                              obj(i,j).dep = obj(i,j);
                              obj(i,j).sigma = [arg1(i,j)- arg2(i,j);arg1(i,j)+ arg2(i,j)];
                                                        
                          end
                      end
                % when arg1 is a column vector and arg2 is a row vector or viceversa      
                elseif isa(arg1,'double')&& isa(arg2, 'double') && ((length(arg1)==length(arg2) && isvector(arg2))) %&& all(all(isreal(arg1))) 
                      obj(size(arg1,1),size(arg1,2))= unc_ut;
                      for i=1:length(arg1)
                          obj(i).nom_value = arg1(i);
                          obj(i).dep = obj(i);
                          obj(i).sigma = [arg1(i)- arg2(i);arg1(i)+ arg2(i)];
                         
                      end
                      
                  %%%%%%%%%%%%%%%% Syntax (3) and (4): (array of unc objects, input covariance matrix, 'cov') or 
                    %                     (array of unc objects, input correlation matrix, 'corr')
                 elseif (isa(arg1,'unc_ut')||isa(arg1,'double')) && isa(arg2, 'double') &&  size(arg2,1)==length(arg1) && size(arg2,2)==length(arg1)
                        
                        %check if the input matrix is symmetric and PSD  
                        [~,p] = chol(arg2);
                        if (isequal(arg2,transpose(arg2)) && (p==0)) 
        
                            if ((nargin==3)&& strcmp(arg3,'cov'))|| (nargin==2) % If user enters a covariance matrix 
                                    
                                if isa(arg1,'double')
                                      arg1 = unc_ut(arg1,sqrt(diag(arg2)).');
                                end
                                    obj = arg1;
  
%                                     perform the karhunen-loeve decomposition. update the dep property 
                                       [K, Coef_mat, eigD] = karhloev( arg1, arg2,'cov');
                                    
                                       % matrix square root
                                       S = sqrtm(arg2);
%                                      % get values of std_unc of the obj
%                                      to check if the values in the Cov
%                                      matrix are inside the range.
                                      [~,std_unc] = eval_std_unc(arg1);   
                                                                                                                   
                                      % %+++ Covariances and correlations assignments ++++
                                      n = length(arg1); % number of unc class objects
                                      for i=1:n
                                         for j=1:length(K)
                                             %check if the values of covariances are inside the range.
                                             
                                                   if (i~=j) && (( arg2(i,j) > (std_unc(i) * std_unc(j)) ) || ( arg2(i,j) < -(std_unc(i) * std_unc(j)) ))
                                                       error('Incorrect value in entry(%d,%d) of the covariance matrix! \n The covariance between two random variables x and y must satisfy the property |cov(x,y)|%s %sx*%sy, where %sx and %sy represent the standard deviation of x and y respectively.'...
                                                           ,i,j,num2str(char(8804)),num2str(char(963)),num2str(char(963)),num2str(char(963)),num2str(char(963)));
                                                   else
                                                      obj(i).dep(j) = K(j);
                                                      obj(i).sigma(1,j) = obj(i).nom_value - S(i,j);
                                                      obj(i).sigma(2,j) = obj(i).nom_value + S(i,j);
                                                   end
                                         end
                                      end 
                                      
                            elseif ((nargin==3)&& strcmp(arg3,'corr')) 
                                
                                [CX_cov] = cor_into_cov_mat(arg1,arg2);
                                
                                if isa(arg1,'double')
                                      arg1 = unc_ut(arg1,diag(CX_cov).');
                                end
                                    obj = arg1;
                                
%                               perform the karhunen-loeve decomposition. update the dep property 
                                       [K, Coef_mat, eigD] = karhloev( arg1, CX_cov,'cov');
                                    
                                       % matrix square root
                                       S = sqrtm(arg2);
%                                                                                                                                                        
                                      % %+++ Covariances and correlations assignments ++++
                                      n = length(arg1); % number of unc class objects
                                      for i=1:n
                                         for j=1:length(K)
                                             
                                                      obj(i).dep(j) = K(j);
                                                      obj(i).sigma(1,j) = obj(i).nom_value - S(i,j);
                                                      obj(i).sigma(2,j) = obj(i).nom_value + S(i,j);
                                         end
                                      end 
                                      

                                     % % End of assignments % %                                                                                     
   
                              else
                                   error('Incorrect type of argument: use cov to enter a covariance matrix, or corr to enter a correlation matrix.');
                              end
                                    
                        else
                            error('The covariance and correlation matrices must be symmetric and positive semidefinite.');
                        end 
                   
                                else
                                        error('Incorrect type of argument(s)');
                        
                end
                       
% % %                     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
          % Case(3): Number of input arguments is equal to 1
                
          elseif (nargin==1)    
              
             % Syntax(1): Value only
             
             if isa(arg1,'double') %&& all(all(isreal(arg1)))
                 
                       obj(size(arg1,1),size(arg1,2))=unc_ut;
                       for i=1:size(arg1,1)
                          for j=1:size(arg1,2)
                              obj(i,j).nom_value = arg1(i,j);
                              obj(i,j).dep=unc_ut.empty(1,0);
                              obj(i,j).sigma=double.empty(2,0);
                              
                          end
                       end
                       

             % Syntax(2): unc_ut object only
             elseif isa(arg1,'unc_ut')
                 
                   obj = arg1;
  
             end
             
          % Case(4): No input arguments
          elseif (nargin==0) 
          
                 obj.nom_value = [];
                 obj.dep=[];
                 obj.sigma=[];
                 
          else
              
              error('Incorrect type of argument(s)');
                          
          end
  
       end  
       
 %% Karhunen-Loeve decomposition %%%%%%%%%%%%%%%%%%%%     
       function [K, Coef_mat, eigD]=karhloev( X,CX,arg3 )
            
            % This method generates uncorrelated unc objects K1,K2,...,Kn 
            % from correlated unc objects x1,x2,...,xn using Karhunen-Loeve decomposition.
            % K : vector of uncorrelated unc_ut objects.
            % Coef_mat : matrix of eigenvectors 
            % eigD: is diagonal matrix of eigenvalues.
            % Input arguments:
            % X: array of unc objects [1*n] (correlated unc objects)
            % CX: Input covariance or correlation matrix
            % arg3: use 'cov' for covariance matrix and 'corr' for correlation matrix
        
        
            
        if (nargin==1) || (nargin==3)    
            
                arg=0;
                
                if nargin==1
                    
                    Z=X;
                    [CX]=get_cov_mat(X);
                    arg3='cov';
                    arg=1;
                elseif nargin==3
                     Z=X;
                end
                
                if iscolumn(Z)
                    Z = Z.';
                    
                end
                     

                if (strcmp(arg3,'cov')) && isa(CX,'double') ...
                        && ( ( size(CX,1)==length(Z) && size(CX,2)==length(Z) ) || (arg==1) ) && size(Z,1)==1 
                       
                    %%% check if CX is diagonal, then the objects are
                    %%% uncorrelated
                    if isdiag(CX) 
                        K = X';
                        Coef_mat = eye(length(X));
                        eigD = CX;
                    else

                   % Eigenvectors 

                    [ E , eigD] = eig(CX);
                    
                    %check if eigenvalues = 0(consider eig < eps as 0) . 
                    eig0 = find(diag(eigD) < eps);
                    %remove the eigenvectors corresponding to eigenvalues 0
                    E(:,eig0)= [];
                    eigD(:,eig0)= [];

                    % Generating the uncorrelated unc objects K1,K2,...Kn                    
                    K = E' * Z';

                    % Coefficients matrix
                    Coef_mat = E;
                                                     
                    
                    % updating dep property for uncorrelated variables
                       for k = 1:length(K)
                           K(k).sigma = [K(k).nom_value - sqrt(eigD(k,k)); K(k).nom_value + sqrt(eigD(k,k))];
                           K(k).dep = K(k);                        
                       end
                    % updating dep property for initial variables
% % % % % % % % % % % % Subrouting: if Coef_mat(i,j)=0, K(j) it is not included in the dependency list of X(i).      
% %                        for i=1:length(X) % index the original variables                        
% %                            for j = 1:size(Coef_mat,2) % index the uncorrelated variables
% %                                if Coef_mat(i,j)~=0
% %                                   X(i).dep(j) = K(j);
% %                                   X(i).sigma(j) = [K(j)- Coef_mat(j,j); K(j)+ Coef_mat(j,j)];                              
% %                                end 
% %                            end
% % %                            s = eval_std_unc(X(i));
% % %                            X(i).std_unc = s.std_unc;
% %                        end
% % % % % % % % % % % % % % % % % % % % 
                    end 
                else
                    error (' Incorrect input arguments!');

                end
                
        else
            
            error('Incorrect number of input arguments!');
        end

        
       end

     %% Converting correlation into Covariance matrix
     
     function [CX_cov]=cor_into_cov_mat(X,CX_cor)
          
          % This method transforms a correlation matrix into a covariance matrix 
          %
          % X: array of unc objects
          % CX_cor: correlation matrix of X
          % CX_cov: covariance matrix of X
          %
          [~,p] = chol(CX_cor);
          n=size(CX_cor,2); 
          if (n~= length(X))
               error('Dimension between the input parameters must agree');
               
          %%check if the correlation matrix CX_cor is symmetric and PSD 
          elseif (isequal(CX_cor,transpose(CX_cor)) && (p==0))
              
              CX_cov=zeros(n,n);
              for i=1:n
                    for j=1:n
                       %check if the values of covariances
                       %are inside the range. 
                       if (i~=j) && ((CX_cor(i,j) > 1 ) || ( CX_cor(i,j) < -1 ))
                         error('Incorrect value in entry(%d,%d) of the correlation matrix! \n The correlation coefficient between two random variables x and y must satisfy the property |corr(x,y)|%s 1.'...
                         ,i,j,num2str(char(8804)));
                     
                       elseif (i==j) && (CX_cor(i,j) ~= 1 )
                          error('The main diagonal elements of a correlation matrix must be equal to 1.'); 
                       end
                       % get values of std_unc of the obj to obtain the cov
                       % matrix
                       [~,std_unc] = eval_std_unc(arg1);
                       CX_cov(i,j)= CX_cor(i,j) * std_unc(i) * std_unc(j);
                       
                    end 
                                                      
              end
              
          else
              error('Correlation matrix must be symmetric and positive semidefinite');
          end 
           
          
       end

%% %%%%%%%%% Basic Operations %%%%%%%%%
           
       function obj = plusscal(obj1,obj2) 
          
       % Method used by operator overloading methods
              
            obj1= unc_ut((obj1));
            obj2= unc_ut((obj2));
 
            obj = unc_ut;
            % Updating nominal value property
            obj.nom_value = obj1.nom_value+obj2.nom_value; 
            
            % Updating dep property        
            [~,ia1,ib1]=intersect(obj1.dep,obj2.dep); % in both
            [~,ia2,ib2]=setxor(obj1.dep,obj2.dep);  % only in one
            
             obj.dep=[obj1.dep(ia1),obj1.dep(ia2),obj2.dep(ib2)];         
                       
             obj.sigma=zeros(2,length(obj.dep));  
             
             % Defining the sigma points 
             lboth=length(ia1);
             la=length(ia2);
             lb=length(ib2);
             j=0;
             for i=1:lboth
                 j=j+1;
                 obj.sigma(1,i)=obj1.sigma(1,ia1(i))+obj2.sigma(1,ib1(i));
                 obj.sigma(2,i)=obj1.sigma(2,ia1(i))+obj2.sigma(2,ib1(i));
             end
                         
             
             for i=1:la
                 j=j+1;
                 obj.sigma(1,j)=obj1.sigma(1,ia2(i))+obj2.nom_value;
                 obj.sigma(2,j)=obj1.sigma(2,ia2(i))+obj2.nom_value;
             end
                          
             for i=1:lb
                 j=j+1;
                 obj.sigma(1,j)=obj1.nom_value+obj2.sigma(1,ib2(i));
                 obj.sigma(2,j)=obj1.nom_value+obj2.sigma(2,ib2(i));
             end

       end
           
            
  %%
  
   function obj = minusscal(obj1,obj2) 
   
   % Method used by operator overloading methods    
       
            obj1= unc_ut((obj1));
            obj2= unc_ut((obj2));
            obj = unc_ut;
            % Updating value property
            obj.nom_value = obj1.nom_value-obj2.nom_value; 
            
            % Updating dep property        
            [~,ia1,ib1]=intersect(obj1.dep,obj2.dep); % in both
            [~,ia2,ib2]=setxor(obj1.dep,obj2.dep);  % only in one
            
            obj.dep=[obj1.dep(ia1),obj1.dep(ia2),obj2.dep(ib2)];         
            
            % Defining the sigma points
             obj.sigma=zeros(2,length(obj.dep));  
      
             lboth=length(ia1);
             la=length(ia2);
             lb=length(ib2);
             j=0;
             for i=1:lboth
                 j=j+1;
                 obj.sigma(1,i)=obj1.sigma(1,ia1(i))-obj2.sigma(1,ib1(i));
                 obj.sigma(2,i)=obj1.sigma(2,ia1(i))-obj2.sigma(2,ib1(i));
             end
                          
             
             for i=1:la
                 j=j+1;
                 obj.sigma(1,j)=obj1.sigma(1,ia2(i))-obj2.nom_value;
                 obj.sigma(2,j)=obj1.sigma(2,ia2(i))-obj2.nom_value;
             end
                          
             for i=1:lb
                 j=j+1;
                 obj.sigma(1,j)=obj1.nom_value-obj2.sigma(1,ib2(i));
                 obj.sigma(2,j)=obj1.nom_value-obj2.sigma(2,ib2(i));
             end

       end         
        
  %%
  
  function obj = multscal(obj1,obj2)
           
       % Method used by operators overloading methods
       
            obj1= unc_ut((obj1));
            obj2= unc_ut((obj2));
            obj = unc_ut;
            % Updating value property
            obj.nom_value = obj1.nom_value*obj2.nom_value ; 
            
            % Updating dep property        
            [~,ia1,ib1]=intersect(obj1.dep,obj2.dep); % in both
            [~,ia2,ib2]=setxor(obj1.dep,obj2.dep);  % only in one
           
            obj.dep=[obj1.dep(ia1),obj1.dep(ia2),obj2.dep(ib2)];         
            
            % Defining Sigma Points
             obj.sigma=zeros(2,length(obj.dep));  
    
             lboth=length(ia1);
             la=length(ia2);
             lb=length(ib2);
             j=0;
             for i=1:lboth
                 j=j+1;
                 obj.sigma(1,i)=obj1.sigma(1,ia1(i))*obj2.sigma(1,ib1(i));
                 obj.sigma(2,i)=obj1.sigma(2,ia1(i))*obj2.sigma(2,ib1(i));
             end
                         
             
             for i=1:la
                 j=j+1;
                 obj.sigma(1,j)=obj1.sigma(1,ia2(i))*obj2.nom_value;
                 obj.sigma(2,j)=obj1.sigma(2,ia2(i))*obj2.nom_value;
             end
                          
             for i=1:lb
                 j=j+1;
                 obj.sigma(1,j)=obj1.nom_value*obj2.sigma(1,ib2(i));
                 obj.sigma(2,j)=obj1.nom_value*obj2.sigma(2,ib2(i));
             end

  end
  
  
  %%     
 
       function obj = divscal(obj1,obj2)
           
           % Method used by operator overloading methods
           
            obj1= unc_ut((obj1));
            obj2= unc_ut((obj2));
            obj = unc_ut;
            % Updating value property
            obj.nom_value = obj1.nom_value/obj2.nom_value ; 
            
            % Updating dep property        
            [~,ia1,ib1]=intersect(obj1.dep,obj2.dep); % in both
            [~,ia2,ib2]=setxor(obj1.dep,obj2.dep);  % only in one
            
            
            obj.dep=[obj1.dep(ia1),obj1.dep(ia2),obj2.dep(ib2)];
            
            % Defining Sigma Points            
            obj.sigma=zeros(2,length(obj.dep));  
          
             lboth=length(ia1);
             la=length(ia2);
             lb=length(ib2);
             j=0;
             for i=1:lboth
                 j=j+1;
                 obj.sigma(1,i)=obj1.sigma(1,ia1(i))/obj2.sigma(1,ib1(i));
                 obj.sigma(2,i)=obj1.sigma(2,ia1(i))/obj2.sigma(2,ib1(i));
             end
                          
             
             for i=1:la
                 j=j+1;
                 obj.sigma(1,j)=obj1.sigma(1,ia2(i))/obj2.nom_value;
                 obj.sigma(2,j)=obj1.sigma(2,ia2(i))/obj2.nom_value;
             end
                          
             for i=1:lb
                 j=j+1;
                 obj.sigma(1,j)=obj1.nom_value/obj2.sigma(1,ib2(i));
                 obj.sigma(2,j)=obj1.nom_value/obj2.sigma(2,ib2(i));
             end
               
       end
           
           
                      
         %%
         
          function obj = plus(obj1,obj2)
              
              % Overloading of the Addition Operator
              obj1= unc_ut((obj1));
              obj2= unc_ut((obj2));
              
              %check the dim
              %%% one of the objects is scalar 
              if length(obj1)==1
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i = 1:size(obj2,1)
                      for j = 1:size(obj2,2)
                          obj(i,j)= plusscal(obj1,obj2(i,j)); 
                      end
                  end
              elseif length(obj2)==1
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  
                  for i = 1:size(obj1,1)
                      for j = 1:size(obj1,2)
                          obj(i,j)= plusscal(obj1(i,j),obj2); 
                      end
                  end
                   %%% vector-vector/ matrix-matrix addition. Objects have the same dim   
              elseif size(obj1)==size(obj2)
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= plusscal(obj1(i,j),obj2(i,j));                        
                      end                    
                  end
                  %%% vector-matrix addition 
              elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= plusscal(obj1(i),obj2(i,j));                        
                      end                    
                  end
              elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= plusscal(obj1(j),obj2(i,j));                        
                      end                    
                  end
              elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= plusscal(obj2(i),obj1(i,j));                        
                      end                    
                  end
              elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= plusscal(obj2(j),obj1(i,j));                        
                      end                    
                  end
                  %%% vector-vector addition
               elseif (isrow(obj1) && iscolumn(obj2))
                  obj(length(obj2),length(obj1))=unc_ut;
                  for i=1:length(obj2)
                      for j=1:length(obj1)
                          obj(i,j)= plusscal(obj2(i),obj1(j));                        
                      end                    
                  end
               elseif (isrow(obj2) && iscolumn(obj1))
                  obj(length(obj1),length(obj2))=unc_ut;
                  for i=1:length(obj1)
                      for j=1:length(obj2)
                          obj(i,j)= plusscal(obj1(i),obj2(j));                        
                      end                    
                  end
              else
                  error('Invalid operation. Matrix dimension must agree')
              end
          end
          

            %% 
            
           function obj= minus(obj1,obj2)
               
           % Overloading of the Subtraction Operator
           
           %check the dim
              if length(obj1)==1
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i = 1:size(obj2,1)
                      for j = 1:size(obj2,2)
                          obj(i,j)= minusscal(obj1,obj2(i,j)); 
                      end
                  end
              elseif length(obj2)==1
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  
                  for i = 1:size(obj1,1)
                      for j = 1:size(obj1,2)
                          obj(i,j)= minusscal(obj1(i,j),obj2); 
                      end
                  end
                      
              elseif size(obj1)==size(obj2)
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= minusscal(obj1(i,j),obj2(i,j));                        
                      end                    
                  end
                %%% vector-matrix substraction 
              elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= minusscal(obj1(i),obj2(i,j));                        
                      end                    
                  end
              elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= minusscal(obj1(j),obj2(i,j));                        
                      end                    
                  end
              elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= minusscal(obj1(i,j),obj2(i));                        
                      end                    
                  end
              elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= minusscal(obj1(i,j),obj2(j));                        
                      end                    
                  end
                  %%% vector-vector substraction
               elseif (isrow(obj1) && iscolumn(obj2))
                  obj(length(obj2),length(obj1))=unc_ut;
                  for i=1:length(obj2)
                      for j=1:length(obj1)
                          obj(i,j)= minusscal(obj1(j),obj2(i));                        
                      end                    
                  end
               elseif (isrow(obj2) && iscolumn(obj1))
                  obj(length(obj1),length(obj2))=unc_ut;
                  for i=1:length(obj1)
                      for j=1:length(obj2)
                          obj(i,j)= minusscal(obj1(i),obj2(j));                        
                      end                    
                  end
              else
                  error('Invalid operation. Matrix dimension must agree')
              end 
               
           end
           
           
                   %% 
          
          function obj=times(obj1,obj2)
              
%               This method implements the overloading of the element-wise 
%               multiplication operator called for the syntax ' .*'
                           
              
              % check if the objects are complex and they are not unc obj
               
              
              issc=false; % issc: stands for 'is scalar'
              
              % if there is a scalar, then it should be obj1
              if isscalar(obj2) 
                  obh=obj1;
                  obj1=obj2;
                  obj2=obh;
              end
              
              if isscalar(obj1) 
                  issc=true;
              end              
              
%               
              %%%case with both elements have the same size or one of them is
              %%%an scalar 
              if min(size(obj1)==size(obj2)) || issc
                  obj = unc_ut(zeros(size(obj1)));
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                      if issc==1
                          obj(i,j)= multscal(obj1,obj2(i,j));
                      else
                          obj(i,j)= multscal(obj1(i,j),obj2(i,j));
                      end
                      end
                  end
                  
              %%%case when both element are vector     
              elseif (isvector(obj1)&& isvector(obj2)) && (iscolumn(obj1)&& isrow(obj2))
                  
                     obj = unc_ut(zeros(length(obj1),length(obj2)));
                      for i=1:length(obj1)
                        for j=1:length(obj2)
                           obj(i,j)= multscal(obj1(i),obj2(j));
                        end
                      end
              elseif (isvector(obj1)&& isvector(obj2)) && (iscolumn(obj2)&& isrow(obj1)) 
                     obj = unc_ut(zeros(length(obj2),length(obj1)));
                      for i=1:length(obj2)
                        for j=1:length(obj1)
                           obj(i,j)= multscal(obj2(i),obj1(j));
                        end
                      end
                      
              %%% case matrix-vector multiplication 
              elseif (isvector(obj1)&& ismatrix(obj2)) && (iscolumn(obj1) && length(obj1)==size(obj2,1))
                  obj = unc_ut(zeros(size(obj2))); 
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= multscal(obj1(i),obj2(i,j));
                      end 
                  end 
              elseif (isvector(obj1)&& ismatrix(obj2)) && (isrow(obj1) && length(obj1)==size(obj2,2))
                  obj = unc_ut(zeros(size(obj2))); 
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= multscal(obj1(j),obj2(i,j));
                      end
                  end
                 
              elseif (isvector(obj2)&& ismatrix(obj1)) && (iscolumn(obj2) && length(obj2)==size(obj1,1))
                  obj = unc_ut(zeros(size(obj1)));
                      for i=1:size(obj1,1)
                          for j=1:size(obj1,2)
                              obj(i,j)= multscal(obj2(i),obj1(i,j));
                          end
                      end 
              elseif (isvector(obj2)&& ismatrix(obj1)) && (isrow(obj2) && length(obj2)==size(obj1,2))
                  obj = unc_ut(zeros(size(obj1)));
                      for i=1:size(obj1,1)
                          for j=1:size(obj1,2)
                                obj(i,j)= multscal(obj2(j),obj1(i,j));
                          end
                      end
                                 
              else
                    error('Matrix dimension must agree')
              end
          end                       
              
                     %% 
          
          function obj =mldivide(obj1,obj2)
              
%               This method implements the overloading of matrix left division  
%               for two arrays obj1 and obj2 of uncertainty objects.
               
               %two objects are scalar
              if isscalar(obj1) && isscalar(obj2) 
                  obj = divscal(obj2,obj1);
                  
               %obj1 is scalar and obj2 is an array
              elseif isscalar(obj1) && (ismatrix(obj2)|| isvector(obj2))
                  obj = unc_ut(zeros(size(obj2)));
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= divscal(obj2(i,j),obj1);                    
                      end 
                  end
                  
                %obj1 is an square matrix
              elseif ismatrix(obj1) && (size(obj1,1)== size(obj1,2)) && (size(obj1,1)== size(obj2,1))
                 obj3=inv(obj1);
                 obj = mtimes(obj3,obj2);
                 
                %obj1 is an mxn matrix
              elseif ismatrix(obj1) &&(size(obj1,1)== size(obj2,1))
                  obj = inv(obj1'*obj1)*obj1'*obj2;
                  
              else 
                  error('Matrix dimension must agree');
                  
              end 
          end
          
          
       %% 
          
       function obj =mrdivide(obj1,obj2)
           
           % Overloading of the division operator
            
               %two objects are scalar
               if isscalar(obj1) && isscalar(obj2) 
                      
                  obj = divscal(obj1,obj2);
                  
               %obj2 is scalar and obj1 is an array
               elseif isscalar(obj2) && (ismatrix(obj1)|| isvector(obj1))
                  obj = unc_ut(zeros(size(obj1)));
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= divscal(obj1(i,j),obj2);                    
                      end 
                  end
                  
                %obj2 is an square matrix
                elseif ismatrix(obj2) && (size(obj2,1)== size(obj2,2)) && (size(obj2,2)== size(obj1,2))
                 obj3=inv(obj2);
                 obj = mtimes(obj1,obj3);
                 
                %obj2 is an mxn matrix
                elseif ismatrix(obj2) &&(size(obj2,2)== size(obj1,2))
                  obj = obj1*obj2'*inv(obj2*obj2');
                                                   
                else 
                  error('Matrix dimension must agree');
                  
                end 
                                              
     end
              
  
           
                     %% 
          
          function obj=rdivide(obj1,obj2)
              
%               This method implements the overloading of the right element-wise 
%               division operator (./) for two arrays obj1 and obj2 of uncertainty objects.

               if size(obj1)==size(obj2) 
                   obj=unc_ut(zeros(size(obj1)));
                   for i=1:size(obj1,1)
                       for j=1:size(obj1,2)
                           obj(i,j)= divscal(obj1(i,j),obj2(i,j));
                       end
                   end
               elseif length(obj2(:))==1
                   obj=unc_ut(zeros(size(obj1)));
                   for i=1:size(obj1,1)
                       for j=1:size(obj1,2)
                           obj(i,j)= divscal(obj1(i,j),obj2);
                       end
                   end
               elseif length(obj1(:))==1
                   obj=unc_ut(zeros(size(obj2)));
                   for i=1:size(obj2,1)
                       for j=1:size(obj2,2)
                           obj(i,j)= divscal(obj1,obj2(i,j));
                       end
                   end
             %%% vector-matrix element wise division 
              elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                  obj=unc_ut(zeros(size(obj2)));
                   for i=1:size(obj2,1)
                       for j=1:size(obj2,2)
                           obj(i,j)= divscal(obj1(i),obj2(i,j));
                       end
                   end
                  
              elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                  obj=unc_ut(zeros(size(obj2)));
                   for i=1:size(obj2,1)
                       for j=1:size(obj2,2)
                           obj(i,j)= divscal(obj1(j),obj2(i,j));
                       end
                   end
              elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                  obj=unc_ut(zeros(size(obj1)));
                   for i=1:size(obj1,1)
                       for j=1:size(obj1,2)
                           obj(i,j)= divscal(obj1(i,j),obj2(i));
                       end
                   end
              elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                  obj=unc_ut(zeros(size(obj1)));
                   for i=1:size(obj1,1)
                       for j=1:size(obj1,2)
                           obj(i,j)= divscal(obj1(i,j),obj2(j));
                       end
                   end
                  %%% vector-vector division
              elseif (isrow(obj1) && iscolumn(obj2))
                  obj=unc_ut(zeros(length(obj2),length(obj1)));
                   for i=1:length(obj2)
                       for j=1:length(obj1)
                           obj(i,j)= divscal(obj1(j),obj2(i));
                       end
                   end
              elseif (isrow(obj2) && iscolumn(obj1))
                  obj=unc_ut(zeros(length(obj1),length(obj2)));
                   for i=1:length(obj1)
                       for j=1:length(obj2)
                           obj(i,j)= divscal(obj1(i),obj2(j));
                       end
                   end
                   
              else
                   error('Arrays must have compatible size')
              end
          end
          
          
          %% 
          
          function obj=ldivide(obj1,obj2)
              
%               This method implements the overloading of the left 
%               element-wise division operator (.\) for two arrays obj1 and obj2 of uncertainty objects.

              obj = rdivide(obj2,obj1);
          end  
          
          
          
             
                    %%
         
           function obj = mtimes(obj1,obj2)

               % ------------------------------------------------------------ %
               %  mtimes(obj1,obj2): Overloading of matrix multiplication operator
               %  for uncertainty objects obj1 and obj2
               % ------------------------------------------------------------ %


                  obj1=unc_ut(obj1);
                  obj2=unc_ut(obj2);

                  if (isscalar(obj1))||(isscalar(obj2))
                      obj = obj1.*obj2;

                  elseif size(obj1,2)==size(obj2,1)  % X*Y: the number of columns of X must
                                                     %      equal the number of rows of Y
                      obj=unc_ut(zeros(size(obj1,1),size(obj2,2)));
                      for i=1:size(obj1,1)                 
                          for j=1:size(obj2,2)
                              for k=1:size(obj2,1)
                                 temp=multscal(obj1(i,k),obj2(k,j));
                                 obj(i,j)=plusscal(obj(i,j),temp);  

                              end 
                          end

                      end
                  else
                       error('Invalid operation. Matrix dimension must agree')                 
                  end
           end
       
           
         %% 
         
          function obj = prod(obj1)
              
               % obj=prod(obj1): Product of the elements of the object obj1. If obj1 is a
               % matrix of uncertainty objects, obj is a row vector with the product over each column. 
             
             if isvector(obj1) %obj1 is a vector
                 obj = obj1(1);
                 for i= 2:length(obj1)
                      obj = obj * obj1(i);
                 end
                
             elseif ismatrix(obj1)
                 obj(1,size(obj1,2))=unc_ut;
                 for j=1: size(obj1,2) % obj1 is a matrix
                     obj3 = obj1(1,j);                                       
                     for i = 2:size(obj1,1)
                         obj3 = obj3 * obj1(i,j);
                     end
                     obj(j) = obj3;
                 end
             end 
             
          end
         
          
          %% 
                    
        function obj = inv(obj1)
              
              %matrix inversion is possible only if:
             %   -) it is a square matrix
             %   -) its determinant is not zero.
             

                  % Find dimensions of input matrix
                    [r1,c1] = size(obj1);
                  
                    % If input matrix is not square, stop function
                    if r1 ~= c1
                       error('Error using inv. Matrix must be Square')                     
                    end
                    if det(gmv(obj1))== 0
                        error ('Matrix is singular.')
                    end 
                    if rank(gmv(obj1))< c1
                       warning('Matrix is close to singular or badly scaled.')
                    end

                    % Target identity matrix to be transformed into the output 
                    % inverse matrix
                   
                    obj=obj1.*0+eye(r1);

                    % matrix inversion
                    for j = 1 : r1 % r.value the nbr of rows of matrix
                        for i = j : r1
                            if gmv(obj1(i,j))~= 0
                                for k = 1 : r1
                                     % rows permutation
                                    s = obj1(j,k); obj1(j,k) = obj1(i,k); obj1(i,k) = s;
                                    s = obj(j,k); obj(j,k) = obj(i,k); obj(i,k) = s;
                                end
                                
%                                 if imag(obj1(j,j))~=0
%                                             t=1/obj1(j,j);
%                                 else
                                            t=divscal(1,obj1(j,j));
                                            
%                                 end
                                
                                
                                for k = 1 : r1
                                    obj1(j,k) = t * obj1(j,k);
                                    obj(j,k) = t * obj(j,k);
                                end
                                for L = 1 : r1
                                    if L ~= j
                                        t = -1*obj1(L,j);
                                        for k = 1 : r1
                                           
                                            obj1(L,k) = obj1(L,k) + t * obj1(j,k);
                                            obj(L,k) = obj(L,k) + t * obj(j,k);
                                        end
                                    end
                                end            
                            end
                            break
                        end
                    end
                      %
              
        end
          
        
        
        %% 
        
        function obj = diag(obj1)
            
            % diag(obj1): Return the diagonal elements of obj1
            
          [r1,c1] = size(obj1);
          obj = unc_ut;
                  
           % If input matrix is not square, stop function
           if r1 ~= c1
                  error('Matrix must be Square.')                     
           end  
           
           for i=1:r1
               obj(i) = obj1(i,i);
                 
           end
              
        end

        
         %% 
        
        function obj = trace(obj1)
            
            % trace(obj1): returns the sum of the diagonal elements of the square matrix obj1 of
            % uncertainty objects
            
          [r1,c1] = size(obj1);
                  
           % If input matrix is not square, stop function
           if r1 ~= c1
                  error('Matrix must be Square.')                     
           end  
              
           obj = sum(diag(obj1));
           
        end
 
        
        %% Overloading the ctranspose() function
        
        function obj = ctranspose(obj1)
            
          % This method overloads the complex conjugate transpose function.
          
          [r1,c1] = size(obj1);
          
          for i=1:r1
              for j=1:c1
                   
                  obj2(j,i)=obj1(i,j);  
                      
              end
          end  
           
          obj = conj(obj2);  
          
           
        end       
        
        

        %% 
        
        function [determ]=det(obj)
            
            % det(obj): Determinant of the uncertainty object obj as the sum of the determinants 
            % of 2x2 matrices. Considering that
            % det(obj) = sum(obj(i,j)*((-1)^(i+j))*Mij), for i=1,n and j=1 (we
            % compute the det throught the first column). 
            % obj(i,j) represents the element of the row i and the column j and 
            % Mij is the minor of the element in row i and column j(determinant of the 
            % reduced matrix obtained removing the row i
            % and the column j of the matrix obj.
            j = 1; % 1st column 
            [m,n]=size(obj);
            determ = 0;
            if(m~=n)
                disp('Matrix needs to be square');
            elseif(n==1)
                determ = obj(1,1);
            elseif(n==2)
                determ = obj(1,1)*obj(2,2) - obj(1,2)*obj(2,1);
            else
                M = obj ;
                for i=1:m
                    Mi = M;
                    Mi(:,j) = []; % matrix resulting of eliminating the column j
                    Mi(i,:) = []; % matrix resulting of eliminating the row i
                    M_ij = Mi;
                    detM_ij = det(M_ij);
                    determ = determ + M(i,j)*((-1)^(i+j))*detM_ij;
                    M = obj;            
                
                end
            end
        end
          
 
          %% 
          
          function obj=uplus(obj1)
              % Overloading of the unary plus operator
              obj=times(1,obj1);
          end
          
          
           %%            
           
          function obj=uminus(obj1)
              
              % Overloading of the unary minus operator
              
              obj=times(-1,obj1);
          end
          
          
          %% 
          
          function value=gmv(obj)
              
              % gmv(obj): Get the matrix of value properties of the
              % array obj of uncertainty objects 
              
              [value,~] = eval_std_unc(obj);
              
          end
          

          %% 
          
          function std_unc=gmu(obj)
              
              % gmu(obj): Get the matrix of std_unc proprties of the
              % array obj of uncertainty objects
              
              [~,std_unc] = eval_std_unc(obj); 
              
          end

          
          %%

           function [Value,Std_unc] = eval_std_unc(obj1)   
               
               % eval_std_unc(obj1): Evaluate the standard uncertainty of the uncertainty
               % object obj1
               
               Value = zeros(size(obj1));
               Std_unc = zeros(size(obj1));
               global kappa beta 
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                        
                        if ~isempty(obj1(i,j).dep)
                        % scaled UT according to http://ais.informatik.uni-freiburg.de/teaching/ws13/mapping/pdf/slam06-ukf-4.pdf
                        % UT parameters k, alpha influence how far the sigma points are away from the mean
                        % alpha, lambda, beta involved in computing the weight  
                        n=size(obj1(i,j).dep,2);  % dimensionality
                        alpha = 1/sqrt(n+kappa);      
                        lambda=alpha^2*(n+kappa)-n;% scaling parameter
                        W0=lambda/(n+lambda);  % weight sigma point 0
                        Wi=1/2/(n+lambda); 
                        value_transformed=sum(obj1(i,j).sigma(:))*Wi+obj1(i,j).nom_value*W0;
                        Pyy=Wi*sum((obj1(i,j).sigma(:)-value_transformed).^2)+(W0+1-alpha^2+beta)*(obj1(i,j).nom_value-value_transformed)^2;
                         %Pyy=Wi*sum((obj1(i,j).sigma(:)-obj1(i,j).nom_value).^2)+(W0+1-alpha^2+beta)*(value_transformed-obj1(i,j).nom_value)^2;
                                      
                        Std_unc(i,j) = sqrt(Pyy);
                        Value(i,j) = value_transformed;
                        
                        else
                            Std_unc(i,j) = 0;
                            Value(i,j) = obj1.nom_value;
                        end    
                        
                        % ------------------------------------------------ 
                      end
                  end       
           end   


%% %%%%%%%%%%% Overloading general Functions %%%%%%%%%%%
          %% sinus
          function obj=sin(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = sin(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = sin(obj1(i,j).sigma);
                  end 
              end
                            
          end
          
          function obj=sind(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = sind(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = sind(obj1(i,j).sigma);
                  end 
              end
                            
          end
                  
          %% cosine
          function obj=cos(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = cos(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = cos(obj1(i,j).sigma);
                  end 
              end
          end
          
          function obj=cosd(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = cosd(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = cosd(obj1(i,j).sigma);
                  end 
              end
          end
          
          %% tan
          function obj=tan(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = tan(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = tan(obj1(i,j).sigma);
                  end 
              end
          end
          
          function obj=tand(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = tand(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = tand(obj1(i,j).sigma);
                  end 
              end
          end
                               
          
          %% Inverse cosine: acos
          function obj=acos(obj1)
              
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = acos(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = acos(obj1(i,j).sigma);
                  end 
              end
          end
          
          function obj=acosd(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = acosd(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = acosd(obj1(i,j).sigma);
                  end 
              end
          end
          
          %% Inverse sine: asin
          function obj=asin(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = asin(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = asin(obj1(i,j).sigma);
                  end 
              end
          end
          
          function obj=asind(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = asind(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = asind(obj1(i,j).sigma);
                  end 
              end
          end
          
          %% Inverse tangent: atan
          function obj=atan(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = atan(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = atan(obj1(i,j).sigma);
                  end 
              end
          end
          
          
          function obj= atand(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = atand(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = atand(obj1(i,j).sigma);
                  end 
              end
          end
          
          %% Hyperbolic cosine: cosh
          function obj = cosh(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = cosh(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = cosh(obj1(i,j).sigma);
                  end 
              end
          end
          
          %% Hyperbolic sine: sinh
          function obj = sinh(obj1)
              
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = sinh(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = sinh(obj1(i,j).sigma);
                  end 
              end
              
          end
          
          %% Hyperbolic tangent:tanh
          function obj = tanh(obj1)
              
             obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = tanh(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = tanh(obj1(i,j).sigma);
                  end 
              end
 
          end
          
          %% Hyperbolic cosine inverse: acosh
          function obj = acosh(obj1)
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = acosh(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = acosh(obj1(i,j).sigma);
                  end 
              end

          end
          
          %% Hyperbolic sine inverse: asinh
          function obj = asinh(obj1)
                     
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = asinh(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = asinh(obj1(i,j).sigma);
                  end 
              end

          end
          
          %% Hyperbolic tangent inverse: atanh
          function obj = atanh(obj1)
             
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = atanh(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = atanh(obj1(i,j).sigma);
                  end 
              end
          end
          
          %% Cotangent : cot
          function obj = cot(obj1)
              
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = cot(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = cot(obj1(i,j).sigma);
                  end 
              end

          end
          
          function obj = cotd(obj1)
              
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = cotd(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = cotd(obj1(i,j).sigma);
                  end 
              end

          end
          
          %% Inverse cotangent: acot
          function obj = acot(obj1)
             
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = acot(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = acot(obj1(i,j).sigma);
                  end 
              end

          end
          
          function obj = acotd(obj1)
             
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = acotd(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = acotd(obj1(i,j).sigma);
                  end 
              end

          end
          
          %% Hyperbolic cotangent: coth
          function obj = coth(obj1)
              
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = coth(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = coth(obj1(i,j).sigma);
                  end 
              end

          end
          
          %% Natural logarithm
          function obj = log(obj1)
                         
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = log(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = log(obj1(i,j).sigma);
                  end 
              end

          end
          
          %% Decimal logarithm 
          
          function obj = log10(obj1)
                           
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = log10(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = log10(obj1(i,j).sigma);
                  end 
              end

          end
              
          
          %% exponential
          
          function obj = exp(obj1)
              
              obj = unc_ut;
              
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
              
                  % Update nominal value property
                  obj(i,j).nom_value = exp(obj1(i,j).nom_value);
              
                  % Update dep list 
                  obj(i,j).dep = obj1(i,j).dep;
                  
                  % Compute sigma points         
                  obj(i,j).sigma = exp(obj1(i,j).sigma);
                  end 
              end

          end
          
                    
          %% 
          
           function obj = atan2(obj1,obj2)
               
               % Overloading the Four-quadrant inverse tangent
               % function.Inputs must have compatible size 
               % Declaring both inputs as unc objects
              obj1 = unc_ut(obj1);
              obj2 = unc_ut(obj2);
              
                           
              %check the dim
              %%% one of the objects is scalar 
              if length(obj1)==1
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i = 1:size(obj2,1)
                      for j = 1:size(obj2,2)
                          obj(i,j)= atan2scal(obj1,obj2(i,j)); 
                      end
                  end
              elseif length(obj2)==1
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  
                  for i = 1:size(obj1,1)
                      for j = 1:size(obj1,2)
                          obj(i,j)= atan2scal(obj1(i,j),obj2); 
                      end
                  end
                   %%% vector-vector/ matrix-matrix addition. Objects have the same dim   
              elseif size(obj1)==size(obj2)
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= atan2scal(obj1(i,j),obj2(i,j));                        
                      end                    
                  end
                  %%% vector-matrix addition 
              elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= atan2scal(obj1(i),obj2(i,j));                        
                      end                    
                  end
              elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                  obj(size(obj2,1),size(obj2,2))=unc_ut;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= atan2scal(obj1(j),obj2(i,j));                        
                      end                    
                  end
              elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= atan2scal(obj2(i),obj1(i,j));                        
                      end                    
                  end
              elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                  obj(size(obj1,1),size(obj1,2))=unc_ut;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= atan2scal(obj2(j),obj1(i,j));                        
                      end                    
                  end
                  %%% vector-vector addition
               elseif (isrow(obj1) && iscolumn(obj2))
                  obj(length(obj2),length(obj1))=unc_ut;
                  for i=1:length(obj2)
                      for j=1:length(obj1)
                          obj(i,j)= atan2scal(obj2(i),obj1(j));                        
                      end                    
                  end
               elseif (isrow(obj2) && iscolumn(obj1))
                  obj(length(obj1),length(obj2))=unc_ut;
                  for i=1:length(obj1)
                      for j=1:length(obj2)
                          obj(i,j)= atan2scal(obj1(i),obj2(j));                        
                      end                    
                  end
              else
                  error('Invalid operation. Matrix dimension must agree')
              end
            
           end 
           
           function obj = atan2d(obj1,obj2)
                obj11 = atan2(obj1,obj2);
                obj = obj11/pi*180;
           end 
           
           function obj = atan2scal(obj1,obj2)
              % Method used by operator overloading methods
            obj = unc_ut;            
            % Updating nominal value property
            obj.nom_value = atan2(obj1.nom_value, obj2.nom_value); 
            
            % Updating dep property        
            [~,ia1,ib1]=intersect(obj1.dep,obj2.dep); % in both
            [~,ia2,ib2]=setxor(obj1.dep,obj2.dep);  % only in one
            
             obj.dep=[obj1.dep(ia1),obj1.dep(ia2),obj2.dep(ib2)];         
                       
             obj.sigma=zeros(2,length(obj.dep));  
             
             % Defining the sigma points 
             lboth=length(ia1);
             la=length(ia2);
             lb=length(ib2);
             j=0;
             for i=1:lboth
                 j=j+1;
                 obj.sigma(1,i)=atan2(obj1.sigma(1,ia1(i)),obj2.sigma(1,ib1(i)));
                 obj.sigma(2,i)=atan2(obj1.sigma(2,ia1(i)),obj2.sigma(2,ib1(i)));
             end
                         
             
             for i=1:la
                 j=j+1;
                 obj.sigma(1,j)=atan2(obj1.sigma(1,ia2(i)),obj2.nom_value);
                 obj.sigma(2,j)=atan2(obj1.sigma(2,ia2(i)),obj2.nom_value);
             end
                          
             for i=1:lb
                 j=j+1;
                 obj.sigma(1,j)=atan2(obj1.nom_value,obj2.sigma(1,ib2(i)));
                 obj.sigma(2,j)=atan2(obj1.nom_value,obj2.sigma(2,ib2(i)));
             end
             

           end
          
                   
          %% 
          
          function [obj,ind] = max(obj1)
              
              % max(obj1): Return the uncertainty object that has the
              % maximum nominal value of all the elements of obj1
              
              obj = unc_ut;
              
              if isvector(obj1) %obj1 is a vector
                 [obj.nom_value,ind] = max([obj1.nom_value]);
                 obj = obj1(ind);
                
              elseif ismatrix(obj1)
                  ind = size(obj1,2);
                 for j=1:size(obj1,2) % obj1 is a matrix
                     [obj(j).nom_value,ind(j)] = max([obj1(:,j).nom_value]);
                     obj(j) = obj1(ind(j),j); 
                      
                 end
               end 
                 
          end
         
          %% 
          
          function obj = min(obj1)
              
              % min(obj1): Return the uncertainty object that has the
              % minimum nominal value of all the elements of obj1
              
           obj = unc_ut;
              
              if isvector(obj1) %obj1 is a vector
                 [obj.nom_value,ind] = min([obj1.nom_value]);
                 obj = obj1(ind);
                
              elseif ismatrix(obj1)
                  ind = size(obj1,2);
                 for j=1:size(obj1,2) % obj1 is a matrix
                     [obj(j).nom_value,ind(j)] = min([obj1(:,j).nom_value]);
                     obj(j) = obj1(ind(j),j); 
                      
                 end
               end 
                 
          end
          
          
          %%        
          
          function obj=cumsum(obj1)
              
              % cumsum(obj1): Return the cumulative sum of the elements of obj1  

%               obj=unc_ut; % Preallocation. obj must be of the same size as obj1.
              if isvector(obj1) %obj1 is a vector
                 obj = obj1(1);
                 for i= 2:length(obj1)
                      obj(i) = obj(i-1) + obj1(i);
                 end
                
               elseif ismatrix(obj1)
                 obj=unc_ut;
                 for j=1: size(obj1,2) % obj1 is a matrix
                     obj(1,j) = obj1(1,j);                                       
                     for i = 2:size(obj1,1)
                         obj(i,j) = obj(i-1,j) + obj1(i,j);
                     end
                     
                 end
               end 
              
              
          end 
          
                    
          %% 
          
          function obj=sum(obj1)
              
              % sum(obj1): Sum of elements of obj1
               if isvector(obj1) %obj1 is a vector
                 obj = obj1(1);
                 for i= 2:length(obj1)
                      obj = obj + obj1(i);
                 end
                
               elseif ismatrix(obj1)
                 obj(1,size(obj1,2))=unc_ut;
                 for j=1: size(obj1,2) % obj1 is a matrix
                     obj3 = obj1(1,j);                                       
                     for i = 2:size(obj1,1)
                         obj3 = obj3 + obj1(i,j);
                     end
                     obj(j) = obj3;
                 end
               end 
                               
          end 
 
          
          %% 
          
           function obj = power(obj1,obj2)
               
               % Overloading of Element-wise powers
              
              obj1=unc_ut(obj1);
              obj2=unc_ut(obj2);
              
             %case obj1 is complex and obj2 are real        
                 if any(imag([obj1.nom_value])~=0) && all(imag([obj2.nom_value])==0)  
                        % obtain the absolute value and angle of obj1
                        obj1_r = abs(obj1);
                        obj1_p = angle(obj1);
                        % Moivre's formula for complex numbers         
                        obj=exp(obj2 .* log(obj1_r)) .* ( cos(obj2.*obj1_p) + 1i * sin(obj2.*obj1_p));
                   
                 %case obj1 and obj2 are real 
                 elseif all(imag([obj1.nom_value])==0) && all(imag([obj2.nom_value])==0)
                    
                             %check the dim
                      %%% one of the objects is scalar 
                      if length(obj1)==1
                          obj(size(obj2,1),size(obj2,2))=unc_ut;
                          for i = 1:size(obj2,1)
                              for j = 1:size(obj2,2)
                                obj(i,j) = power_check(obj1,obj2(i,j));   
                              end
                          end
                      elseif length(obj2)==1
                          obj(size(obj1,1),size(obj1,2))=unc_ut;

                          for i = 1:size(obj1,1)
                              for j = 1:size(obj1,2)
                                  obj(i,j) = power_check(obj1(i,j),obj2); 
                              end
                          end
                           %%% vector-vector/ matrix-matrix addition. Objects have the same dim   
                      elseif size(obj1)==size(obj2)
                          obj(size(obj2,1),size(obj2,2))=unc_ut;
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                  obj(i,j)= power_check(obj1(i,j),obj2(i,j));                        
                              end                    
                          end
                          %%% vector-matrix addition 
                      elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                          obj(size(obj2,1),size(obj2,2))=unc_ut;
                          for i=1:size(obj2,1)
                              for j=1:size(obj2,2)
                                  obj(i,j)= power_check(obj1(i),obj2(i,j));                        
                              end                    
                          end
                      elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                          obj(size(obj2,1),size(obj2,2))=unc_ut;
                          for i=1:size(obj2,1)
                              for j=1:size(obj2,2)
                                  obj(i,j)= power_check(obj1(j),obj2(i,j));                        
                              end                    
                          end
                      elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                          obj(size(obj1,1),size(obj1,2))=unc_ut;
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                  obj(i,j)= power_check(obj2(i),obj1(i,j));                        
                              end                    
                          end
                      elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                          obj(size(obj1,1),size(obj1,2))=unc_ut;
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                  obj(i,j)= power_check(obj2(j),obj1(i,j));                        
                              end                    
                          end
                          %%% vector-vector addition
                       elseif (isrow(obj1) && iscolumn(obj2))
                          obj(length(obj2),length(obj1))=unc_ut;
                          for i=1:length(obj2)
                              for j=1:length(obj1)
                                  obj(i,j)= power_check(obj2(i),obj1(j));                        
                              end                    
                          end
                       elseif (isrow(obj2) && iscolumn(obj1))
                          obj(length(obj1),length(obj2))=unc_ut;
                          for i=1:length(obj1)
                              for j=1:length(obj2)
                                  obj(i,j)= power_check(obj1(i),obj2(j));                        
                              end                    
                          end
                      else
                      error('Invalid operation. Matrix dimension must agree')
                    
                      end 
    
                     
                 elseif any(any(imag(obj2)~=0)) 
                     error('Error using .^ Operation with complex exponent is not available.')
                     
                 end 
                        
                  
           end
          
           %%%this function is only used by power. Check is the base is positive, neg or 0. 
           function obj = power_check(obj1,obj2)
               % check if obj1 is positive
               if (obj1.nom_value > 0 )
                   obj = exp(obj2 .* log(obj1));
               elseif (obj1.nom_value < 0 )
                   obj =((-1).^gmv(obj2)) .* exp(obj2 .* log(abs(obj1)));
               elseif (obj1.nom_value == 0 )  
                   obj = 0;
               end
                   
           end
          
           %%%%%%%
          

          %% 
          
          function obj = mpower(obj1,obj2)
              
              % This method implements the overloading of matrix power called for the syntax 'obj1 ^ obj2'
              
              if (obj2 == 0) && (size(obj1,1) == size(obj1,2))
                  obj = obj1*inv(obj1);
                  
              elseif  isa(obj2,'double') && ((abs(round(obj2)-obj2))<eps*2)
                obj=unc_ut(obj1);
                for i=2:int16(obj2)
                      obj=obj*obj1;
                end             
              else
                  error('Sorry - this operation is currently not supported')
              end 
          end
          
          
          %% 
          
          function obj = sqrt(obj1)
              
              % sqrt(obj1): Return square root of obj1
              
              obj = power(obj1,0.5);
                 
          end
          
          
          %% 
          
          function obj = realpow(obj1,obj2)
              
              % This method implements the overloading of the element-by-element powers
              obj = unc_ut;
              L = max(size(obj1),size(obj2));
              if ismatrix(obj1) || ismatrix(obj2)
                   for i=1:L(1)
                       if size(obj1)==size(obj2)
                           obj(i,:)= power(obj1(i,:),obj2(i,:));
                           if ~isreal(obj(i).nom_value)
                               error('Error using realpow. Realpow produced complex result.')
                           end
                       elseif isscalar(obj1) && ~isscalar(obj2)
                           obj(i,:)=power(obj1,obj2(i,:));
                           if ~isreal(obj(i).nom_value)
                               error('Error using realpow. Realpow produced complex result.')
                           end
                       elseif ~isscalar(obj1) && isscalar(obj2)
                           
                           obj(i,:)= power(obj1(i,:),obj2);
                           if ~isreal(obj(i).nom_value)
                               error('Error using realpow. Realpow produced complex result.')
                           end
                       else
                           error('Error using realpow. Matrix dimensions must agree.')
                       end
                   end
                       
              end
          end
          
          %% Get Covariance and Correlation matrices 
          function [CX_cov] = get_cov_mat(X)

           % Return the covariance matrix CX_cov of the uncertainty 
           % objects contained in the array X. 
              global kappa beta  
              if size(X,1)==1
                   CX_cov = zeros(length(X));

                      for i=1:length(X)
                           for j=1:length(X)

                              obj1=unc_ut(X(i));
                              obj2=unc_ut(X(j));

                              % find if obj1 and obj2 have dependencies in comum 
                              [c1,ia1,ib1] = intersect(obj1.dep,obj2.dep);
                              if isempty(c1)
                                 CX_cov(i,j) = 0;
                              else
                              n = length(c1);  % dimensionality
                              alpha = 1/sqrt(n+kappa);      
                              lambda=alpha^2*(n+kappa)-n;% scaling parameter
                              W0=lambda/(n+lambda);  % weight sigma point 0
                              Wi=1/2/(n+lambda); 
                              
                              value_transformed1 = sum(obj1.sigma(:,ia1),'all')*Wi + obj1.nom_value*W0;
                              value_transformed2 = sum(obj2.sigma(:,ib1),'all')*Wi + obj2.nom_value*W0;
                              
                       
                              CX_cov(i,j) = Wi*sum((obj1.sigma(:,ia1)-value_transformed1).*(obj2.sigma(:,ib1)- value_transformed2),'all') + (W0+1-alpha^2+beta)*(obj1.nom_value-value_transformed1)*(obj2.nom_value-value_transformed2);
                              
                             end 
                           end
                      end
                      
               else 
                   error('Input argument must be a one-row array of uncertainty objects!')
               end
          end
          
          
          %% FFT and IFFT
           function obj = fft(obj1)
                
                N = numel(obj1);
                if isrow(obj1)
                    obj1 = obj1.';
                end 
                % check if the length of the input vector is a power of 2
                if ~(floor(log2(N))==log2(N)) 
                    N1 = power(2,ceil(log2(N))); % new dimension( power of 2)
                    obj1 = [obj1;unc_ut(zeros((N1-N),1))];
                    N = N1; % define the new dimension of N
                end 

                xp = obj1(1:2:end);
                xpp = obj1(2:2:end);
                if N>=8
                     Xp = myfft(xp);
                     Xpp = myfft(xpp);
                     obj = unc_ut(zeros(N,1));
                     Wn = exp(-1i*2*pi*((0:N/2-1)')/N);
                     tmp = Wn .* Xpp;

                     obj = [(Xp + tmp);(Xp - tmp)];
                else
                     switch N
                         case 2
                             obj = [1 1;1 -1]*obj1;
                         case 4
                             obj = [1 0 1 0; 0 1 0 -1i; 1 0 -1 0;0 1 0 1i]*[1 0 1 0;1 0 -1 0;0 1 0 1;0 1 0 -1]*obj1;

                     end
                end
           end
           
           function obj = ifft(obj1)
                
                N = numel(obj1);
                Xc = conj(obj1); % complex conjugate of X 
                xc = myfft(Xc); % compute fft of N points of the complex conj
                obj = 1/N * conj(xc);
           end
           
                %% 1D Linear interpolation
        function Yp = interp1(X,Y,Xp)
            % overloading 1D linear interpolation 
            N = length(X);
            X = unc_ut(X);
            Y = unc_ut(Y);
            Xp = unc_ut(Xp); 

            if isvector(Y) && length(Y)== N
                Yp = unc_ut(zeros(size(Xp)));
                ind_0 = zeros(size(Xp));
                ind_1 = zeros(size(Xp));
                for i = 1:length(Xp)

                    if isempty(max(X(Xp(i).nom_value > [X(:).nom_value]))) || isempty( min(X(Xp(i).nom_value< [X(:).nom_value])))
                       Yp(i) = NaN; 
                    else
                     % find the index of the lower and upper points 
                    ind_0(i) = find(X(:)== max(X(Xp(i).nom_value> [X(:).nom_value])));
                    ind_1(i) = find(X(:)== min(X(Xp(i).nom_value< [X(:).nom_value])));   

                    % Find the coordinates of the points. xo < Xp(i)< x1
                        x0 = X(ind_0(i));
                        x1 = X(ind_1(i));
                        y0 = Y(ind_0(i));
                        y1 = Y(ind_1(i));

                        Yp(i) = (y0*(x1-Xp(i)) + y1*(Xp(i)- x0))/(x1-x0);
                    end 
                end 

            else 
                error('X and Y must have the same length');

            end
        end
                 
         
          %% 
          
           function obj = gt(obj1,obj2)
           
              % This method implements the overloading of the greater 
              % than operator called for the syntax 'obj1 > obj2'
                       
           obj=unc_ut(0,0);
           
           if isa(obj1,'unc_ut') && isa(obj2,'unc_ut')
               obj1=unc_ut(obj1);
               obj2=unc_ut(obj2);
               
               h=obj1-obj2;
               
               obj.nom_value=double(h.nom_value>0); % obj.value=1 if h.value>0 
                                            % obj.value=0 if h.value<0
               
               obj.std_unc=normcdf(h.value,0,h.std_unc); % obj.std_unc= area under the standard
                                                         % normal distribution curve from
                                                         % -infinity to (h.value - 0)/h.std_unc
               
           
               % obj.value represents the logical decision (in double)
               % obj.std_unc represents the probability of obj1 beeing larger than obj2
            
           else 
               % Comparison with a constant
               if (isa(obj1,'unc_ut') && isa(obj2,'double')) || (isa(obj1,'double') && isa(obj2,'unc_ut'))
                   
                   if isa(obj1,'unc_ut')
                       
                       obj.value=double(obj1.value>obj2);
                       obj.std_unc=1-normcdf(obj2,obj1.value,obj1.std_unc);
                       
                   else
                       if isa(obj2,'unc_ut')
                          
                          obj.value=double(obj2.value<obj1);
                          obj.std_unc=normcdf(obj1,obj2.value,obj2.std_unc);
                          
                       end
                
                   end
               end
                       
           end         
                   
                           
           end    

        
           %% obj1 less than obj2
          
           function obj = lt(obj1,obj2)
               
%                This method implements the overloading of the less than operator 
%                called for the syntax 'obj1 < obj2'
               
             obj=unc_ut(0,0);
           
             if isa(obj1,'unc_ut') && isa(obj2,'unc_ut')
               obj1=unc_ut(obj1);
               obj2=unc_ut(obj2);
             
               h=obj1-obj2;
               
               obj.value=double(h.value<0); % obj.value=1 if h.value<0 
                                            % obj.value=0 if h.value>0
                             
               obj.std_unc=normcdf(-h.value,0,h.std_unc); % obj.std_unc= area under the standard
                                                         % normal distribution curve from
                                                         % -infinity to (h.value - 0)/h.std_unc                                 
                                                         
               % obj.value represents the logical decision (in double)
               % obj.std_unc represents the probability of obj1 beeing less than obj2
               
             else 
               % Comparison with a constant
                if (isa(obj1,'unc_ut') && isa(obj2,'double')) || (isa(obj1,'double') && isa(obj2,'unc_ut'))
                   
                   if isa(obj1,'unc_ut')
                       
                       obj.value=double(obj1.value<obj2);
                       obj.std_unc=normcdf(obj2,obj1.value,obj1.std_unc);
                       
                   else
                       if isa(obj2,'unc_ut')
                          
                          obj.value=double(obj2.value>obj1);
                          obj.std_unc=1-normcdf(obj1,obj2.value,obj2.std_unc);
                         
                       end
                
                   end
                end
                                                       
             end  
           end

          
          %%
          %%%%%%%%%%%%%% disp_dep() method modified to accept n-dimensional array %%%%%%%%%%%%%%%%%%%
          %display an array containing the names and values of the variable on which an uncertainty object depends%%   
          function  disp_dep(obj)
                 
                 % This method displays an array containing the variables
                 % on which obj depends
                
                depend = cell(size(obj)); % for saving the dependency list of each element in obj    
%                 fprintf('\n')
                compdep = {}; % cell array used to compare the dependent objects.
                for i=1:size(obj,1) % index the rows
                    for j=1:size(obj,2) % index the colums   
                       
                       for k = 1:length(obj(i,j).dep)

                            if isempty(obj(i,j).dep(1,k).name)
                                dep_name = 'not named';
                                dep_value = obj.GenerateDispStr(obj(i,j).dep(1,k).value,obj(i,j).dep(1,k).std_unc);
                                dep_inf = [dep_name,' = ' ,dep_value];
                                                                 
                            else

                                dep_name = obj(i,j).dep(1,k).name{1,1};                         
                                dep_value = obj.GenerateDispStr(obj(i,j).dep(1,k).value,obj(i,j).dep(1,k).std_unc);
                                dep_inf = [dep_name,' = ' ,dep_value];
                                                                                       
                            end 
                        
                        % evaluate if the current object is repeated                     
                            if isempty(compdep)| ~strcmp(dep_inf,compdep)
                                
                                compdep{end+1} = dep_inf;  %store the dep
                                Dep_inf{k,:} = dep_inf;
                                
                                depend{i,j} = char(Dep_inf{k,:});
                                disp(depend{i,j})
                                fprintf('\n')
                            end 
                       end
                       
                      
                    end
               end
                        
          end

                                                               
          %%

          function disp(obj)
              
              % Overloading of the Matlab disp() function
                    
             for i =1 :size(obj,1)
                for j=1 : size(obj,2)
                           
                     if imag(obj(i,j).nom_value)~=0 
                           
                            obj_X = real(obj(i,j));
                            obj_Y= imag(obj(i,j));
                            [objX_val,objX_stdUnc] = eval_std_unc(obj_X); 
                            [objY_val,objY_stdUnc] = eval_std_unc(obj_Y);
                            
                            if obj_Y.nom_value >= 0
                           
                               fprintf(['    \t',obj_X.GenerateDispStr(objX_val,objX_stdUnc) ,' + '...
                                         obj_Y.GenerateDispStr(objY_val,objY_stdUnc) , ' * i ']);
                           
                            else
                              fprintf(['    \t',obj_X.GenerateDispStr(objX_val,objX_stdUnc) ,' - '...
                                        obj_Y.GenerateDispStr(abs(objY_val),objY_stdUnc) , ' * i ']);
                            end
                                
                   
                       
                     elseif imag(obj(i,j).nom_value)==0
                                  
                                         [objX_val,objX_stdUnc] = eval_std_unc(obj(i,j)); 
                                         fprintf(['    \t',obj.GenerateDispStr(objX_val,objX_stdUnc)]);
                               
                               
                            
                     end
                            
                end  
                
                fprintf('\n')
                
              end
             
          end
            
                     
          
%% 

function Z = complex(x,y)
    
    % Z = complex(x,y): Create complex uncertainty objects Z out of real uncertainty objects
    % x and y
    
   [r1,c1]=size(x); 
   
   [r2,c2]=size(y);
   
   if (r1==r2) && (c1==c2)  
    Z = unc_ut; 
       for i=1:r1
           for j=1:c1
               
               if (imag(x(i,j).nom_value)==0) && (imag(y(i,j).nom_value)==0)
       
                Z(i,j)=unc_ut; 
                
                Z(i,j).nom_value= x(i,j).nom_value + 1i * y(i,j).nom_value; % Mean value assignment
                % Updating dep property        
                [~,ia1,ib1]=intersect(x(i,j).dep,y(i,j).dep); % in both
                [~,ia2,ib2]=setxor(x(i,j).dep,y(i,j).dep);  % only in one

                Z(i,j).dep=[x(i,j).dep(ia1),x(i,j).dep(ia2),y(i,j).dep(ib2)];         

                Z(i,j).sigma=zeros(2,length(Z(i,j).dep));  

                 % Defining the sigma points 
                 lboth=length(ia1);
                 la=length(ia2);
                 lb=length(ib2);
                 j1=0;
                 for i1=1:lboth
                     j1=j1+1;
                     Z(i,j).sigma(1,i1)=x(i,j).sigma(1,ia1(i1))+ 1i * y(i,j).sigma(1,ib1(i1));
                     Z(i,j).sigma(2,i1)=x(i,j).sigma(2,ia1(i1))+ 1i * y(i,j).sigma(2,ib1(i1));
                 end


                 for i1=1:la
                     j1=j1+1;
                     Z(i,j).sigma(1,j1)=x(i,j).sigma(1,ia2(i1))+ 1i * y(i,j).nom_value;
                     Z(i,j).sigma(2,j1)=x(i,j).sigma(2,ia2(i1))+ 1i * y(i,j).nom_value;
                 end

                 for i1=1:lb
                     j1=j1+1;
                     Z(i,j).sigma(1,j1)=x(i,j).nom_value+ 1i * y(i,j).sigma(1,ib2(i1));
                     Z(i,j).sigma(2,j1)=x(i,j).nom_value+ 1i * y(i,j).sigma(2,ib2(i1));
                 end
                
                
               else
                   error('Input arrays must be real objects.')
               end
                   
           end
       end
   
   else
       
       error('Input arrays must have the same size.')
       
   end
    
    
    
    
end


%% 
          
 function X=real(Z)
      
     % real(Z): Return the real part of Z
     X = unc_ut;
     for i=1:size(Z,1)
         for j=1:size(Z,2)
             
             % Update nominal value property
             X(i,j).nom_value = real(Z(i,j).nom_value);
              
             % Update dep list 
             X(i,j).dep = Z(i,j).dep;
                  
             % Compute sigma points         
             X(i,j).sigma = real(Z(i,j).sigma);
             
         end
     end
             
 end

 
%% 
          
 function Y=imag(Z)
     
     % imag(Z): Return the imaginary part of Z
     Y = unc_ut;
     for i=1:size(Z,1)
         for j=1:size(Z,2)
             
             % Update nominal value property
             Y(i,j).nom_value = imag(Z(i,j).nom_value);
              
             % Update dep list 
             Y(i,j).dep = Z(i,j).dep;
                  
             % Compute sigma points         
             Y(i,j).sigma = imag(Z(i,j).sigma);

         end
     end
   
 end
   
 
%% 
          
 function obj=conj(Z)
       
     % This method takes the complex conjugate of each
     % entry of the matrix Z
 
       obj=unc_ut(Z);
       for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                 obj(i,j).nom_value= conj(Z(i,j).nom_value);
                 obj(i,j).sigma(1,:)=conj(Z(i,j).sigma(1,:));
                 obj(i,j).sigma(2,:)=conj(Z(i,j).sigma(2,:));
                 obj(i,j).dep = Z(i,j).dep;
         
         end
       end
       
             
 end  
 
 
 
 %% 
          
 function obj=abs(Z)
   % overloading the function abs() that returns the absolute value of the uncertainty object Z. When
   % Z is complex it returns the complex modulus (magnitude) of the elements of Z.
   obj = unc_ut();
       for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                 obj(i,j).nom_value = abs(Z(i,j).nom_value);
                 obj(i,j).sigma(1,:) = abs(Z(i,j).sigma(1,:));
                 obj(i,j).sigma(2,:) = abs(Z(i,j).sigma(2,:));
                 obj(i,j).dep = Z(i,j).dep;
                 
                
         end
       end
       
             
 end  
 
 function obj=angle(Z)
      % overloading the function angle() that returns the phase angles, in radians, of the uncertainty object Z.
      obj = unc_ut();
       for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                obj(i,j).nom_value = angle(Z(i,j).nom_value);
                obj(i,j).sigma(1,:) = angle(Z(i,j).sigma(1,:));
                obj(i,j).sigma(2,:) = angle(Z(i,j).sigma(2,:));
                obj(i,j).dep = Z(i,j).dep;
                
         end
       end
       
             
 end  
 
function disp_contribution(obj1)
     warning('The function disp_contribution is not implemented in the Unscented Transform approach!')
 
 end  
 
end

          %%
          methods (Static)
              %%%%%%%%%%%%%%%%%function GenerateDispStr modified%%%%%%
          function DispString = GenerateDispStr(value,std_unc)
               
               % This method is used by the disp() method to customize 
               % the display of uncertainty objects
               
                global Number_of_Digits_to_Display;
                d=Number_of_Digits_to_Display-1; 
                
                
                % 1) Clip the standard uncertainty to max 3 digits
                if (std_unc==0) || (log10(abs(value)/abs(std_unc))>12)
                    digit=-12;
                else
                    digit = floor(log10(abs(std_unc))); % take the integer part no matter is the value of the decimal part
                end
                
                
                
                disp_std_unc = round(std_unc / 10^(digit-d)); % round the value into the closest integer.
                
                
                % 2) Clip the value to show the last d+1 values 
                disp_value = round(value / 10^(digit-d))*10^(digit-d);
                
                
                if disp_std_unc==0
                        eval_str = '[num2str(disp_value),''(0)'']';
                        DispString = eval (eval_str);
                else
                     
                    if (digit<-4) && (abs(disp_value)<0.001) && (disp_value~=0)
                        % show the result in scientific notation if there
                        % are more than 2 leading zeros
                        eval_str = sprintf ('[sprintf(''%%0.0%ie'',disp_value),''('',num2str(disp_std_unc),'')'']',abs(digit)-abs(floor(log10(abs(value))))+d);
                        DispString = eval (eval_str);
                    elseif (digit<0) 
                        eval_str = sprintf ('[sprintf(''%%0.0%if'',disp_value),''('',num2str(disp_std_unc),'')'']',abs(digit)+d);
                        DispString = eval (eval_str);
                    else
                        if digit<d
                            disp_std_unc = disp_std_unc*10^(digit-d); %#ok<NASGU>
                             eval_str = sprintf ('[sprintf(''%%0.0%if'',disp_value),''('',sprintf(''%%0.0%if'',disp_std_unc),'')'']',d-digit,d-digit);
                            DispString = eval (eval_str); 
                        else
                            disp_std_unc = disp_std_unc*10^(digit-d);
                            DispString = [num2str(disp_value),'(',num2str(disp_std_unc),')'];
                        end
                    end
                end
          end
           
          
             
  
       end
end

         
          
