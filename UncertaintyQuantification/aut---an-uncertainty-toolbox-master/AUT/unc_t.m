classdef unc_t < handle
  % AUT - An Uncertainty Toolbox V 0.91
  % Date: 27/11/2019
  % Type "help disclaimer" for more information

   
  
  %% Properties of the uncertainty class
   
       properties  
       % Define the properties of the unc_t object to be created by the Class Constructor 
          values;  % samples 
         
       end
       
                   
  %% Methods of the uncertainty class    
       
       methods           
  
           %% ------ Class constructor ------- %
  
       function obj = unc_t(arg1,arg2,arg3)
                           
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
          
          global MCSAMPLES;
          if isempty(MCSAMPLES)
               MCSAMPLES=10000;
               fprintf(['Using %d MC Samples\n',13],MCSAMPLES);
          end
          
          
          if (nargin==0)

                          obj.values=zeros(1,MCSAMPLES);
                                  
          elseif (nargin==1)    
              
             % Syntax(1): Value only
             
             if isa(arg1,'double') 
                 
                       obj(size(arg1,1),size(arg1,2))=unc_t;
                       for i=1:size(arg1,1)
                          for j=1:size(arg1,2)

                          obj(i,j).values=arg1(i,j)*ones(1,MCSAMPLES);
                          
                          end
                       end
                       
             % Syntax(2): unc_t object only
             elseif isa(arg1,'unc_t')
                       
                       obj = arg1;
                       
             end
               % Syntax (2):  unc_t(Value, uncertainty)
          elseif (nargin >= 2)      
                 if isa(arg1,'double')&& isa(arg2, 'double') && size(arg1,1)==size(arg2,1) && size(arg1,2)==size(arg2,2) && all(all(isreal(arg1)))
                      obj(size(arg1,1),size(arg1,2))=unc_t;
                      for i=1:size(arg1,1)
                          for j=1:size(arg1,2) 

                               obj(i,j).values = normrnd(arg1(i,j),arg2(i,j),1,MCSAMPLES);
                                                              
                          end
                      end
                 elseif isa(arg1,'double')&& isa(arg2, 'double') && ismatrix(arg2)&& size(arg2,1)==length(arg1) && size(arg2,2)==length(arg2) && all(all(isreal(arg1)))
                      obj(size(arg1,1),size(arg1,2))=unc_t;
                      
                        R=mvnrnd(arg1,arg2,MCSAMPLES);

                        for i=1:length(arg1)
                                obj(i).values=R(:,i)';

                        end
                                       
                 end
           else                 
                 error('Incorrect type of argument(s)');    

          end
          

       end

       
    %% Overloading basic operations 
 
       function obj = plus(obj1,obj2)
           
            % Overloading of the Addition Operator
           global MCSAMPLES;
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj = unc_t;
                       
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));

           for kk =1:MCSAMPLES
 
               for i1 = 1:size(obj1,1)
                   for j1 = 1:size(obj1,2)
                       obj11(i1,j1) = obj1(i1,j1).values(kk);
                   end 
               end 
               
               for i2 = 1:size(obj2,1)
                   for j2 = 1:size(obj2,2)
                       obj22(i2,j2) = obj2(i2,j2).values(kk);
                   end 
               end 
              
                    
               obj_value = obj11 + obj22;
               
               
                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 
           end 
                   
       end 
           
                              
            %% 
            
       function obj = minus(obj1,obj2)
           
           % Overloading of the Subtraction Operator
           global MCSAMPLES;
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj = unc_t;
                       
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));

           for kk =1:MCSAMPLES
 
               for i1 = 1:size(obj1,1)
                   for j1 = 1:size(obj1,2)
                       obj11(i1,j1) = obj1(i1,j1).values(kk);
                   end 
               end 
               
               for i2 = 1:size(obj2,1)
                   for j2 = 1:size(obj2,2)
                       obj22(i2,j2) = obj2(i2,j2).values(kk);
                   end 
               end 
              
                    
               obj_value = obj11 - obj22;
               
               
                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 
           end 
               
                     
       end 
            
           
           %% 
           
       function obj = times(obj1,obj2)
           
           % Overloading times Operator
            
           global MCSAMPLES;
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj = unc_t;
                       
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));

           for kk =1:MCSAMPLES
 
               for i1 = 1:size(obj1,1)
                   for j1 = 1:size(obj1,2)
                       obj11(i1,j1) = obj1(i1,j1).values(kk);
                   end 
               end 
               
               for i2 = 1:size(obj2,1)
                   for j2 = 1:size(obj2,2)
                       obj22(i2,j2) = obj2(i2,j2).values(kk);
                   end 
               end 
              
                    
               obj_value = obj11 .* obj22;
               
               
                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 
           end 
           
       end 
          

          %% 
          
       function obj =mrdivide(obj1,obj2)
           
           % Overloading of the division operator
           
           global MCSAMPLES;
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj = unc_t;
                       
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));

           for kk =1:MCSAMPLES
 
               for i1 = 1:size(obj1,1)
                   for j1 = 1:size(obj1,2)
                       obj11(i1,j1) = obj1(i1,j1).values(kk);
                   end 
               end 
               
               for i2 = 1:size(obj2,1)
                   for j2 = 1:size(obj2,2)
                       obj22(i2,j2) = obj2(i2,j2).values(kk);
                   end 
               end 
              
                    
               obj_value = obj11 / obj22;
               
               
                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 
           end 
                  
     end
           
          %%  matrix operations 
          
          function obj =mldivide(obj1,obj2)
              
%               This method implements the overloading of matrix left division  
%               for two arrays obj1 and obj2 of uncertainty objects.

          global MCSAMPLES;
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj = unc_t;
                       
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));

           for kk =1:MCSAMPLES
 
               for i1 = 1:size(obj1,1)
                   for j1 = 1:size(obj1,2)
                       obj11(i1,j1) = obj1(i1,j1).values(kk);
                   end 
               end 
               
               for i2 = 1:size(obj2,1)
                   for j2 = 1:size(obj2,2)
                       obj22(i2,j2) = obj2(i2,j2).values(kk);
                   end 
               end 
              
                    
               obj_value = obj11 \ obj22;
               
               
                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 
           end 
         end
          
           
            %% 
          
          function obj=rdivide(obj1,obj2)
              
%               This method implements the overloading of the right element-wise 
%               division operator (./) for two arrays obj1 and obj2 of uncertainty objects.

           global MCSAMPLES;
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj = unc_t;
                       
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));

           for kk =1:MCSAMPLES
 
               for i1 = 1:size(obj1,1)
                   for j1 = 1:size(obj1,2)
                       obj11(i1,j1) = obj1(i1,j1).values(kk);
                   end 
               end 
               
               for i2 = 1:size(obj2,1)
                   for j2 = 1:size(obj2,2)
                       obj22(i2,j2) = obj2(i2,j2).values(kk);
                   end 
               end 
              
                    
               obj_value = obj11 ./ obj22;
               
               
                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 
           end 
% %             

          end
          
          
          %% 
          
          function obj=ldivide(obj1,obj2)
              
%               This method implements the overloading of the left 
%               element-wise division operator (.\) for two arrays obj1 and obj2 of uncertainty objects.

           global MCSAMPLES;
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj = unc_t;
                       
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));

           for kk =1:MCSAMPLES
 
               for i1 = 1:size(obj1,1)
                   for j1 = 1:size(obj1,2)
                       obj11(i1,j1) = obj1(i1,j1).values(kk);
                   end 
               end 
               
               for i2 = 1:size(obj2,1)
                   for j2 = 1:size(obj2,2)
                       obj22(i2,j2) = obj2(i2,j2).values(kk);
                   end 
               end 
              
                    
               obj_value = obj11 .\ obj22;
               
               
                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 
           end
            
          end         
          
                    
          
           %%
            
        function obj = mtimes(obj1,obj2)
            
            % ------------------------------------------------------------ %
            %  mtimes(obj1,obj2): Overloading of matrix multiplication operator
            %  for uncertainty objects obj1 and obj2
            % ------------------------------------------------------------ %
            
            global MCSAMPLES ;
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj11 = zeros(size(obj1,1),size(obj1,2));
           obj22 = zeros(size(obj2,1),size(obj2,2));
           obj = unc_t;
           for kk =1 : MCSAMPLES
           
               for i1 = 1:size(obj1,1)
                   for j1 = 1:size(obj1,2)
                       obj11(i1,j1) = obj1(i1,j1).values(kk);
                   end 
               end 
               
               for i2 = 1:size(obj2,1)
                   for j2 = 1:size(obj2,2)
                       obj22(i2,j2) = obj2(i2,j2).values(kk);
                   end 
               end 
              
                    
               obj_value = obj11 * obj22;
               
               for i = 1:size(obj_value,1)
                   for j = 1:size(obj_value,2)
                       obj(i,j).values(kk) = obj_value(i,j);
                   end 
               end 
           end
               
       end 
        
       
         %% 
         
          function obj = prod(obj1)
              
               % obj=prod(obj1): Product of the elements of the object obj1. If obj1 is a
               % matrix of uncertainty objects, obj is a row vector with the product over each column.
               
               global MCSAMPLES;
               obj1=unc_t(obj1);
               obj11 = zeros(size(obj1,1),size(obj1,2));
               obj = unc_t;
               
                for kk = 1:MCSAMPLES


                   for i1 = 1:size(obj1,1)
                       for j1 = 1:size(obj1,2)
                           obj11(i1,j1) = obj1(i1,j1).values(kk);
                       end 
                   end 


                   obj_value = prod(obj11);

                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 
                   
                end 
        
          end
         
          
          %% 
                    
        function obj = inv(obj1)
              
%              
                    
                    global MCSAMPLES
                    obj=unc_t(zeros(size(obj1,1)));
                    
                    matSamples = zeros(size(obj1,1)); % matrix of samples
                    
                     for kk = 1: MCSAMPLES
                        for i = 1 : size(obj1,1) 
                         for j = 1: size(obj1,1)
                             matSamples(i,j) = obj1(i,j).values(kk);                         
             
                         end 
                        end 
                        
                        InvMatrix = inv(matSamples);
                        
                      
                        for i1 = 1 : size(obj1,1) 
                         for j1 = 1: size(obj1,1)
                        
                            obj(i1,j1).values(kk)= InvMatrix(i1,j1);
                            
                         end
                        end
                        
                     end 
 
        end
          
        
        
        %% 
        
        function obj = diag(obj1)
            
            % diag(obj1): Return the diagonal elements of obj1
            
          [r1,c1] = size(obj1);
          obj = unc_t;
                  
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
          obj2 = unc_t;
          for i=1:r1
              for j=1:c1
                   
                  obj2(j,i)=obj1(i,j);  
                      
              end
          end  
           
          obj = conj(obj2);  
          
           
        end       
        
        

        %% 
        
        function [determ]=det(obj)
            
                    determ=unc_t();
                    
                    matSamples = zeros(size(obj)); % matrix of samples
                    global MCSAMPLES
                    for kk = 1: MCSAMPLES
                        for i = 1 : size(obj,1) 
                         for j = 1: size(obj,2)
                             matSamples(i,j) = obj(i,j).values(kk);                         
             
                         end 
                        end 
                        
                        determ.values(kk) = det(matSamples);
  
                    end 
            
         end
          
                   
         
          %% 
          
          

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
              
              % gmv(obj): Get the matrix of value proprties of the
              % array obj of uncertainty objects 
              value = zeros(size(obj));
            
              for i=1:size(obj,1)
                  for j=1:size(obj,2)
                    value(i,j) = mean(obj(i,j).values);
                  end 
              end
              
          end
          

          %% 
          
          function std_unc_t=gmu(obj)
              
              % gmu(obj): Get the matrix of std_unc_t proprties of the
              % array obj of uncertainty objects
              std_unc_t = zeros(size(obj));
              
              for i=1:size(obj,1)
                  for j=1:size(obj,2)
                    std_unc_t(i,j) = std(obj(i,j).values);
                  end 
              end
              
                              
          end

          
                    

 
          %% sinus
          function obj=sin(obj1)
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = sin(obj1(i,j).values);
                  end 
              end
              
                            
          end
          
          function obj=sind(obj1)
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = sind(obj1(i,j).values);
                  end 
              end
              
                            
          end
                  
          %% cosine
          function obj=cos(obj1)
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = cos(obj1(i,j).values);
                  end 
              end
              
          end
          
          function obj=cosd(obj1)
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = cosd(obj1(i,j).values);
                  end 
              end
              
          end
          
          %% tan
          function obj = tan(obj1)
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = tan(obj1(i,j).values);
                  end 
              end
              

          end
          
          function obj = tand(obj1)
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = tand(obj1(i,j).values);
                  end 
              end
              

          end
                               
          
          %% Inverse cosine: acos
          function obj=acos(obj1)
              
             obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = acos(obj1(i,j).values);
                  end 
              end
              
          end
          
          function obj=acosd(obj1)
              
             obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values = acosd(obj1(i,j).values);
                  end 
              end
               
          end
          
          %% Inverse sine: asin
          function obj=asin(obj1)
             
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = asin(obj1(i,j).values);
                  end 
              end
                                     
          end
          
          function obj=asind(obj1)
             
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = asind(obj1(i,j).values);
                  end 
              end
                                     
          end
          %% Inverse tangent: atan
          function obj=atan(obj1)
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = atan(obj1(i,j).values);
                  end 
              end
 
          end
          
          function obj=atand(obj1)
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = atand(obj1(i,j).values);
                  end 
              end
 
          end
          
          %% Hyperbolic cosine: cosh
          function obj = cosh(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = cosh(obj1(i,j).values);
                  end 
              end
 
          end
          
          %% Hyperbolic sine: sinh
          function obj = sinh(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = sinh(obj1(i,j).values);
                  end 
              end
 
              
          end
          
          %% Hyperbolic tangent:tanh
          function obj = tanh(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = tanh(obj1(i,j).values);
                  end 
              end
 
          end
          
          %% Hyperbolic cosine inverse: acosh
          function obj = acosh(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = acosh(obj1(i,j).values);
                  end 
              end

          end
          
          %% Hyperbolic sine inverse: asinh
          function obj = asinh(obj1)
                     
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = asinh(obj1(i,j).values);
                  end 
              end

          end
          
          %% Hyperbolic tangent inverse: atanh
          function obj = atanh(obj1)
             
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = atanh(obj1(i,j).values);
                  end 
              end

          end
          
          %% Cotangent : cot
          function obj = cot(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = cot(obj1(i,j).values);
                  end 
              end

          end
          
          function obj = cotd(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = cotd(obj1(i,j).values);
                  end 
              end

          end
          
          %% Inverse cotangent: acot
          function obj = acot(obj1)
             
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = acot(obj1(i,j).values);
                  end 
              end

          end
          
          function obj = acotd(obj1)
             
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = acotd(obj1(i,j).values);
                  end 
              end

          end
          
          %% Hyperbolic cotangent: coth
          function obj = coth(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = coth(obj1(i,j).values);
                  end 
              end

          end
          
          %% Natural logarithm
          function obj = log(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = log(obj1(i,j).values);
                  end 
              end

          end
          
          %% Decimal logarithm 
          
          function obj = log10(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = log10(obj1(i,j).values);
                  end 
              end

          end
              
          
          %% exponential
          
          function obj = exp(obj1)
              
              obj = unc_t;
              for i=1:size(obj1,1)
                  for j=1:size(obj1,2)
                    obj(i,j).values  = exp(obj1(i,j).values);
                  end 
              end

          end
          
                    
          %% 
          
           function obj = atan2(obj1,obj2)
               
               % Overloading the Four-quadrant inverse tangent function
               % obj1 and obj2 must have compatible size
               
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));
           obj = unc_t;
           
               global MCSAMPLES 
               for kk =1:MCSAMPLES

                   for i1 = 1:size(obj1,1)
                       for j1 = 1:size(obj1,2)
                           obj11(i1,j1) = obj1(i1,j1).values(kk);
                       end 
                   end 

                   for i2 = 1:size(obj2,1)
                       for j2 = 1:size(obj2,2)
                           obj22(i2,j2) = obj2(i2,j2).values(kk);
                       end 
                   end 


                   obj_value = atan2(obj11,obj22);

                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 

               end 
           end
             
           function obj = atan2d(obj1,obj2)
               
               % Overloading the Four-quadrant inverse tangent
               % function(degree).
               % obj1 and obj2 must have compatible size
               
           obj1=unc_t(obj1);
           obj2=unc_t(obj2);
           
           obj11 = zeros(size(obj1));
           obj22 = zeros(size(obj2));
           obj = unc_t;
           
               global MCSAMPLES 
               for kk =1:MCSAMPLES

                   for i1 = 1:size(obj1,1)
                       for j1 = 1:size(obj1,2)
                           obj11(i1,j1) = obj1(i1,j1).values(kk);
                       end 
                   end 

                   for i2 = 1:size(obj2,1)
                       for j2 = 1:size(obj2,2)
                           obj22(i2,j2) = obj2(i2,j2).values(kk);
                       end 
                   end 


                   obj_value = atan2d(obj11,obj22);

                   for i = 1:size(obj_value,1)
                       for j = 1:size(obj_value,2)
                           obj(i,j).values(kk) = obj_value(i,j);
                       end 
                   end 

               end 
             end
          
                   
          %% 
          
          function [obj, ind] = max(obj1)
              
              % max(obj1): Return the uncertainty object that has the
              % maximum nominal value of all the elements of obj1
              
              obj= unc_t;
              obj1_value = gmv(obj1);
              
              [~,ind]= max(obj1_value);
              obj.values = obj1(ind).values;
              
          end
          
          function [obj,ind_prob] = max_RV(obj1)
             %% Ind_prob is an array containing the probability for each element of obj1 be the maximum.
             % Obj is a random variable for the maximum.
             
              obj = unc_t;
              ind_prob = zeros(1,length(obj1));
                             
              global MCSAMPLES
              obj11 = zeros(length(obj1),MCSAMPLES);
  
              for i = 1:length(obj1)
                 [obj11(i,:)] = obj1(i).values;
              end
        
              [obj_rv, I] = max(obj11);
              
              for ii = 1:length(obj1)
                 [ind_prob(ii)] = length(find(I==ii))/MCSAMPLES;
              end
              obj.values = obj_rv;
 
          end 

          
          %% 
          
          function [obj,ind] = min(obj1)
              
              % min(obj1): Return the uncertainty object that has the
              % minimum nominal value of all the elements of obj1
              
              obj= unc_t;
              obj1_value = gmv(obj1);
              
              [~,ind]= min(obj1_value);
              obj.values = obj1(ind).values;        
              
          end
          
          function [obj,ind_prob] = min_RV(obj1)
             %% Ind_prob is an array containing the probability for each element of obj1 be the minimum.
             % Obj is the object with the highest probability to be the minimum.
              obj = unc_t;
              ind_prob = zeros(1,length(obj1));
                             
              global MCSAMPLES
              obj11 = zeros(length(obj1),MCSAMPLES);
  
              for i = 1:length(obj1)
                 [obj11(i,:)] = obj1(i).values;
              end
        
              [obj_rv, I] = min(obj11);
              
              
              for ii = 1:length(obj1)
                 [ind_prob(ii)] = length(find(I==ii))/MCSAMPLES;
              end
                 
              
              obj.values = obj_rv;
 
          end 
          
          
          %%        
          
          function obj=cumsum(obj1)
              
              % cumsum(obj1): Return the cumulative sum of the elements of obj1  

              obj=unc_t(size(obj1)); % Preallocation. obj must be of the same size as obj1.
              obj11 = zeros(size(obj1));
           
              global MCSAMPLES ;
              for kk =1:MCSAMPLES
               
                   for i1 = 1:size(obj1,1)
                       for j1 = 1:size(obj1,2)
                           obj11(i1,j1) = obj1(i1,j1).values(kk);
                       end 
                   end
                   
                   obj_result = cumsum(obj11);

                  for i = 1: size(obj1,1) 
                      for j = 1: size(obj1,2) 
                        obj(i,j).values(kk) = obj_result(i,j);
                      end 
                  end 
              
             end 
              
              
              
          end 
          
                    
          %% 
          
          function obj=sum(obj1,arg2)
              
              % sum(obj1): Sum of elements of obj1
              obj11 = zeros(size(obj1));
              obj = unc_t;
              global MCSAMPLES
              if isvector(obj1) && nargin == 1
                  for kk = 1:MCSAMPLES
                      
                      for i = 1:length(obj1)
                          obj11(i) = obj1(i).values(kk);
                      end
                      
                      obj.values(kk)= sum(obj11); 
                  end 
              elseif ismatrix(obj1) && nargin == 1
                  for kk = 1:MCSAMPLES
                      
                      for i = 1:size(obj1,1)
                          for j = 1:size(obj1,2)
                            obj11(i,j) = obj1(i,j).values(kk);
                          end
                      end
                      
                      obj_vector = sum(obj11);
                    
                      for ii = 1: size(obj1,2)
                          obj(ii).values(kk) = obj_vector(ii);
                      end
                      
                  end 
              elseif ismatrix(obj1) && nargin == 2 && strcmp(arg2,'all')
                  for kk = 1:MCSAMPLES
                      
                      for i = 1:size(obj1,1)
                          for j = 1:size(obj1,2)
                            obj11(i,j) = obj1(i,j).values(kk);
                          end
                      end
            
                          obj.values(kk) = sum(obj11,'all');
                  end 
        
              end 
       
          end 
 
          
          %% 
          
           function obj = power(obj1,obj2)
               
               % Overloading of Element-wise powers (.^)
              
              obj1=unc_t(obj1);
              obj2=unc_t(obj2);
              obj11 = zeros(size(obj1));
              obj22 = zeros(size(obj2));
              obj = unc_t;
              global MCSAMPLES
              for kk= 1: MCSAMPLES
                  
                      for i1 = 1:size(obj1,1)
                          for j1 = 1:size(obj1,2)
                            obj11(i1,j1) = obj1(i1,j1).values(kk);
                          end
                      end
                      
                      for i2 = 1:size(obj2,1)
                          for j2 = 1:size(obj2,2)
                            obj22(i2,j2) = obj2(i2,j2).values(kk);
                          end
                      end
                      
                      obj_value = obj11.^obj22;
                      
                      for i = 1:size(obj_value,1)
                          for j = 1:size(obj_value,2)
                            obj(i,j).values(kk) = obj_value(i,j);
                          end
                      end
               end 
 
                             
          
           end
          
          
          %% 
          
          function obj = mpower(obj1,obj2)
              
              % This method implements the overloading of matrix power called for the syntax 'obj1 ^ obj2'
              obj1=unc_t(obj1);
              obj2=unc_t(obj2);
              obj11 = zeros(size(obj1));
              obj22 = zeros(size(obj2));
              obj = unc_t;
              global MCSAMPLES
              for kk= 1: MCSAMPLES
                  
                      for i1 = 1:size(obj1,1)
                          for j1 = 1:size(obj1,2)
                            obj11(i1,j1) = obj1(i1,j1).values(kk);
                          end
                      end
                      
                      for i2 = 1:size(obj2,1)
                          for j2 = 1:size(obj2,2)
                            obj22(i2,j2) = obj2(i2,j2).values(kk);
                          end
                      end
                      
                      obj_value = obj11^obj22;
                      
                      for i = 1:size(obj_value,1)
                          for j = 1:size(obj_value,2)
                            obj(i,j).values(kk) = obj_value(i,j);
                          end
                      end
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
              obj1=unc_t(obj1);
              obj2=unc_t(obj2);
              obj = unc_t;
              
              obj11 = zeros(size(obj1));
              obj22 = zeros(size(obj2));
    
              global MCSAMPLES
              for kk= 1: MCSAMPLES
                  
                      for i1 = 1:size(obj1,1)
                          for j1 = 1:size(obj1,2)
                            obj11(i1,j1) = obj1(i1,j1).values(kk);
                          end
                      end
                      
                      for i2 = 1:size(obj2,1)
                          for j2 = 1:size(obj2,2)
                            obj22(i2,j2) = obj2(i2,j2).values(kk);
                          end
                      end
                      
                      obj_value = realpow(obj11,obj22);
                      
                      for i = 1:size(obj_value,1)
                          for j = 1:size(obj_value,2)
                            obj(i,j).values(kk) = obj_value(i,j);
                          end
                      end
               end 
              
          end
                   
                 
          
          %% 
          
           function obj = gt(obj1,obj2)
           
              % This method implements the overloading of the greater 
              % than operator called for the syntax 'obj1 > obj2'
              %%% only take into account comparison between scalars         
           obj=unc_t;
           
               if isa(obj1,'unc_t') && isa(obj2,'unc_t')
                   obj1=unc_t(obj1);
                   obj2=unc_t(obj2);

                   obj.values = obj1.values > obj2.values;

                   % Comparison with a constant
               elseif (isa(obj1,'unc_t') && isa(obj2,'double'))  

                    obj.values = obj1.values > obj2;

               elseif (isa(obj1,'double') && isa(obj2,'unc_t'))

                     obj.values = obj1.values > obj2.values;


               end   
           
           end    

        
           %% obj1 less than obj2
          
           function obj = lt(obj1,obj2)
               
%                This method implements the overloading of the less than operator 
%                called for the syntax 'obj1 < obj2'
               
             obj=unc_t;
           
             if isa(obj1,'unc_t') && isa(obj2,'unc_t')
               obj1=unc_t(obj1);
               obj2=unc_t(obj2);
               
                obj.values = obj1.values < obj2.values;
                
                     % Comparison with a constant
             elseif (isa(obj1,'unc_t') && isa(obj2,'double'))  

                    obj.values = obj1.values < obj2;

             elseif (isa(obj1,'double') && isa(obj2,'unc_t'))

                     obj.values = obj1.values < obj2.values;
                
                                            
             end  
           end

          
       
                                                               
          %%
          function [obj,ind]= sort(obj1) 
             %% Sort based on the mean value of the array obj1. 
             %% obj is a sorted array in ascending order.
             
%             %%%overloading sort function using bubble sort algorithm. 
%              flag_sort = 0;
%              k = 0 ;
%               while ~flag_sort
%                   flag_sort =1;
%                   k = k+1;
%                   for i = 1:length(obj1)-k
%                       if mean(obj1(i).values)> mean(obj1(i+1).values)
%                           temp = obj1(i);
%                           obj1(i) = obj1(i+1);
%                           obj1(i+1) = temp ;
%                           flag_sort = 0;
%                       end
%                   end
%               end 
%               
%               obj = obj1;
              
              %%% overloading the sort function using matlab build-in sort (faster for large arrays) 
              obj1_mean = zeros(size(obj1));
              
              for i = 1:length(obj1)
                  obj1_mean(i) = mean(obj1(i).values);
              end 
              [~,ind] = sort(obj1_mean);
                           
              obj = obj1(ind);
          
          end 
          
          function [obj, Ind_prob] = sort_RV(obj1)
              
              %% Ind_prob is a matrix containing the probability for each element. The first row contain the probability of each element of the array to be the first one. 
              % Ex: Ind_prob(1,2)-- probability that the second element of the
              % array is the first one in the sorted array.
              obj = unc_t;
              Ind_prob = zeros(length(obj1));
               
              global MCSAMPLES
              obj11 = zeros(length(obj1),MCSAMPLES);
                            
              for i = 1:length(obj1)
                 [obj11(i,:)] = obj1(i).values;
              end
             
              [obj_rv , I] = sort(obj11);
             
                 for i0 = 1:length(obj1)
                     for j0 = 1 : length(obj1)
                        obj(i0).values = obj_rv(i0,:);
                        Ind_prob(i0,j0) = length(find(I(i0,:)==j0))/MCSAMPLES;
                     end
                 end 
              
          end
          
          %% Get Covariance and correlation matrices 
          function CX_cov = get_cov_mat(X)

           % Return the covariance matrix CX_cov of the uncertainty 
           % objects contained in the array X. 
           global MCSAMPLES

           if size(X,1)==1
               
                  A = zeros(MCSAMPLES,length(X));
                  % matrix whose columns represent random variables and whose rows represent observations
                  for i=1:length(X)
                      A(:,i) = X(i).values;  
                  end

                  CX_cov = cov(A);              
           else 
            error('Input argument must be a one-row array of uncertainty objects!')
           end 
           
          end 
          
          function CX_cor = get_cor_mat(X)

           % Return the covariance matrix CX_cov of the uncertainty 
           % objects contained in the array X. 
           global MCSAMPLES

           if size(X,1)==1
               
                  A = zeros(MCSAMPLES,length(X));
                  % matrix whose columns represent random variables and whose rows represent observations
                  for i=1:length(X)
                      A(:,i) = X(i).values;  
                  end

                  CX_cor = corr(A);              
           else 
            error('Input argument must be a one-row array of uncertainty objects!')
           end 
           
          end 
          
                  
          
          
          %%
             function obj = fft(obj1)
             % Overloading the Discrete Fourier transform 
             
              obj1=unc_t(obj1);
              obj = unc_t;
              
              obj11 = zeros(size(obj1));
                  
              global MCSAMPLES
              for kk= 1: MCSAMPLES
                  
                      for i1 = 1:size(obj1,1)
                          for j1 = 1:size(obj1,2)
                            obj11(i1,j1) = obj1(i1,j1).values(kk);
                          end
                      end
                      
                      obj_value = fft(obj11);
                      
                      for i = 1:size(obj_value,1)
                          for j = 1:size(obj_value,2)
                            obj(i,j).values(kk) = obj_value(i,j);
                          end
                      end
               end 
              
          end 
          
          function obj = myfft(obj1)
                %only works if N = 2^k
                N = numel(obj1);
                xp = obj1(1:2:end);
                xpp = obj1(2:2:end);
                if N>=8
                     Xp = myfft(xp);
                     Xpp = myfft(xpp);
                     obj = zeros(N,1);
                     Wn = exp(-1i*2*pi*((0:N/2-1)')/N);
                     try
                         tmp = Wn .* Xpp;
                     catch
                         i=1;
                     end
                     obj = [(Xp + tmp);(Xp - tmp)];
                else
                     switch N
                         case 2
                             obj = [1 1;1 -1]*obj1;
                         case 4
                             obj = [1 0 1 0; 0 1 0 -1i; 1 0 -1 0;0 1 0 1i]*[1 0 1 0;1 0 -1 0;0 1 0 1;0 1 0 -1]*obj1;
                         otherwise
                             error('N not correct.');
                     end
                end
            end
          
          
          function obj = ifft(obj1)
           % Overloading the inverse discrete Fourier transform of obj1
           
              obj1=unc_t(obj1);
              obj = unc_t;
              
              obj11 = zeros(size(obj1));
                  
              global MCSAMPLES
              for kk= 1: MCSAMPLES
                  
                      for i1 = 1:size(obj1,1)
                          for j1 = 1:size(obj1,2)
                            obj11(i1,j1) = obj1(i1,j1).values(kk);
                          end
                      end
                      
                      obj_value = ifft(obj11);
                      
                      for i = 1:size(obj_value,1)
                          for j = 1:size(obj_value,2)
                            obj(i,j).values(kk) = obj_value(i,j);
                          end
                      end
               end 
              
          end 
          
          %% 1D Linear interpolation
        function Yp = interp1(X,Y,Xp)
            % overloading 1D linear interpolation 
            X = unc_t(X);
            Y = unc_t(Y);
            Xp = unc_t(Xp);
            Yp = unc_t(size(Xp));
            
            X11 = zeros(size(X));
            Y11 = zeros(size(Y));
            Xp11 = zeros(size(Xp));
            
            global MCSAMPLES
              for kk = 1: MCSAMPLES
                  
                      for i1 = 1:length(X)
                          X11(i1) = X(i1).values(kk);
                          Y11(i1) = Y(i1).values(kk);
                      end
                                            
                      for i2 = 1:length(Xp)
                          Xp11(i2) = Xp(i2).values(kk);
                      end
                      
                      Yp_value = interp1(X11,Y11,Xp11);
                      
                      for i = 1:length(Yp_value)
                             Yp(i).values(kk) = Yp_value(i);
                          
                      end
               end 

        end


           
          
%% Calibration. This method implements a general procedure for polynomial calibration, computing the fitting parameters when using polynomial regression algorithm.

    function [Y_cal, cali_param] = calibrate(Y, x_cal, y_cal, n) 
        % Y : array of sensor signal
        % x_cal : column vector of calibration points 
        % y_cal : array with the output of sensors at the selected
        % calibration points.
        % n : order of polynomial calibration(starting with 0)
        % n=0 --> offset calibration
        % n=1 --> linear gain 
        % n=2 --> quadratic polynomial compensation 
        
        numb_cal = length(x_cal); % number of calibration points
        H_cal = unc_t(zeros(numb_cal,n+1)); % information matrix 
        Y_cal = unc_t;
        h = unc_t(zeros(length(Y),n+1));
        if n > numb_cal-1
            
            error('The number of calibration points has to exceed the polinomial order by at least one.');
        else 
            
            % computing the information matrix H_cal
            for j = 1:n+1
                for i = 1: numb_cal
                   H_cal(i,j) = y_cal(i).^(j-1) ;
                   h(:,j) = Y.^(j-1);
                end
            end
            % vector of calibration parameters
            cali_param  = inv(H_cal'*H_cal)*H_cal'* x_cal ;
            
            % calibration calculation
            for ii = 1: length(Y) 
                
                Y_cal(ii) = h(ii,:) * cali_param ;
            end 
            
         end
        
    end         
%% Copula Plot 

    function copula_plot(obj,names)
      % obj is the array of inputs elements 
      % names is an array with the names of the elements in x.

      lobj=length(obj); % length of the input vector  

      global MCSAMPLES
      
      figure

          for ii=1:lobj

              [xs,i]=sort(obj(ii).values);

              for jj=1:lobj

                  [~,j]=sort(obj(jj).values);

                  v = obj(1).values; 
                      for k=1:MCSAMPLES 
                          v(k)=find(j(k)==i);
                      end

                  subplot(lobj+1,lobj+1,(ii-1)*(lobj+1)+jj)
                  o=histogram2((1:MCSAMPLES)/MCSAMPLES,v/MCSAMPLES,'BinMethod','scott');
                  hold on;
                  plot3((1:MCSAMPLES)/MCSAMPLES,v/MCSAMPLES,0.*v+1.1*max(o.Values,[],'all'),'r.') 
                  view(2)

                  if (jj==1)&&(nargin>1)
                      ylabel(names{ii});
                  end


              end

              subplot(lobj+1,lobj+1,(lobj+1)*ii)
              hold on
              histogram(xs,'Normalization','pdf')
              plot(xs,(1:MCSAMPLES)/MCSAMPLES,'r')


              subplot(lobj+1,lobj+1,lobj*lobj+lobj+ii)
              histogram(-xs,'Normalization','pdf','orientation','horizontal')
              hold on
              plot((1:MCSAMPLES)/MCSAMPLES,-xs)
              yticklabels(-yticks)

              if (nargin>1)
                  xlabel(names{ii});
              end



          end



    end
    
%% Correlation Plot

function corr_plot(varargin)
      % correlation plot for unc_t objects
      % refer to help corrplot for more details

      global MCSAMPLES
      obj=varargin(1);
      Vec=zeros(MCSAMPLES,length(obj{1}));  
      for i=1:length(obj{1})
          Vec(:,i)=obj{1}(i).values';
      end
      varargin{1}=Vec;
      figure
      
      if (size(varargin)==1)
          corrplot(varargin{:});
      else
      corrplot(varargin{1},'varNames',varargin{2});
      end
      
      
 end
      

      
    %%
    function [obj] = norm(obj1,arg)
        %% overloading the norm method for unc_t vectors/matrices.
        
        global MCSAMPLES
        
        obj11 = zeros(size(obj1,1),size(obj1,2),MCSAMPLES);
        obj = unc_t;
       
            for i = 1:size(obj1,1)
                for j = 1:size(obj1,2)
                    obj11(i,j,:) = obj1(i,j).values;
                end
            end 
            
            if nargin == 1 
                
                for kk = 1:MCSAMPLES
                   obj.values(kk) = norm(obj11(:,:,kk));
                end
                
            else
                for kk = 1:MCSAMPLES
                   obj.values(kk) = norm(obj11(:,:,kk),arg);
                end
            end 
   end 

%% Handeling complex unc objects. 

function Z = complex(x,y)
    
    % Z = complex(x,y): Create complex uncertainty objects Z out of real uncertainty objects
    % x and y
    
   [r1,c1]=size(x); 
   
   [r2,c2]=size(y);
   
   Z=unc_t(); 
   
   if (r1==r2) && (c1==c2)  
    
       for i=1:r1
           for j=1:c1
               
               if all(imag(x(i,j).values)==0) && all(imag(y(i,j).values)==0)
       
                            
                Z(i,j).values = x(i,j).values + 1i * y(i,j).values;
                
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
     X=unc_t();
     
     for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                X(i,j).values = real(Z(i,j).values);
                
         end
     end
             
 end

 
%% 
          
 function Y=imag(Z)
     
     % imag(Z): Return the imaginary part of Z
     Y=unc_t();
             
         for i=1:size(Z,1)
           for j=1:size(Z,2)
               
               Y(i,j).values = imag(Z(i,j).values);
  
           end
         end
         
         
         
     end
     
       
%% 
          
 function obj=conj(Z)
       
     % This method takes the complex conjugate of each
     % entry of the matrix Z
     obj=unc_t();
     
       for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                 obj(i,j).values = real(Z(i,j).values) - 1i * imag(Z(i,j).values);
                
         end
       end
       
             
 end  
 
 
 
 %% 
          
 function obj=abs(Z)
   % overloading the function abs() that returns the absolute value of the uncertainty object Z. When
   % Z is complex it returns the complex modulus (magnitude) of the elements of Z.
   obj = unc_t();
      for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                 obj(i,j).values = sqrt((real(Z(i,j).values)).^2 + (imag(Z(i,j).values)).^2);
                
         end
       end
             
 end  
 
 function obj=angle(Z)
      % overloading the function angle() that returns the phase angles, in radians, of the uncertainty object Z.
      obj = unc_t();
       for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                 obj(i,j).values = atan2(imag(Z(i,j).values), real(Z(i,j).values));
                
         end
       end
       
             
 end 
 
 function Jac = get_Jac(obj1, obj2)
     % get jacobian for obj1 (vector) with respect. to obj2 (vector)
     % obj2 must comprise primary RVs
     
     lo1=length(obj1);
     lo2=length(obj2);
     Jac=zeros(lo1,lo2);
     for i=1:lo1
         for j=1:lo2
             C=cov(obj1(i).values,obj2(j).values);
             Jac(i,j)=C(1,2)/C(2,2);
         end
     end
 end
 
 function disp_contribution(obj1)
     warning('The function disp_contribution is not implemented in the Monte Carlo approach!')
 
 end
 
 
%% Displaying uncertainty objects 
      function disp(obj)
          
          objValue = zeros(size(obj));
          obj_X = zeros(size(obj));
          obj_Y = zeros(size(obj));
          objXStd_unc = zeros(size(obj));
          objYStd_unc = zeros(size(obj));
          objStd_unc = zeros(size(obj));
          global MCSAMPLES
          
              % Overloading of the Matlab disp() function
                    
             for i =1 :size(obj,1)
                for j=1 : size(obj,2)
                    
                    objValue(i,j) = mean(obj(i,j).values);
                                               
                     if islogical(obj(i,j).values)
                         if nnz(obj(i,j).values) > nnz(~obj(i,j).values)
                            objLogic(i,j) = true;
                            objProb(i,j) = nnz(obj(i,j).values)/MCSAMPLES;
                         elseif nnz(obj(i,j).values) < nnz(~obj(i,j).values)
                            objLogic(i,j) = false;
                            objProb(i,j) = nnz(~obj(i,j).values)/MCSAMPLES;
                             
                         end 
                         % objLogic represents the logical decision 
                         % objProb represents the probability of make the
                         % right desicion.
  
                         fprintf([' \t',num2str(objLogic(i,j)),'(',num2str(objProb(i,j)),')']);
                         
                     elseif imag(objValue(i,j))~=0 
                           
                            obj_X(i,j) = real(objValue(i,j));
                            obj_Y(i,j) = imag(objValue(i,j));
                            
                            % obtaining the std_unc of the real and img
                            % parts 
                            objXStd_unc(i,j) = std(real(obj(i,j).values));
                            objYStd_unc(i,j) = std(imag(obj(i,j).values));
                           
                            if obj_Y(i,j)>= 0
                           
                               fprintf(['    \t',obj.GenerateDispStr(obj_X(i,j),objXStd_unc(i,j)) ,' + '...
                                         obj.GenerateDispStr(obj_Y(i,j),objYStd_unc(i,j)) , ' * i ']);
                           
                            else
                              fprintf(['    \t',obj.GenerateDispStr(obj_X(i,j),objXStd_unc(i,j)) ,' - '...
                                        obj.GenerateDispStr(abs(obj_Y(i,j)),objYStd_unc(i,j)) , ' * i ']);
                            end
                                
                   
                       
                     elseif imag(objValue(i,j))==0
                         
                         objStd_unc(i,j) = std(obj(i,j).values);
                                  
                            
                                         fprintf([' \t',obj.GenerateDispStr(objValue(i,j),objStd_unc(i,j))]);
                               
                            
                            
                     end
                            
                end  
                
                fprintf('\n')
                
             end
             
          end
 
 
end

          %%
          methods (Static)
              %%%%%%%%%%%%%%%%%function GenerateDispStr modified%%%%%%
          function DispString = GenerateDispStr(value,std_unc_t)
               
               % This method is used by the disp() method to customize 
               % the display of uncertainty objects
               
                global Number_of_Digits_to_Display;
                d=Number_of_Digits_to_Display-1; 
                
                
                % 1) Clip the standard uncertainty to max 3 digits
                if (std_unc_t==0) || (log10(abs(value)/abs(std_unc_t))>12)
                    digit=-12;
                else
                    digit = floor(log10(abs(std_unc_t))); % take the integer part no matter is the value of the decimal part
                end
                
                
                
                disp_std_unc_t = round(std_unc_t / 10^(digit-d)); % round the value into the closest integer.
                
                
                % 2) Clip the value to show the last d+1 values 
                disp_value = round(value / 10^(digit-d))*10^(digit-d);
                
                
                if disp_std_unc_t==0
                        eval_str = '[num2str(disp_value),''(0)'']';
                        DispString = eval (eval_str);
                else
                     
                    if (digit<-4) && (abs(disp_value)<0.001) && (disp_value~=0)
                        % show the result in scientific notation if there
                        % are more than 2 leading zeros
                        eval_str = sprintf ('[sprintf(''%%0.0%ie'',disp_value),''('',num2str(disp_std_unc_t),'')'']',abs(digit)-abs(floor(log10(abs(value))))+d);
                        DispString = eval (eval_str);
                    elseif (digit<0) 
                        eval_str = sprintf ('[sprintf(''%%0.0%if'',disp_value),''('',num2str(disp_std_unc_t),'')'']',abs(digit)+d);
                        DispString = eval (eval_str);
                    else
                        if digit<d
                            disp_std_unc_t = disp_std_unc_t*10^(digit-d); %#ok<NASGU>
                             eval_str = sprintf ('[sprintf(''%%0.0%if'',disp_value),''('',sprintf(''%%0.0%if'',disp_std_unc_t),'')'']',d-digit,d-digit);
                            DispString = eval (eval_str); 
                        else
                            disp_std_unc_t = disp_std_unc_t*10^(digit-d);
                            DispString = [num2str(disp_value),'(',num2str(disp_std_unc_t),')'];
                        end
                    end
                end
          end
           
 %%
           
       end
end
