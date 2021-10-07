classdef unc < handle
  % AUT - An Uncertainty Toolbox V 0.91
  % Date: 27/11/2019
  % Type "help disclaimer" for more information

   
  
  %% Properties of the uncertainty class
   
       properties  
       % Define the properties of the unc object to be created by the Class Constructor 
          value ;  % Mean value of the unc object
          std_unc; % Standard deviation of the unc object 
          name;    % Name of the unc object
          dep;     % Array containing all the random variables (unc objects) 
                   % on which the unc object depends
          grad;    % Array containing all gradient values of the unc object
                   % The gradient is calculated with respect to the
                   % elements of the dep property
          OutDataTypeStr; % Data type of an uncertainty object
          rel_obj;  % Cell-array containing the unc objects with which a unc 
                    % class object is correlated. The correlation and covariance 
                    % values are stored in the property called rel_mat.
          rel_mat;  % 2-row matrix in which: -first row contains covariance values
                    %                         of the unc object with respect to other
                    %                         random variables (unc objects) whose names
                    %                         are found in the property rel_mat
                    %                        -second row contains the correlation values
                    %                         of the unc object with respect to other
                    %                         random variables (unc objects) whose names
                    %                         are found in the property rel_mat
          r;   % Real part of the unc object
          img; % Imaginary part of the unc object
          
       end
       
                   
  %% Methods of the uncertainty class    
       
       methods

  %% Proprties handling methods

       function obj = set.name(obj,Name)                    
                           
                 if ~isempty(Name)

                             [obj_names,obj_array] = obj.get_ws_obj();
                             
                             if isa(Name,'char')
                                 
                                 
                                 
                                 % ------------ Check for same variable names ------------ %
                                 
                                     c=0; % counter for number of times the variable name 'Name' is assigned
                                     
                                     for i=1:length(obj_names)   
                                         
                                       temp = obj_array{i}.name ;
                                       
                                       if strcmp(temp,Name)  
                                           c=c+1;
                                       end                                  
                                     end 
                                     
                                     if c>=1
                                                warning(' Variable name ''%s'' is already ''%d'' time(s) assigned!'...
                                                           ,Name,c);
                                     end
                                    
                                      
                                     % ------------ End of check for same variable names ------------ %
                                 
                                     obj.name = {Name}; % convert from char into cell
                                     
                                     
                             else

                                 if isa(Name,'cell')
                                     
                                  if length(Name)==1   
                                     
                                    % ------------ Check for same variable names ------------ %
                                    
                                     c=0; % counter for number of times the variable name 'Name' is assigned
                                     
                                     for i=1:length(obj_names)   
                                         
                                       temp = obj_array{i}.name ;
                                       if strcmp(temp,Name)  
                                           c=c+1;
                                       end                                  
                                     end 
                                     
                                      if c==1
                                          warning(' Variable name ''%s'' already exists in Workspace or assigned to another unc object!',Name{1}); 
                                      else
                                          if c>1
                                                warning(' Variable name ''%s'' is already ''%d'' time(s) assigned to (an)other variable(s)!'...
                                                           ,Name{1},c);
                                          end
                                      end
                                      
                                     % ------------ End of check for same variable names ------------ %
                                     
                                     obj.name = Name;
                                  else
                                      error('One name only is allowed!')
                                  end
                                     
                                 else
                                    error('Invalid datatype! Name must be character or cell data type.')
                                 end
                             end
                 end
                                           
       end
                 
       
       function obj = set.value(obj,Value)
           
                        if imag(Value)~=0
                            
%                             obj_X=real(obj);
%                             obj_Y=imag(obj);
                            
                            Value_X=real(Value);
                            Value_Y=imag(Value);
                            
                                                   
                            if isa(Value,'double')
                               
                                    obj.value = Value;
                                    
                            else
                                
                                error('Value property assignment to complex unc objects is not needed.')
                                             
                            end
                            
                        else
                            if imag(Value)==0
           
%                                     if ~isempty(Value)
                                        if isa(Value,'double')
                                            obj.value = Value;
                                        else
                                            error('Invalid datatype!Value must be a double!')
                                        end
%                                     end
                            end
                            
                        end
                        
                        
       end
                                   
          
       function obj = set.std_unc(obj,Std_Unc)
           
                    if imag(obj)~=0
                            
                            obj_X=real(obj);
                            obj_Y=imag(obj);
                            
                            if isa(Std_Unc,'double') 
                                
                                    obj.std_unc = Std_Unc;
                                    
                            else isa(Std_Unc,'double') %&& ( Std_Unc == sqrt( obj_X.std_unc^2 + obj_Y.std_unc^2))
                                
                                error('Standard uncertainty property assignment to complex unc objects is not needed.')
                                             
                            end
                            
                        else
                            if imag(obj)==0
                                
                                    if ~isempty(Std_Unc)
                                        if isa(Std_Unc,'double')
                                            obj.std_unc = Std_Unc;   
                                        else
                                            error('Invalid datatype!std_unc must be a double!')
                                        end
                                    end
                            end
                    end
       end
                                                
           
           
  %% ------ Class constructor ------- %
  
       function obj = unc(arg1,arg2,arg3)
                           
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
                           
                       
          % Case(1): Number of input arguments is equal to 3:
          
          if (nargin==3) %  
                            
             % Syntax (1) : unc(Value,Uncertainty,'name') where the property name is of type character
                           
             if isa(arg1,'double')&& isa(arg2,'double')&& isa(arg3,'char')&&...
                     length(arg1)==1 && length(arg2)==1 && size(arg3,1)==1 && isreal(arg1)
                            
                               obj.value=arg1;
                               obj.std_unc=arg2;
                               obj.name={arg3}; % convert from char type to cell type 
                               obj.dep=obj;
                               obj.grad=1;
                               obj.OutDataTypeStr='double';
                               obj.rel_obj={};
                               obj.rel_mat=[]; 
                               obj.img=0;
                               obj.r= obj; 
                               
             else              
                 % Syntax (2): unc(Value, Uncertainty, Name) where the property name is of type cell
                
                    if isa(arg1,'double')&& isa(arg2,'double')&& isa(arg3,'cell') && size(arg1,1)==size(arg2,1) ...
                            && size(arg1,2)==size(arg2,2) && size(arg1,1)==size(arg3,1) && size(arg1,2)==size(arg3,2)...
                            && all(all(isreal(arg1)))  
                                
                              obj(size(arg1,1),size(arg1,2))=unc;
         
                              for i=1:size(arg1,1)
                                  for j=1:size(arg1,2)
                                      obj(i,j).value = arg1(i,j);
                                      obj(i,j).std_unc = arg2(i,j);
                                      obj(i,j).name = arg3(i,j);
                                      obj(i,j).dep=obj(i,j);
                                      obj(i,j).grad=1;
                                      obj(i,j).OutDataTypeStr='double';
                                      obj(i,j).rel_obj={};
                                      obj(i,j).rel_mat=[]; 
                                      obj(i,j).img=0;
                                      obj(i,j).r= obj(i,j);
                                  end
                              end
                    else
                        
                    % Syntax (2.1): unc(Value, Uncertainty, Name) where the property name is of type character. 
                    % In the case it is not interested to assign a different name to each array's element.
                    % All the elements of the array are going to have the
                    % same name                                        
                
                        if isa(arg1,'double')&& isa(arg2,'double')&& isa(arg3,'char') && size(arg1,1)==size(arg2,1) ...
                            && size(arg1,2)==size(arg2,2) && size(arg3,1)== 1 && all(all(isreal(arg1)))
                                
                              obj(size(arg1,1),size(arg1,2))=unc;
         
                              for i=1:size(arg1,1)
                                  for j=1:size(arg1,2)
                                      obj(i,j).value = arg1(i,j);
                                      obj(i,j).std_unc = arg2(i,j);
                                      obj(i,j).name = arg3 ;
                                      obj(i,j).dep=obj(i,j);
                                      obj(i,j).grad=1;
                                      obj(i,j).OutDataTypeStr='double';
                                      obj(i,j).rel_obj={};
                                      obj(i,j).rel_mat=[]; 
                                      obj(i,j).img=0;
                                      obj(i,j).r= obj(i,j); 
                                  end
                              end
                    else

                    % Syntax (3) and (4): (array of unc objects, input covariance matrix, 'cov') or 
                    %                     (array of unc objects, input correlation matrix, 'corr')
                    
                       if isa(arg1,'unc') && size(arg1,1)==1 && isa(arg2, 'double') &&  size(arg2,1)==size(arg1,2)...
                            && size(arg2,2)==size(arg1,2) && isa(arg3,'char') 
                        
                        %check if the input matrix is symmetric and PSD  
                        [~,p] = chol(arg2);
                        if (isequal(arg2,transpose(arg2)) && (p==0))
                                                                                                  
                                    if strcmp(arg3,'cov') % If user enters a covariance matrix                                        
                                        
                                      obj = arg1; 
                                      
                                      % %+++ Covariances and correlations assignments ++++
                                      n=size(arg1,2); % number of unc class objects
                                      for i=1:n
                                          k=1; % counter
                                          for j=1:n
                                              if i~=j
                                                   %check if the values of covariances
                                                   %are inside the range.
                                                   if ( arg2(i,j) > (arg1(i).std_unc * arg1(j).std_unc) ) || ( arg2(i,j) < -(arg1(i).std_unc * arg1(j).std_unc) )
                                                       error('Incorrect value in entry(%d,%d) of the covariance matrix! \n The covariance between two random variables x and y must satisfy the property |cov(x,y)|%s %sx*%sy, where %sx and %sy represent the standard deviation of x and y respectively.'...
                                                           ,i,j,num2str(char(8804)),num2str(char(963)),num2str(char(963)),num2str(char(963)),num2str(char(963)));
                                                   else
                                                      arg1(i).rel_obj{1,k}= arg1(j);
                                                      arg1(i).rel_mat(1,k)=arg2(i,j); % get covariance
                                                      arg1(i).rel_mat(2,k)=arg2(i,j)/ sqrt(arg2(i,i) * arg2(j,j)); % get correlation
                                                      k=k+1; % increment counter
                                                      
                                                   end
                                                                                                     
                                              end
                                          end
                                      end
                                      % perform the karhunen-loeve decomposition. update the dep property 
                                       [~, ~] = karhloev( arg1,arg2,'cov' );

                                     % % End of assignments % %
                             
          
                                     else
                                         if strcmp(arg3,'corr') % User enters a correlation matrix
                                             
                                             obj = arg1; 
                                             
                                             % %+++ Covariances and correlations assignments ++++
                                             n=size(arg1,2); % number of unc class objects.
                                             for i=1:n
                                                  k=1; % counter
                                                  for j=1:n
                                                      if i~=j
                                                          %check if the values of correlation
                                                          %are inside the range.
                                                          if ( arg2(i,j) > 1 ) || ( arg2(i,j) < -1 )
                                                          error('Incorrect value in entry(%d,%d) of the correlation matrix! \n The correlation coefficient between two random variables x and y must satisfy the property |corr(x,y)|%s 1.'...
                                                           ,i,j,num2str(char(8804)));
                                                          else
                                                          
                                                              arg1(i).rel_obj{1,k}= arg1(j); 
                                                              arg1(i).rel_mat(2,k)=arg2(i,j); % get correlation.
                                                              arg1(i).rel_mat(1,k)=arg2(i,j) * arg1(i).std_unc * arg1(j).std_unc ; % get covariance.
                                                              k=k+1; % increment counter 
                                                          end
                                                        
                                                      end
                                                  end
                                             end
                                       % perform the karhunen-loeve decomposition. update the dep property 
                                       [~, ~] = karhloev( arg1,arg2,'corr' );
                                             % % End of assignments % %
                                            
   
                                         else
                                                error('Incorrect type of argument: use cov to enter a covariance matrix, or corr to enter a correlation matrix.');
                                         end
                                    end
                        else
                            error('The covariance and correlation matrices must be symmetric and positive semidefinite.');
                        end 
                   
                                else
                                        error('Incorrect type of argument(s)');
                        
                        end
                       
                       end
                     end
               end
          end
          
          
          % Case(2): Number of input arguments is equal to 2
               
          if (nargin==2) 
                          
             % Syntax (1): unc(value,'name') where the property name is of type character
                           
             if isa(arg1,'double')&& isa(arg2,'char')&& length(arg1)==1 && size(arg2,1)==1 && all(all(isreal(arg1))) 
                            
                               obj.value=arg1;
                               obj.std_unc=0;
                               obj.name=arg2;
                               obj.dep=obj;
                               obj.grad=1;
                               obj.OutDataTypeStr='double';
                               obj.img=0;
                               obj.r= obj; 
                               
             elseif isa(arg1,'double')&& isa(arg2,'char')&& length(arg1)==1 && size(arg2,1)==1 && any(any(~isreal(arg1)))
                        
                            obj = complex(unc(real(arg1)), unc(imag(arg1)));
                            obj.name = arg2; 
                               
                 
             else
                 
                    
                 % Syntax (2):  unc(Value, Uncertainty)
                 
                 if isa(arg1,'double')&& isa(arg2, 'double') && size(arg1,1)==size(arg2,1) && size(arg1,2)==size(arg2,2) && all(all(isreal(arg1))) 
                      obj(size(arg1,1),size(arg1,2))=unc;
                      for i=1:size(arg1,1)
                          for j=1:size(arg1,2) 
                              obj(i,j).value   = arg1(i,j);
                              obj(i,j).std_unc = arg2(i,j);
                              obj(i,j).name    = {};
                              obj(i,j).dep=obj(i,j);
                              obj(i,j).grad=1;
                              obj(i,j).OutDataTypeStr='double';
                              obj(i,j).rel_obj={};
                              obj(i,j).rel_mat=[];
                              obj(i,j).img=0;
                              obj(i,j).r= obj(i,j); 
                          end
                      end
                      
                 elseif isa(arg1,'double')&& isa(arg2, 'double') && size(arg1,1)==size(arg2,1) && size(arg1,2)==size(arg2,2) &&...
                         any(any(~isreal(arg1))) && all(all(arg2 == 0))
                     
                     obj(size(arg1,1),size(arg1,2))=unc;
                        for i=1:size(arg1,1)
                          for j=1:size(arg1,2)
                              
                             obj(i,j) = complex(unc(real(arg1(i,j))),unc(imag(arg1(i,j))));  
                           
                           end
                        end
                        
                        
                 elseif isa(arg1,'double')&& isa(arg2, 'double') && size(arg1,1)==size(arg2,1) && size(arg1,2)==size(arg2,2) &&...
                         any(any(~isreal(arg1))) && any(any(arg2 ~= 0))
                     
                     error('Complex uncertainty objects should be definined in terms of its real and imaginary components. Use function Complex.' ); 
                    
                 else
                     
                    % Syntax (3): unc(Value, Name)
                    
                    if isa(arg1,'double')&& isa(arg2, 'cell') && size(arg1,1)==size(arg2,1) && size(arg1,2)==size(arg2,2) && all(all(isreal(arg1))) 
                              obj(size(arg1,1),size(arg1,2))=unc;
                              for i=1:size(arg1,1)
                                  for j=1:size(arg1,2)
                                      obj(i,j).value   = arg1(i,j);
                                      obj(i,j).std_unc = 0;
                                      obj(i,j).name    = arg2(i,j);
                                      obj(i,j).dep=obj(i,j);
                                      obj(i,j).grad=1;
                                      obj(i,j).OutDataTypeStr='double';
                                      obj(i,j).rel_obj={};
                                      obj(i,j).rel_mat=[];
                                      obj(i,j).img=0;
                                      obj(i,j).r= obj(i,j); 
                                  end
                              end
                    elseif isa(arg1,'double')&& isa(arg2, 'cell') && size(arg1,1)==size(arg2,1) && size(arg1,2)==size(arg2,2) && any(any(~isreal(arg1))) 
                        
                        obj(size(arg1,1),size(arg1,2))=unc;
                        for i=1:size(arg1,1)
                          for j=1:size(arg1,2)
                              
                             obj(i,j) = complex(unc(real(arg1(i,j))),unc(imag(arg1(i,j))));  
                           
                           end
                       end                     
                        
                       

                       else
                            
                            
                       % Syntax (4): unc(matrix of values , input covariance matrix)
                       
                       if isa(arg1,'double')&& isa(arg2, 'double') && size(arg2,1)==size(arg2,2) && size(arg2,1)==size(arg1,1)*size(arg1,2) && all(all(isreal(arg1)))
                           
                           %check if the cov matrix is symmetric and PSD       
                           [~,p] = chol(arg2);
                           if (isequal(arg2,transpose(arg2)) && (p==0))
                               
                              obj(1, size(arg1,1) * size(arg1,2))=unc; % Initialization: create a row vector of unc objects
                              
                              % Create row vector arg1_row from matrix arg1
                                      arg1_row=arg1;
                                      arg1_row=arg1_row.';
                                      arg1_row=arg1_row(:);
                                      %
                              
                              %
%                               c=1; % counter
                              diag_arg2 = diag(arg2); % array containing the elements on the diagonal of arg2
                              %
                              
                              for i=1:size(arg1,1) * size(arg1,2)
                                      obj(i).value = arg1_row(i);
                                      %
                                      % get std_unc from the diagonal of arg2
                                      obj(i).std_unc = sqrt(diag_arg2(i));
%                                       c=c+1;
                                      %
                                      obj(i).dep=obj(i);
                                      obj(i).grad=1;
                                      obj(i).OutDataTypeStr='double';
                                      obj(i).img=0;
                                      obj(i).r= obj(i); 
                              end   
                              
                              
                              % %+++ Covariances and correlations assignments ++++ 
                              
                                      n=size(arg2,1); % number of unc class objects
                                      for i=1:n
                                          k=1; % counter
                                          for j=1:n
                                              if i~=j
                                                   %check if the values of covariances
                                                   %are inside the range.
                                                   if ( arg2(i,j) > (obj(i).std_unc * obj(j).std_unc) ) || ( arg2(i,j) < -(obj(i).std_unc * obj(j).std_unc) )
                                                       error('Incorrect value in entry(%d,%d) of the covariance matrix! \n The covariance between two random variables x and y must satisfy the property |cov(x,y)|%s %sx*%sy, where %sx and %sy represent the standard deviation of x and y respectively.'...
                                                           ,i,j,num2str(char(8804)),num2str(char(963)),num2str(char(963)),num2str(char(963)),num2str(char(963)));
                                                   else
                                                      obj(i).rel_obj{1,k}= obj(j);
                                                      obj(i).rel_mat(1,k)=arg2(i,j); % get covariance
                                                      obj(i).rel_mat(2,k)=arg2(i,j)/ sqrt(arg2(i,i) * arg2(j,j)); % get correlation
                                                      k=k+1; % increment counter
                                                   end
                                              end
                                          end
                                      end
                                      
                                      % perform the karhunen-loeve decomposition. update the dep property 
                                       [~, ~] = karhloev( obj,arg2,'cov' );
                                                                           
                                % % End of assignments % %
                           else
                               error('Covariance matrix must be symmetric and positive semidefinite (PSD).');
                           end 
                              
                          else
                           
                          % Syntax (5): unc(array of unc objects, input covariance matrix)
                          
                          if isa(arg1,'unc')&& isa(arg2, 'double') &&  size(arg2,1)==size(arg1,2) && size(arg2,2)==size(arg1,2)  
                             
                                      obj = arg1;
                                %check if the cov matrix is symmetric and PSD       
                                 [~,p] = chol(arg2);
                                if (isequal(arg2,transpose(arg2)) && (p==0))
                                      % %+++ Covariances and correlations assignments ++++
                                      n=size(arg1,2); % number of unc class objects
                                      for i=1:n
                                          k=1; % counter
                                          for j=1:n
                                              if i~=j
                                                  %check if the values of covariances
                                                   %are inside the range.
                                                   if ( arg2(i,j) > (arg1(i).std_unc * arg1(j).std_unc) ) || ( arg2(i,j) < -(arg1(i).std_unc * arg1(j).std_unc) )
                                                       error('Incorrect value in entry(%d,%d) of the covariance matrix! \n The covariance between two random variables x and y must satisfy the property |cov(x,y)|%s %sx*%sy, where %sx*%sy represent the standard deviation of x and y respectively.'...
                                                           ,i,j,num2str(char(8804)),num2str(char(963)),num2str(char(963)),num2str(char(963)),num2str(char(963)));
                                                   else
                                                      arg1(i).rel_obj{k}= arg1(j); 
                                                      arg1(i).rel_mat(1,k)=arg2(i,j); % get covariance
                                                      arg1(i).rel_mat(2,k)=arg2(i,j)/ sqrt(arg2(i,i) * arg2(j,j)); % get correlation
                                                      k=k+1; % increment counter
                                                   end
                                              end
                                          end
                                      end
                                       % perform the karhunen-loeve decomposition. update the dep property 
                                       [~, ~] = karhloev( obj,arg2,'cov' );
                                     % % End of assignments % %
                               else
                                   error('Covariance matrix must be symmetric and positive semidefinite');
                               end
                          else
                              
                                     error('Incorrect type of argument(s)');    
                          end
                       end
                    end
                 end
              end
          end
              
             
          % Case(3): Number of input arguments is equal to 1
                
          if (nargin==1)    
              
             % Syntax(1): Value only
             
             if isa(arg1,'double') && all(all(isreal(arg1)))
                 
                       obj(size(arg1,1),size(arg1,2))=unc;
                       for i=1:size(arg1,1)
                          for j=1:size(arg1,2)
                              obj(i,j).value   = arg1(i,j);
                              obj(i,j).std_unc = 0;
                              obj(i,j).name    = {};
                              obj(i,j).dep=[];
                              obj(i,j).grad=[];
                              obj(i,j).OutDataTypeStr='double';
                              obj(i,j).rel_obj={};
                              obj(i,j).rel_mat=[];
                              obj(i,j).img=0;
                              obj(i,j).r= obj(i,j); 
                          end
                       end
                       
             elseif isa(arg1,'double') && any(any(~isreal(arg1)))                 
                 
                        obj(size(arg1,1),size(arg1,2))=unc;
                        for i=1:size(arg1,1)
                          for j=1:size(arg1,2)
                              
                             obj(i,j) = complex(unc(real(arg1(i,j))),unc(imag(arg1(i,j))));  
                           
                           end
                       end 
             else
                 
                 % Syntax(2): unc object only
                 
                 if isa(arg1,'unc') 
                       
                       obj = arg1;
                       
                 else
                           error('Incorrect type of argument!');
                 end
             end
             
          end
          
               
          % Case(4): No input arguments   
                
          if (nargin==0)  
                          obj.value   = 0;
                          obj.std_unc = 0;
                          obj.name    = {}; 
                          obj.dep=obj;
                          obj.grad=1;
                          obj.OutDataTypeStr='double';
                          obj.rel_obj={};
                          obj.rel_mat=[];
                          obj.img = 0; % obj.img = 0;
                          obj.r= 0; %obj.r= obj; 
          end
               


          end
             
      
  %% 
     
     function set_cov ( x1 , x2 , cov )
         
         % set_cov ( x1 , x2 , cov ): Set the covariance between two uncertainty objects x1 and x2 to the value cov
           
        if (x1.img==0) &&  (x2.img==0)  
            
             if ~(size(x1,1)==size(x1,2)==size(x2,1)==size(x2,2)==1) 
               error(' Input argument must be a 1x1 array!');
             end
           
             if ( cov > (x1.std_unc * x2.std_unc) ) || ( cov < -(x1.std_unc * x2.std_unc) )
                 
                 error('Incorrect covariance value!\nCovariance must have a value between %f and %f .'...
                         , -x1.std_unc * x2.std_unc , x1.std_unc * x2.std_unc );
                     
             end
             
             if (x1 == x2)
                 error('Covariance can only be set for two different unc objects. The variance of this object is %f.',x1.std_unc^2);
             end
             
             
             % If x1 is already correlated with some other unc objects, then search if x2 is one of them
             ind1 = zeros(1,length(x1.rel_obj));
             for i =1:length(x1.rel_obj)
              ind1(i) = (x2 == x1.rel_obj{i});
             end 
             indx = find(ind1);
             %check if exists objects in comun in the dependency list
             [c1,ia1,ib1] = intersect(x1.dep,x2.dep);
             
             if isempty(indx) && all(x1.dep == x1) && all(x2.dep==x2) % case where x1 and x2 have not yet been assigned any covariance value
                 
                 % Set the covariance property for x1
                    x1.rel_obj{length(x1.rel_obj) + 1} = x2 ;  
                    x1.rel_mat(1, length(x1.rel_obj)) = cov ; 
               
                 % Update the correlation property for x1 
                 x1.rel_mat(2, length(x1.rel_obj)) = cov/(x1.std_unc * x2.std_unc) ;
                 
                 % Set the covariance property for x2                 
                     x2.rel_obj{ length(x2.rel_obj) + 1 } = x1 ;  
                     x2.rel_mat(1, length(x2.rel_obj) ) = cov ;
                 
                 % Update the correlation property for x2 
                 x2.rel_mat(2, length(x2.rel_obj) ) = cov/(x1.std_unc * x2.std_unc) ;
                 
                 % update the dep property.
                 CX = get_cov_mat([x1 x2]);
                 [~,~]=karhloev([x1 x2],CX,'cov');
                                  
              % case x1 and x2 have been correlated with others unc objects, but they dont have any comun dep    
%              elseif ((length(x1.dep)>=1) || (length(x2.dep)>=1)) && isempty(c1)
%                  % Set the covariance property for x1
%                     x1.rel_obj{length(x1.rel_obj) + 1} = x2 ;  
%                     x1.rel_mat(1, length(x1.rel_obj)) = cov ; 
%                
%                  % Update the correlation property for x1 
%                  x1.rel_mat(2, length(x1.rel_obj)) = cov/(x1.std_unc * x2.std_unc) ;
%                  
%                  % Set the covariance property for x2                 
%                      x2.rel_obj{ length(x2.rel_obj) + 1 } = x1 ;  
%                      x2.rel_mat(1, length(x2.rel_obj) ) = cov ;
%                  
%                  % Update the correlation property for x2 
%                  x2.rel_mat(2, length(x2.rel_obj) ) = cov/(x1.std_unc * x2.std_unc) ;
%                  
%                  
%                  %% update the dep property.
%                  % save the dependency list of x1 and x2 (when dep list is
%                  % not the obj itself.)
%                  if all(x1.dep ~= x1) && all(x2.dep == x2)
%                      x1_dep = x1.dep;
%                      x1_grad = x1.grad;
%                      x1_flag = 1; % this flag indicate that the dep of x1 have to be restored
%                      x2_flag = 0; % this flag indicate that the dep of x2 doesnt have to be restored
%                  elseif all(x2.dep ~= x2) && all(x1.dep == x1)
%                      x2_dep = x2.dep;
%                      x2_grad = x2.grad;
%                      x2_flag = 1;
%                      x1_flag = 0;
%                  elseif all(x2.dep ~= x2) && all(x1.dep ~= x1)
%                      x1_dep = x1.dep;
%                      x1_grad = x1.grad;
%                      x2_dep = x2.dep;
%                      x2_grad = x2.grad;
%                      x2_flag = 1;
%                      x1_flag = 1; 
%                  end
%                  %  
%                  CX = get_cov_mat([x1 x2]);
%                  [~,~]=karhloev([x1 x2],CX,'cov');
%                  
%                  %restore dep and grad property
%                  if x1_flag ==1
%                  x1.dep = [x1.dep x1_dep];
%                  x1.grad = [x1.grad x1_grad];
%                  elseif x2_flag ==1
%                  x2.dep = [x2.dep x2_dep];
%                  x2.grad = [x2.grad x2_grad];
%                  end
             elseif ((length(x1.dep)>=1) || (length(x2.dep)>=1)) && isempty(c1)
                 warning ('The input variables are already correlated with others uncertainty objects.');
                 
             else % case where x1 and x2 have already been assigned a covariance value
                                 
                warning('The input variables are already correlated!');
                 
                 
             end

        else
            error('For complex unc objects, covariance must be set for real and imaginary parts.')
                
        end
        
     end
     
     
 %% 
     
     function set_correl ( x1 , x2 , cor )
         
         % set_correl ( x1 , x2 , cor ): Set the correlation coefficient between two uncertainty objects x1 and x2 to the value cor
         
        if (x1.img==0) &&  (x2.img==0) 
            
             if ~(size(x1,1)==size(x1,2)==size(x2,1)==size(x2,2)==1) 
               error(' Input argument must be a 1x1 array!');
             end             
            
             if (cor > 1) || (cor < -1)
                 error('Correlation value must be between -1 and +1!');
             end
             
             if (x1 == x2)
                 error('Correlation can only be set for two different unc objects.');
             end
             
             % If x1 is already correlated with some other variables, then search if x2 is one of them:
             ind1 = zeros(1,length(x1.rel_obj));
             for i =1:length(x1.rel_obj)
              ind1(i) = (x2 == x1.rel_obj{i});
             end 
             indx = find(ind1); 
             %check if exists objects in comun in the dependency list
             [c1,~,~] = intersect(x1.dep,x2.dep);
             
             if isempty(indx) && all(x1.dep == x1) && all(x2.dep==x2) % case where x1 and x2 have not yet been assigned any correlation.
                 
                 % Set the correlation property for x1 variable                 
                     x1.rel_obj{ length(x1.rel_obj) + 1 } = x2 ;  
                     x1.rel_mat(2, length(x1.rel_obj)) = cor ;                
                 
                 % Update the covariance property for x1 variable
                 x1.rel_mat(1, length(x1.rel_obj)) = cor * x1.std_unc * x2.std_unc ;
                 
                 % Set the correlation property for x2 variable
                     x2.rel_obj{ length(x2.rel_obj) + 1 } = x1 ;  
                     x2.rel_mat(2, length(x2.rel_obj) ) = cor ;
                                  
                 % Update the covariance property for x2 variable
                 x2.rel_mat(1, length(x2.rel_obj) ) = cor * x1.std_unc * x2.std_unc ;
                 
                 % update the dep property.
                 CX = get_cov_mat([x1 x2]);
                 [~,~]=karhloev([x1 x2],CX,'cov');
                 

                 
             elseif ((length(x1.dep)>=1) || (length(x2.dep)>=1)) && isempty(c1)
                 warning ('The input variables are already correlated with others uncertainty objects.');
                 
             else % case where x1 and x2 have already been assigned a covariance value
                                 
                warning('The input variables are already correlated!');
             end
                       
        else
            error('For complex unc objects, correlation must be set for real and imaginary parts.')
                
        end
   
     
     end
     
     
   
  
   %% 
     
       function [CX_cov]=cor_into_cov_mat(X,CX_cor)
          
          % This method transforms a correlation matrix into a covariance matrix 
          %
          % X: array of unc objects
          % CX_cor: correlation matrix of X
          % CX_cov: covariance matrix of X
          %
          if ~isreal(gmv(X)) && ~isreal(CX_cor)
              error('The elements of the covariance matrix must be real.');
          else
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
                       CX_cov(i,j)=CX_cor(i,j) * X(i).std_unc * X(j).std_unc ;
                       
                    end 
                                                      
              end
              
          else
              error('Correlation matrix must be symmetric and positive semidefinite');
          end 
          end 
          
       end
           

    %% 
     
       function [CX_cor]=cov_into_cor_mat(X,CX_cov)
          
          % This method transforms a covariance matrix into a correlation matrix 
          %
          % X: array of unc objects
          % CX_cov: covariance matrix of X
          % CX_cor: correlation matrix of X
          % Cx_real: covariance matrix of the real part of X
          % Cx_img: covariance matrix of the imaginary part of X 
          % Cx_ir: covariance matrix between the imaginary and real part of X
          % Cx_ri: covariance matrix between the real and imaginary part of X
          [~,p] = chol(CX_cov);
          n=size(CX_cov,2); 
          CX_cor=zeros(n,n);
          
          m = length(X);
% %           Cx_corComplex = zeros(m);
                   
           % check if CX_cov and X have compatible dimension
           if (n~= length(X)) && (n~= 2*length(X))
               error('Dimension between the input parameters must agree.');
           
           %%check if the correlation matrix CX_cor is symmetric and PSD 
           elseif isequal(CX_cov,transpose(CX_cov)) && (isreal(gmv(X))) && (p==0)  
                       
              for i=1:n
                    for j=1:n
                       %check if the values of covariances
                       %are inside the range. 
                       if (i~=j) && (( CX_cov(i,j) > (X(i).std_unc * X(j).std_unc) ) || ( CX_cov(i,j) < -(X(i).std_unc * X(j).std_unc) ))
                       error('Incorrect value in entry(%d,%d) of the covariance matrix! \n The covariance between two random variables x and y must satisfy the property |cov(x,y)|%s %sx*%sy, where %sx and %sy represent the standard deviation of x and y respectively.'...
                       ,i,j,num2str(char(8804)),num2str(char(963)),num2str(char(963)),num2str(char(963)),num2str(char(963)));
                       end
                       CX_cor(i,j)= CX_cov(i,j) / (X(i).std_unc * X(j).std_unc);

                    end
              end
           elseif ~isreal(gmv(X)) && isreal(CX_cov) && isequal(CX_cov,transpose(CX_cov)) 
               X1 = real(X);
               X2 = imag(X);
               Xext = [X1; X2]; 
               Xext = Xext(:)'; 
            
               %check if the values of covariances are inside the range. 
                 for i=1:n
                    for j=1:n
                       if (i~=j) && (( CX_cov(i,j) > (Xext(i).std_unc * Xext(j).std_unc) ) || ( CX_cov(i,j) < -(Xext(i).std_unc * Xext(j).std_unc) ))
                       error('Incorrect value in entry(%d,%d) of the covariance matrix! \n The covariance between two random variables x and y must satisfy the property |cov(x,y)|%s %sx*%sy, where %sx and %sy represent the standard deviation of x and y respectively.'...
                       ,i,j,num2str(char(8804)),num2str(char(963)),num2str(char(963)),num2str(char(963)),num2str(char(963)));
                       end
                       CX_cor(i,j)= CX_cov(i,j) / (Xext(i).std_unc * Xext(j).std_unc);
                    end
                 end
                 
% %               if ~isreal(CX_complex) && nargin>=2
% %               % complex correlation coefficient is used as a measure for the degree of impropriety of X.
% %               % properness iimplies that cov(x_real,x_img)=0 and
% %               % var(x_real)= var(x_img)
% %               for i=1:m
% %                   for j=1:m
% %                             corr_real(i,j) = (Cx_real(i,j)^2 - Cx_img(i,j)^2 - Cx_ir(i,j)^2 + Cx_ri(i,j)^2)/...
% %                                 ((Cx_real(i,j) + Cx_img(i,j))^2 + (Cx_ir(i,j) - Cx_ri(i,j))^2);
% %                             
% %                             corr_img(i,j) = (2*(Cx_real(i,j)*Cx_ri(i,j) + Cx_img(i,j)*Cx_ir(i,j)))/...
% %                                 ((Cx_real(i,j) + Cx_img(i,j))^2 + (Cx_ir(i,j) - Cx_ri(i,j))^2);
% %                             
% %                             Cx_corComplex(i,j)= complex(corr_real(i,j),corr_img(i,j));
% %                             
% %                            % just for check if the matrix are equivalent 
% %                             corrNum(i,j) = complex( Cx_real(i,j) - Cx_img(i,j), Cx_ir(i,j) + Cx_ri(i,j));
% %                             corrDen(i,j) = complex( Cx_real(i,j) + Cx_img(i,j), Cx_ir(i,j) - Cx_ri(i,j));
% %                             CX_cor1(i,j) = corrNum(i,j) / corrDen(i,j);
% %                       
% %                   end 
% %               end
% %               
% %               end 
           else
              error('Correlation matrix must be symmetric and positive semidefinite'); 
           end
           
  
         end  
   
        
  %%

     function [K, Coef_mat]=karhloev( X,CX,arg3 )
            
            % This method generates uncorrelated unc objects K1,K2,...,Kn 
            % from correlated unc objects x1,x2,...,xn using Karhunen-Loeve decomposition.
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
 
                else
                    if nargin==3
                        Z=X;
                    end
                end
                
                

                if (strcmp(arg3,'cov') || strcmp(arg3,'corr')) && isa(CX,'double') ...
                        && ( ( size(CX,1)==length(Z) && size(CX,2)==length(Z) ) || (arg==1) ) && size(Z,1)==1 
                       
                    
                    if strcmp(arg3,'corr')
                        CX=cor_into_cov_mat(Z,CX);
                    end
                    
                    %%% check if CX is diagonal, then the objects are
                    %%% uncorrelated
                    if isdiag(CX) 
                        K = X';
                        Coef_mat = eye(length(X));
                    else

%                     % Covariance assignment
%                     unc(Z,CX,'cov');

                    % Eigenvectors 

                    [ E , eigD] = eig(CX);
                    
                    %check if eigenvalues = 0(consider eig < eps as 0) . 
                    eig0 = find(diag(eigD) < eps);
                    %remove the eigenvectors corresponding to eigenvalues 0
                    E(:,eig0)= []; 

                    % Generating the uncorrelated unc objects K1,K2,...Kn                    
                    K = E' * Z';

                    % Coefficients matrix
                    Coef_mat= E;
                    
                    % updating dep property for uncorrelated variables
                       for k = 1:length(K)
                          K(k).dep = K(k);
                          K(k).grad = 1;
                       end
                    % updating dep property for initial variables
% % % % % % % % % % Subrouting: if Coef_mat(i,j)=0, K(j) it is not included in the dependency list of X(i).      
                       for i=1:length(X) % index the original variables                        
                           for j = 1:size(Coef_mat,2) % index the uncorrelated variables
                               if Coef_mat(i,j)~=0
                                  X(i).dep(j) = K(j);
                                  X(i).grad(j) = Coef_mat(i,j);                              
                               end 
                           end
                           s = eval_std_unc(X(i));
                           X(i).std_unc = s.std_unc;
                       end
% % % % % % % % % % % % % % % % % % 
                    end 
                else
                    error (' Incorrect input arguments!');

                end
                
        else
            
            error('Incorrect number of input arguments!');
        end

        
     end
            
            
 
%% 
            
function [cov,CY_cov] = covr( y1 , y2 )
    
    % Covariance between two unc objects y1 and y2

    if (y1.img~=0) && (y2.img~=0) && (y1~=y2)
        cov=zeros(4,4);
        A=[y1.r y1.img y2.r y2.img];

        for i=1:4
            for j=1:4

           cov(i,j)=covr(A(i),A(j));     

            end
        end

        CY_cov=cov;

    elseif (y1.img~=0) && (y2.img~=0) && (y1==y2)

            cov=zeros(2,2);
            A=[y1.r y1.img];

            for i=1:2
                for j=1:2

                    cov(i,j)=covr(A(i),A(j));     

                end
            end

            CY_cov=cov;

    elseif ( ( (y1.img~=0) && (y2.img==0) ) ||  ( (y1.img==0) && (y2.img~=0) ) ) && (y1~=y2)

                if y1.img==0
                    a=y1;
                    b=y2.r;
                    c=y2.img;
                elseif y2.img==0
                        a=y2;
                        b=y1.r;
                        c=y1.img;
                  
                end

                cov(1)=covr(a,b);
                cov(2)=covr(a,c);

                CY_cov=cov;

    elseif (y1.img==0) && (y2.img==0)

                    %-----------------------------------------------
                    cov = 0; %initialize value of covariance
                    % Case where y1 and y2 are identical
                   if (y1==y2)
                       cov=(y1.std_unc)^2; % cov(y1,y1)=variance(y1)
                       CY_cov=[];
                   
                   % Case where y1 and y2 have been just created for the first time
                   elseif length(y1.grad)<=1 && length(y2.grad)<=1 && y1.grad ==1 && y2.grad==1  
                             ind1 = zeros(length(y1.rel_obj));
                             for i =1:length(y1.rel_obj)
                                ind1(i) = (y2 == y1.rel_obj{i});
                             end 
                            indx = find(ind1);

    %                          indx = find( strcmp( y2.name ,y1.rel_obj.name ) );

                             if isempty(indx)

                                 cov = 0;
                                 CY_cov=[];

                             else

                                 cov = y1.rel_mat(1,indx);
                                 CY_cov=[];

                             end
                     
% % %                     % Case y1 and y2 have already been assigned a covariance
% % %                     
% % %                     % check if y2 is one of the related objects with y1
% % %                      
% % %                     ind1 = zeros(length(y1.rel_obj));
% % %                     for i =1:length(y1.rel_obj)
% % %                         ind1(i) = (y2 == y1.rel_obj{i});
% % %                     end
% % %                    indx = find(ind1);
% % % 
% % %                    % find if obj1 and obj2 have dependencies in comum 
% % %                     [c1,ia1,ib1] = intersect(y1.dep,y2.dep);
% % % 
% % %                    %assign the value of cov to the output  
% % %                    if ~isempty(indx)  
% % % 
% % %                        cov=y1.rel_mat(1, indx );
% % %                        CY_cov=[];               
% % % 
% % %                    elseif ~isempty(c1) && isempty(indx) 
% % % 
% % %                           for k=1:length(c1)
% % %                                cov = cov + y1.grad(ia1(k))*y2.grad(ib1(k))*gmu(y1.dep(ia1(k)))^2;
% % %                           end
% % %                           CY_cov=[];                   

                                              
                   else
                       
                    
                        % Case: y1 or y2 is output function
                        var_objects= [ y1.dep  y2.dep ]; 
%                         var_names = [ y1.dep.name  y2.dep.name ];

                        [ C , IA , ~] = unique( var_objects ); % C contains all the variable on which y1 and y2 depend. 


                        Grad1 = zeros(1, length(C) );

                        Grad2 = zeros(1, length(C) );



                        % ---------- Generating Grad1 ------------ %

                        % Create an array that contains all the variables of y1:
%                         var_y1=cell(1, length(y1.dep));
%                         for i=1:length(y1.dep)
%                             var_y1 (i) = y1.dep(1,i).name;
%                         end

                        %  Create Grad1

                        for i=1:length(y1.dep)
                            for j=1:length(C)
                                ind(j) = ( y1.dep(i) == C(j));
                                
                            end
                            indx = find(ind);
                            Grad1(indx) = y1.grad(i);
                        end


                        % ---------- Generating Grad2 ------------ %

                        % Create an array that contains all the variables of y2:
%                         var_y2=cell(1, length(y2.dep));
%                         for i=1:length(y2.dep)
%                             var_y2 (i) = y2.dep(1,i).name;
%                         end 

                        %  Create Grad2

                        for i=1:length(y2.dep)
                            for j=1:length(C)
                                ind(j) = ( y2.dep(i) == C(j));
                                
                            end
                            indx = find(ind);
                            Grad2(indx) = y2.grad(i);
                        end


                        % ------- Get the input unc objects ----------%

                        for i=1:length(C)

                            X(i)= var_objects(IA(i));

                        end

                        % ---- Reconstruction of the Input Covariance Matrix CX ---- %

                        [CX]=get_input_cov_mat([y1 y2]);            

                        % -------- Perform KLD ------%

                        [K, Coef_mat]=karhloev( X , CX ,'cov');

                        % --------- Get the Jacobian -------%

                        J= [ Grad1 ; Grad2 ] * Coef_mat ;


                        % Get the input covariance matrix of the
                        % uncorrelated unc objects

                        CK= zeros(length(K));

                        for i=1:size(CK,1)
                            for j=1:size(CK,2)
                                if i==j
                                    CK(i,j)=(K(i).std_unc)^2; 
                                end

                            end
                        end


                        % output covariance matrix

                        CY_cov= J * CK * J' ; 

                        % covariance between y1 and y2

                        cov = CY_cov(1,2); 

                   end
                                                
% % %                    end
              %--------------------------------------------------      
    
    end

end


%% 
            
function [ cor , CY_corr ] = correl( y1 , y2 )
    
    % Correlation between two unc objects (Pearson's correlation coefficient)

    if (y1.img~=0) && (y2.img~=0) && (y1~=y2)
    
        cor=zeros(4,4);
        A=[y1.r y1.img y2.r y2.img];

        for i=1:4
            for j=1:4

           cor(i,j)=correl(A(i),A(j));     

            end
        end

        CY_corr=cor;
    
    elseif (y1.img~=0) && (y2.img~=0) && (y1==y2)
        
        cor=zeros(2,2);
        A=[y1.r y1.img];
    
        for i=1:2
            for j=1:2
            
                cor(i,j)=correl(A(i),A(j));     
            
            end
        end
        
        CY_corr=cor;
        
    elseif ( ( (y1.img~=0) && (y2.img==0) ) ||  ( (y1.img==0) && (y2.img~=0) ) ) && (y1~=y2)
            
            if y1.img==0
                a=y1;
                b=y2.r;
                c=y2.img;
            elseif y2.img==0
                    a=y2;
                    b=y1.r;
                    c=y1.img;
                
            end
                    
            cor(1)=correl(a,b);
            cor(2)=correl(a,c);
            
            CY_corr=cor;
            
    elseif (y1.img==0) && (y2.img==0)
                
                %-----------------------------------------------
                % Case where y1 and y2 are identical
                if (y1==y2)
                       cor = 1; 
                       CY_corr=[];
                % Case 1: y1 and y2 have already been assigned a correlation
                elseif length(y1.grad)<=1 && length(y2.grad)<=1 && y1.grad ==1 && y2.grad==1
                      ind1 = zeros(length(y1.rel_obj));
                      for i =1:length(y1.rel_obj)
                         ind1(i) = (y2 == y1.rel_obj{i});
                      end 
                      indx = find(ind1);

                      if isempty(indx)

                          cor = 0;
                          CY_corr=[];

                      else 

                          cor=y1.rel_mat(2, indx );
                          CY_corr=[];

                      end 
                      
                else 
                    % Case: y1 or y2 is output function

                [cov,CY]=covr(y1,y2);
                             
                  m=size(CY,1);
                  n=size(CY,2);
                  CY_corr=zeros(m,n);
                      for i=1:m
                         for j=1:n
                              CY_corr(i,j)=CY(i,j)/sqrt( CY(i,i)* CY(j,j) ); 
                         end
                      end
                             
                  cor = CY_corr(1,2);
                  
                end
                 
                             
    end
end
                
          %--------------------------------------------------      
       
   %% 
   
    function [CX_cov,CX_complex,Cx_real,Cx_img,Cx_ir,Cx_ri]=get_cov_mat(X)

       % Return the covariance matrix CX_cov of the uncertainty 
       % objects contained in the array X. When X is complex CX_Cov is the
       % 2n X 2n covariance matrix, containing the covariance values between 
       % the real and imag components and CX_complex is the nXn complex matrix resulting for 
       % CX_complex = E{(X - meanX)(X - meanX)^H}
       
       if size(X,1)==1
           if isreal(gmv(X)) && (nargout<=1)
              CX_cov = zeros(length(X));
          
              for i=1:length(X)
                   for j=1:length(X)

                      obj1=unc(X(i));
                      obj2=unc(X(j));
                      
                      % check if a specific value of covariance have been assigned between obj1 and obj2
                      ind = zeros(1,length(obj1.rel_obj));
                      for ii =1:length(obj1.rel_obj)
                        ind(ii) = (obj2 == obj1.rel_obj{ii});
                      end 
                      indx = find(ind);             
                      if (~isempty(indx)) && (i~=j)
                         CX_cov(i,j) = obj1.rel_mat(1,indx); 
                      else
                          
                      % find if obj1 and obj2 have dependencies in comum 
                      [c1,ia1,ib1] = intersect(obj1.dep,obj2.dep);
                      
                      for k=1:length(c1)
                           CX_cov(i,j)= CX_cov(i,j)+obj1.grad(ia1(k))*obj2.grad(ib1(k))*gmu(obj1.dep(ia1(k)))^2;
                      end            
                      end
                   end 
               end
           elseif (~isreal(gmv(X))) && ( nargout <=6 ) % the input vector X is complex
               % separe the vector X into its real and imaginary part,
               % which must have the same dimension of X.    
               % So when the real or the imaginary part are empty, i.e. 
               % when the elements are purely imaginary or real, we set it to 0.  
%                CX_complex = zeros(length(X));   
             CX_cov = zeros(2*length(X)); 
               X_real = real(X);
               X_img = imag(X);
               Cx_real = zeros(length(X));   % covariance matrix of the real part 
               Cx_img = zeros(length(X));    % covariance matrix of the imaginary part 
               Cx_ir = zeros(length(X));     % covariance matrix of imaginary-real
               Cx_ri = zeros(length(X));    % covariance matrix of real-imaginary
              
               for i=1:length(X)
                   for j=1:length(X)

                      obj1_r=unc(X_real(i));
                      obj2_r=unc(X_real(j)); 
                      
                      obj1_i=unc(X_img(i));
                      obj2_i=unc(X_img(j));
                      
                      % check if an specific value of covariance have been
                      % assigned to obj1_r and obj2_r
% %                       ind_r = zeros(1,length(obj1_r.rel_obj));
% %                       for ii =1:length(obj1_r.rel_obj)
% %                         ind_r(ii) = (obj2_r == obj1_r.rel_obj{ii});
% %                       end 
% %                       indx_r = find(ind_r);             
% %                       if (~isempty(indx_r)) && (i~=j)
% %                          Cx_real(i,j) = obj1_r.rel_mat(1,indx_r); 
% %                       else
                      % find the dependencies between the real parts of obj1 and obj2 
                      [c1_r,ia1_r,ib1_r] = intersect(obj1_r.dep,obj2_r.dep);

                      for k1=1:length(c1_r)
                           Cx_real(i,j)= Cx_real(i,j)+obj1_r.grad(ia1_r(k1))*obj2_r.grad(ib1_r(k1))*gmu(obj1_r.dep(ia1_r(k1)))^2;
                      end
% %                       end 
                      
                      % check if an specific value of covariance have been
                      % assigned to obj1_i and obj2_i
% %                       ind_i = zeros(1,length(obj1_i.rel_obj));
% %                       for ii =1:length(obj1_i.rel_obj)
% %                         ind_i(ii) = (obj2_i == obj1_i.rel_obj{ii});
% %                       end 
% %                       indx_i = find(ind_i);             
% %                       if (~isempty(indx_i)) && (i~=j)
% %                          Cx_img(i,j) = obj1_i.rel_mat(1,indx_i); 
% %                       else
                      % find the dependencies between the imaginaries parts of obj1 and obj2 
                      [c1_i,ia1_i,ib1_i] = intersect(obj1_i.dep,obj2_i.dep);

                      for k2=1:length(c1_i)
                           Cx_img(i,j)= Cx_img(i,j)+obj1_i.grad(ia1_i(k2))*obj2_i.grad(ib1_i(k2))*gmu(obj1_i.dep(ia1_i(k2)))^2;
                      end 
% %                       end
                      
                      % check if an specific value of covariance have been
                      % assigned to obj1_i and obj2_r
% %                       ind_ir = zeros(1,length(obj1_i.rel_obj));
% %                       for ii =1:length(obj1_i.rel_obj)
% %                         ind_ir(ii) = (obj2_r == obj1_i.rel_obj{ii});
% %                       end 
% %                       indx_ir = find(ind_ir);             
% %                       if (~isempty(indx_ir)) && (i~=j)
% %                          Cx_ir(i,j) = obj1_r.rel_mat(1,indx_ir); 
% %                       else
                      
                      % find the dependencies between the imaginary part of obj1 and real part of obj2 
                      [c1_ir,ia1_ir,ib1_ir] = intersect(obj1_i.dep,obj2_r.dep);

                      for k3=1:length(c1_ir)
                           Cx_ir(i,j)= Cx_ir(i,j)+obj1_i.grad(ia1_ir(k3))*obj2_r.grad(ib1_ir(k3))*gmu(obj1_i.dep(ia1_ir(k3)))^2;
                      end 
% %                       end
                      
                      % check if an specific value of covariance have been
                      % assigned to obj1_r and obj2_i
% %                       ind_ri = zeros(1,length(obj1_r.rel_obj));
% %                       for ii =1:length(obj1_r.rel_obj)
% %                         ind_ri(ii) = (obj2_i == obj1_r.rel_obj{ii});
% %                       end 
% %                       indx_ri = find(ind_ri);             
% %                       if (~isempty(indx_ri)) && (i~=j)
% %                          Cx_ri(i,j) = obj1_r.rel_mat(1,indx_ri); 
% %                       else
                      % find the dependencies between the imaginary part of obj1 and real part of obj2 
                      [c1_ri,ia1_ri,ib1_ri] = intersect(obj1_r.dep,obj2_i.dep);

                      for k4=1:length(c1_ri)
                          Cx_ri(i,j)= Cx_ri(i,j)+obj1_r.grad(ia1_ri(k4))*obj2_i.grad(ib1_ri(k4))*gmu(obj1_r.dep(ia1_ri(k4)))^2;
                      end 
% %                       end
                      
                      CX_cov(2*i-1,2*j-1) = Cx_real(i,j);
                      CX_cov(2*i-1,2*j) = Cx_ri(i,j);
                      CX_cov(2*i,2*j-1) = Cx_ir(i,j);
                      CX_cov(2*i,2*j) = Cx_img(i,j);
                                         
                      
                   end
               end 
               
                         
             if ( nargout>=2 ) 
               CX_complex = complex((Cx_real + Cx_img), (Cx_ir - Cx_ri));
             end
            
           end 
       else 
            error('Input argument must be a one-row array of uncertainty objects!')
       end
  
   end
   
   %% 
   
   function [CX_cor]=get_cor_mat(Y)
       
       % get_cor_mat(Y): Return the correlation matrix of the uncertainty objects contained in the array Y
       if isreal(gmv(Y))
        [CX_cov]=get_cov_mat(Y);
        CX_cor=cov_into_cor_mat(Y,CX_cov);
       elseif ~isreal(gmv(Y))
        [CX_cov,CX_complex,Cx_real,Cx_img,Cx_ir,Cx_ri] = get_cov_mat(Y);    
        [CX_cor]=cov_into_cor_mat(Y,CX_cov);
       end 
   end
       
   
    %%    
    
          function [J,Names]= jacobian(obj)
              
              % This method returns:
              % J: Jacobian matrix
              % Names: names of the variables with respect to which the
              % partial derivatives are taken
                   
              J=zeros();
              

              % If obj is just a single unc object then there is no use to
              % perform ordering operation
              
              if (size(obj,1)==size(obj,2)) &&  (size(obj,2)==1)
                  for i=1:length(obj.grad)
                        J(1,i)= obj.grad(i);
                        
                        if isempty(obj.dep(1,i).name)
                            Names(1,i)={'not named'};
                        else
                            Names(1,i)=obj.dep(1,i).name;
                        end
                  end
                  
              else
                      % If obj is an array of unc object then 
                      % performing ordering operation is necessary
                      % and also the unc objects must have names assigned
                      % to them.
                      if size(obj,2)>1 && size(obj,1)==1
                          
                          
                          
                          % ----------------------------
                          % get all variable names on which depend the unc objects
                              c=0;
                              for i=1:length(obj)
                                  for j=1:length(obj(i).dep)

                                        if isempty(obj(1,i).dep(1,j).name)
                                            c=c+1;
                                            var_names(c)={'not named'};
                                        else
                                            c=c+1;
                                            var_names(c)=obj(1,i).dep(1,j).name;
                                        end
                                  
                                  end
                              end
                              
                              
                              [ C , ~ , ~] = unique( var_names ); % C contains all the variable names but not-rpeated
                          
                          % ----------------------------
                      
                              for i=1:length(obj)
                                  for j=1:length(obj(i).grad)
                                      
                                       % ------------------
                                       % get the variables names of the unc object obj(1,i)  
                                       
                                            for n=1:length(obj(i).grad)
                                                 obj_Names(i,n)=obj(1,i).dep(1,n).name;
                                            end
                                        
                                       % now  obj_Names(i,:) contains the variables names of the unc object obj(1,i)       
                                       % ------------------


                                        % -------------------

                                         %  Create: 
                                         % [1] ordered Grad(i,j)which
                                         %     is the gradient of the unc
                                         %     object obj(1,i)
                                         % [2] ordered variables names of 
                                         %     the unc object obj(1,i) with
                                         %     respect to which the partial
                                         %     derivatives are taken
                                              

                                        for ii=1:size(obj_Names(i,:),2)
                                            for jj=1:length(C)
                                                indx = find( strcmp(obj_Names(i,ii),C));
                                                Grad(i,indx) = obj(1,i).grad(ii); 
                                                Names(i,indx) = obj(1,i).dep(1,ii).name;
                                            end
                                        end

                                 end
                              end
                              
                              % The Jacobian matrix is created as follows:
                              %   J=[ Grad(1,:); Grad(2,:);... ;Grad(n,:)]
                              J=Grad;
                              
                      else
                          error(' Number of rows of the input array must be equal to one.');
                      end
              
              end
          end
   
            
 %%
           
       function obj = plusscal(obj1,obj2) 
           
       % Method used by operator overloading methods
       
          if (imag(obj1)~=0 && imag(obj2)~=0) || (imag(obj1)==0 && imag(obj2)~=0) || (imag(obj1)~=0 && imag(obj2)==0)
              
            obj1_1= unc(real(obj1));
            obj2_1= unc(real(obj2));
            
            obj1_2=unc(imag(obj1));
            obj2_2=unc(imag(obj2));

            % Updating value property
            obj = complex( unc(obj1_1.value+obj2_1.value,0) , unc(obj1_2.value+obj2_2.value,0) ); 
            obj_X = real(obj);
            obj_Y = imag(obj);
            
            % Updating dep property
            [C_X,IA_X,IC_X]=unique([obj1_1.dep,obj2_1.dep]);
            obj_X.dep=C_X;
            
            [C_Y,IA_Y,IC_Y]=unique([obj1_2.dep,obj2_2.dep]);
            obj_Y.dep=C_Y;

            % End of dep property update %

             % Updating grad property
             % Real part
             
             obj_X.grad=zeros(1,length(IA_X));  
             Grad=[obj1_1.grad,obj2_1.grad];

             for i=1:length(IC_X)
                  obj_X.grad(IC_X(i))=obj_X.grad(IC_X(i))+Grad(i);
             end
             
             % Imaginary part
         
             obj_Y.grad=zeros(1,length(IA_Y));  
             Grad=[obj1_2.grad,obj2_2.grad];

             for i=1:length(IC_Y)
                  obj_Y.grad(IC_Y(i))=obj_Y.grad(IC_Y(i))+Grad(i);
             end
             
             


            % Updating std_unc property
            % Real part
             s_X = eval_std_unc(obj_X);
             obj_X.std_unc = s_X.std_unc;
             
             % Imaginary part
             s_Y = eval_std_unc(obj_Y);
             obj_Y.std_unc = s_Y.std_unc;
              
             % Final result
             obj = complex(obj_X , obj_Y);
              
              
          else
              if imag(obj1)==0 && imag(obj2)==0
               
                   obj1=unc(obj1);
                   obj2=unc(obj2);  

                   % Updating value property
                   obj = unc(obj1.value+obj2.value,0); 

                   % Updating dep property
                   [C,IA,IC]=unique([obj1.dep,obj2.dep]);
                   obj.dep=C;

                   % End of dep property update %

                   % Updating grad property
                   obj.grad=zeros(1,length(IA));  
                   Grad=[obj1.grad,obj2.grad];

                   for i=1:length(IC)
                       obj.grad(IC(i))=obj.grad(IC(i))+Grad(i);
                   end


                   % Updating std_unc property
                   s = eval_std_unc(obj);
                   obj.std_unc = s.std_unc;
              end
      
              
          end
       end
           
            
  %%
  
   function obj = minusscal(obj1,obj2) 
   
   % Method used by operator overloading methods    
       
         if imag(obj1)~=0 && imag(obj2)~=0 || (imag(obj1)==0 && imag(obj2)~=0) || (imag(obj1)~=0 && imag(obj2)==0)
              
            obj1_1= unc(real(obj1));
            obj2_1= unc(real(obj2));
            
            obj1_2=unc(imag(obj1));
            obj2_2=unc(imag(obj2));

            % Updating value property
            obj = complex( unc(obj1_1.value-obj2_1.value,0) , unc(obj1_2.value-obj2_2.value,0) ); 
            obj_X = real(obj);
            obj_Y = imag(obj);
            
            % Updating dep property
            [C_X,IA_X,IC_X]=unique([obj1_1.dep,obj2_1.dep]);
            obj_X.dep=C_X;
            
            [C_Y,IA_Y,IC_Y]=unique([obj1_2.dep,obj2_2.dep]);
            obj_Y.dep=C_Y;

            % End of dep property update %

             % Updating grad property
             % Real part
             
             obj_X.grad=zeros(1,length(IA_X));  
             Grad=[obj1_1.grad,-obj2_1.grad];

             for i=1:length(IC_X)
                  obj_X.grad(IC_X(i))=obj_X.grad(IC_X(i))+Grad(i);
             end
             
             % Imaginary part
         
             obj_Y.grad=zeros(1,length(IA_Y));  
             Grad=[obj1_2.grad,-obj2_2.grad];

             for i=1:length(IC_Y)
                  obj_Y.grad(IC_Y(i))=obj_Y.grad(IC_Y(i))+Grad(i);
             end
             
             


            % Updating std_unc property
            % Real part
             s_X = eval_std_unc(obj_X);
             obj_X.std_unc = s_X.std_unc;
             
             % Imaginary part
             s_Y = eval_std_unc(obj_Y);
             obj_Y.std_unc = s_Y.std_unc;
              
             % Final result
             obj = complex(obj_X , obj_Y);
              
              
          else
              if imag(obj1)==0 && imag(obj2)==0
               
                   obj1=unc(obj1);
                   obj2=unc(obj2);  

                   % Updating value property
                   obj = unc(obj1.value-obj2.value,0); 

                   % Updating dep property
                   [C,IA,IC]=unique([obj1.dep,obj2.dep]);
                   obj.dep=C;

                   % End of dep property update %

                   % Updating grad property
                   obj.grad=zeros(1,length(IA));  
                   Grad=[obj1.grad,-obj2.grad];

                   for i=1:length(IC)
                       obj.grad(IC(i))=obj.grad(IC(i))+Grad(i);
                   end


                   % Updating std_unc property
                   s = eval_std_unc(obj);
                   obj.std_unc = s.std_unc;
              end
      
              
          end
               
              
       end         
           
           
  %%
  
  function obj = multscal(obj1,obj2)
           
       % Method used by operators overloading methods
       
       if (imag(obj1)~=0 && imag(obj2)~=0) || (imag(obj1)==0 && imag(obj2)~=0) || (imag(obj1)~=0 && imag(obj2)==0)
           
          
            obj1_1= unc(real(obj1));
            obj2_1= unc(real(obj2));
            
            obj1_2=unc(imag(obj1));
            obj2_2=unc(imag(obj2));

            % Updating value property
            obj = complex( unc( obj1_1.value * obj2_1.value - obj1_2.value * obj2_2.value , 0 ) ,...
                           unc( obj1_1.value * obj2_2.value + obj1_2.value * obj2_1.value , 0 ) ); 
                       
            obj_X = real(obj);
            obj_Y = imag(obj);
            
            % Updating dep property
            [C_X,IA_X,IC_X]=unique([obj1_1.dep,obj1_2.dep,obj2_1.dep,obj2_2.dep]);
            obj_X.dep=C_X;
            
            [C_Y,IA_Y,IC_Y]=unique([obj1_1.dep,obj1_2.dep,obj2_1.dep,obj2_2.dep]);
            obj_Y.dep=C_Y;

            % End of dep property update %

             % Updating grad property
             % Real part
             
             obj_X.grad=zeros(1,length(IA_X));  
             Grad=[ obj1_1.grad * obj2_1.value  ,-obj1_2.grad * obj2_2.value...
                     , obj2_1.grad * obj1_1.value , -obj2_2.grad * obj1_2.value ];

             for i=1:length(IC_X)
                  obj_X.grad(IC_X(i))=obj_X.grad(IC_X(i))+Grad(i);
             end
             
             % Imaginary part
         
             obj_Y.grad=zeros(1,length(IA_Y));  
             Grad=[ obj1_1.grad * obj2_2.value  , obj1_2.grad * obj2_1.value...
                     , obj2_1.grad * obj1_2.value , obj2_2.grad * obj1_1.value ];

             for i=1:length(IC_Y)
                  obj_Y.grad(IC_Y(i))=obj_Y.grad(IC_Y(i))+Grad(i);
             end
             
             


            % Updating std_unc property
            % Real part
             s_X = eval_std_unc(obj_X);
             obj_X.std_unc = s_X.std_unc;
             
             % Imaginary part
             s_Y = eval_std_unc(obj_Y);
             obj_Y.std_unc = s_Y.std_unc;
              
             % Final result
             obj = complex(obj_X , obj_Y);
              
               
       else
           if imag(obj1)==0 && imag(obj2)==0    
                                
              obj1=unc(obj1);
              obj2=unc(obj2); 
              
              obj = unc(obj1.value*obj2.value,0);
                
              [C,~,IC]=unique([obj1.dep,obj2.dep]);
                                      
              obj.dep=C;
              obj.grad = zeros(1,length(C));

              Grad=[obj1.grad*obj2.value,obj2.grad*obj1.value];
                   
              for i=1:length(IC)
                   obj.grad(IC(i))=obj.grad(IC(i))+Grad(i);
              end
              s = eval_std_unc(obj);
              obj.std_unc = s.std_unc;
              
           end
       end
             
              

  end
  
  
  %%     
 
       function obj = divscal(obj1,obj2)
           
           % Method used by operator overloading methods
           
           if (imag(obj1)~=0 && imag(obj2)~=0) || (imag(obj1)==0 && imag(obj2)~=0) || (imag(obj1)~=0 && imag(obj2)==0)
               
               mrdivide(obj1,obj2)
               
           else
               if imag(obj1)==0 && imag(obj2)==0
               
           
               obj1=unc(obj1);
               obj2=unc(obj2);
               
               obj = unc(obj1.value/obj2.value,0);
                            
               [C,IA,IC]=unique([obj1.dep,obj2.dep]);
               obj.dep=C;
               obj.grad=zeros(1,length(IA));
               
               Grad=[obj1.grad/obj2.value,-obj2.grad*obj1.value/(obj2.value)^2];
                          
               for i=1:length(IC)
                    obj.grad(IC(i))=obj.grad(IC(i))+Grad(IC(i));
               end
               s = eval_std_unc(obj);
               obj.std_unc = s.std_unc;     
               
               end
           end
           
       end
           
           
         %%
         
       function obj = mtimes(obj1,obj2)
           
           % ------------------------------------------------------------ %
           %  mtimes(obj1,obj2): Overloading of matrix multiplication operator
           %  for uncertainty objects obj1 and obj2
           % ------------------------------------------------------------ %
            
           % Handling the case of the product of a complex number and
           % a complex unc object
             
              if isa(obj1,'unc')==0 &  (~isreal(obj1))
                  
                    obj1=complex( unc(real(obj1)) , unc(imag(obj1)) );
               
              end
              
              if isa(obj2,'unc')==0 & (~isreal(obj2))
                  
                    obj2=complex( unc(real(obj2)) , unc(imag(obj2)) );
               
              end
          
              % End of handling routine

              obj1=unc(obj1);
              obj2=unc(obj2);
              
              if (length(obj1(:))==1)||(length(obj2(:))==1)
                  obj=obj1.*obj2;
              
              elseif size(obj1,2)==size(obj2,1)  % X*Y: the number of columns of X must
                                                 %      equal the number of rows of Y
                  obj=unc(zeros(size(obj1,1),size(obj2,2)));
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
         
          function obj = plus(obj1,obj2)
              
              % Overloading of the Addition Operator
              
              
              %check the dim
              %%% one of the objects is scalar 
              if length(obj1)==1
                  obj(size(obj2,1),size(obj2,2))=unc;
                  for i = 1:size(obj2,1)
                      for j = 1:size(obj2,2)
                          obj(i,j)= plusscal(obj1,obj2(i,j)); 
                      end
                  end
              elseif length(obj2)==1
                  obj(size(obj1,1),size(obj1,2))=unc;
                  
                  for i = 1:size(obj1,1)
                      for j = 1:size(obj1,2)
                          obj(i,j)= plusscal(obj1(i,j),obj2); 
                      end
                  end
                   %%% vector-vector/ matrix-matrix addition. Objects have the same dim   
              elseif size(obj1)==size(obj2)
                  obj(size(obj2,1),size(obj2,2))=unc;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= plusscal(obj1(i,j),obj2(i,j));                        
                      end                    
                  end
                  %%% vector-matrix addition 
              elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                  obj(size(obj2,1),size(obj2,2))=unc;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= plusscal(obj1(i),obj2(i,j));                        
                      end                    
                  end
              elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                  obj(size(obj2,1),size(obj2,2))=unc;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= plusscal(obj1(j),obj2(i,j));                        
                      end                    
                  end
              elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                  obj(size(obj1,1),size(obj1,2))=unc;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= plusscal(obj2(i),obj1(i,j));                        
                      end                    
                  end
              elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                  obj(size(obj1,1),size(obj1,2))=unc;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= plusscal(obj2(j),obj1(i,j));                        
                      end                    
                  end
                  %%% vector-vector addition
               elseif (isrow(obj1) && iscolumn(obj2))
                  obj(length(obj2),length(obj1))=unc;
                  for i=1:length(obj2)
                      for j=1:length(obj1)
                          obj(i,j)= plusscal(obj2(i),obj1(j));                        
                      end                    
                  end
               elseif (isrow(obj2) && iscolumn(obj1))
                  obj(length(obj1),length(obj2))=unc;
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
                  obj(size(obj2,1),size(obj2,2))=unc;
                  for i = 1:size(obj2,1)
                      for j = 1:size(obj2,2)
                          obj(i,j)= minusscal(obj1,obj2(i,j)); 
                      end
                  end
              elseif length(obj2)==1
                  obj(size(obj1,1),size(obj1,2))=unc;
                  
                  for i = 1:size(obj1,1)
                      for j = 1:size(obj1,2)
                          obj(i,j)= minusscal(obj1(i,j),obj2); 
                      end
                  end
                      
              elseif size(obj1)==size(obj2)
                  obj(size(obj2,1),size(obj2,2))=unc;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= minusscal(obj1(i,j),obj2(i,j));                        
                      end                    
                  end
                %%% vector-matrix substraction 
              elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                  obj(size(obj2,1),size(obj2,2))=unc;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= minusscal(obj1(i),obj2(i,j));                        
                      end                    
                  end
              elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                  obj(size(obj2,1),size(obj2,2))=unc;
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= minusscal(obj1(j),obj2(i,j));                        
                      end                    
                  end
              elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                  obj(size(obj1,1),size(obj1,2))=unc;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= minusscal(obj1(i,j),obj2(i));                        
                      end                    
                  end
              elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                  obj(size(obj1,1),size(obj1,2))=unc;
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= minusscal(obj1(i,j),obj2(j));                        
                      end                    
                  end
                  %%% vector-vector substraction
               elseif (isrow(obj1) && iscolumn(obj2))
                  obj(length(obj2),length(obj1))=unc;
                  for i=1:length(obj2)
                      for j=1:length(obj1)
                          obj(i,j)= minusscal(obj1(j),obj2(i));                        
                      end                    
                  end
               elseif (isrow(obj2) && iscolumn(obj1))
                  obj(length(obj1),length(obj2))=unc;
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
         
          function obj = prod(obj1)
              
               % obj=prod(obj1): Product of the elements of the object obj1. If obj1 is a
               % matrix of uncertainty objects, obj is a row vector with the product over each column. 
             
             if size(obj1,1)==1 %obj1 is a vector
                 obj = unc(1);
                 for i= 1:length(obj1)
                      obj = obj * obj1(i);
                 end
                
             else
                 obj(1,size(obj1,2))=unc;
                 for j=1: size(obj1,2) % obj1 is a matrix
                     obj3=unc(1);                                       
                     for i=1:size(obj1,1)
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
                   
                    obj=unc(eye(r1));

                    % matrix inversion
                    for j = 1 : r1 % r.value the nbr of rows of matrix
                        for i = j : r1
                            if gmv(obj1(i,j))~= 0
                                for k = 1 : r1
                                     % rows permutation
                                    s = obj1(j,k); obj1(j,k) = obj1(i,k); obj1(i,k) = s;
                                    s = obj(j,k); obj(j,k) = obj(i,k); obj(i,k) = s;
                                end
                                
                                if imag(obj1(j,j))~=0
                                            t=1/obj1(j,j);
                                else
                                            t=divscal(1,obj1(j,j));
                                            
                                end
                                
                                
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
                  
           % If input matrix is not square, stop function
           if r1 ~= c1
                  error('Matrix must be Square.')                     
           end  
            
           
           for i=1:r1
               for j=1:c1  
                 if i==j 
                    obj(i) = obj1(i,j);
                 end     
               end
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
          
          function obj=rdivide(obj1,obj2)
              
%               This method implements the overloading of the right element-wise 
%               division operator (./) for two arrays obj1 and obj2 of uncertainty objects.

               if size(obj1)==size(obj2) 
                   obj=unc(zeros(size(obj1)));
                   for i=1:size(obj1,1)
                       for j=1:size(obj1,2)
                           obj(i,j)=obj1(i,j)/obj2(i,j);
                       end
                   end
               elseif length(obj2(:))==1
                   obj=unc(zeros(size(obj1)));
                   for i=1:size(obj1,1)
                       for j=1:size(obj1,2)
                           obj(i,j)=obj1(i,j)/obj2;
                       end
                   end
               elseif length(obj1(:))==1
                   obj=unc(zeros(size(obj2)));
                   for i=1:size(obj2,1)
                       for j=1:size(obj2,2)
                           obj(i,j)=obj1/obj2(i,j);
                       end
                   end
             %%% vector-matrix element wise division 
              elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                  obj=unc(zeros(size(obj2)));
                   for i=1:size(obj2,1)
                       for j=1:size(obj2,2)
                           obj(i,j)=obj1(i)/obj2(i,j);
                       end
                   end
                  
              elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                  obj=unc(zeros(size(obj2)));
                   for i=1:size(obj2,1)
                       for j=1:size(obj2,2)
                           obj(i,j)=obj1(j)/obj2(i,j);
                       end
                   end
              elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                  obj=unc(zeros(size(obj1)));
                   for i=1:size(obj1,1)
                       for j=1:size(obj1,2)
                           obj(i,j)=obj1(i,j)/obj2(i);
                       end
                   end
              elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                  obj=unc(zeros(size(obj1)));
                   for i=1:size(obj1,1)
                       for j=1:size(obj1,2)
                           obj(i,j)=obj1(i,j)/obj2(j);
                       end
                   end
                  %%% vector-vector division
              elseif (isrow(obj1) && iscolumn(obj2))
                  obj=unc(zeros(length(obj2),length(obj1)));
                   for i=1:length(obj2)
                       for j=1:length(obj1)
                           obj(i,j)=obj1(j)/obj2(i);
                       end
                   end
              elseif (isrow(obj2) && iscolumn(obj1))
                  obj=unc(zeros(length(obj1),length(obj2)));
                   for i=1:length(obj1)
                       for j=1:length(obj2)
                           obj(i,j)=obj1(i)/obj2(j);
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
          
          function obj =mldivide(obj1,obj2)
              
%               This method implements the overloading of matrix left division  
%               for two arrays obj1 and obj2 of uncertainty objects.

               %two objects are scalar
              if isscalar(obj1) && isscalar(obj2) 
                  obj = obj2/obj1;
                  
               %obj1 is scalar and obj2 is an array
              elseif isscalar(obj1) && (ismatrix(obj2)|| isvector(obj2))
                  obj = unc(zeros(size(obj2)));
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= obj2(i,j)/obj1;                    
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
            
         if (any(imag(obj1)~=0) & any(imag(obj2)~=0)) | (any(imag(obj1)==0) & any(imag(obj2)~=0)) | (any(imag(obj1)~=0) & any(imag(obj2)==0))
         
              
            obj1_1= unc(real(obj1));
            obj2_1= unc(real(obj2));
            
            obj1_2=unc(imag(obj1));
            obj2_2=unc(imag(obj2));

            % Updating value property
            a= (obj1_1.value * obj2_1.value + obj1_2.value * obj2_2.value )/( obj2_1.value^2 + obj2_2.value^2);
            b= (-obj1_1.value * obj2_2.value + obj2_1.value * obj1_2.value )/( obj2_1.value^2 + obj2_2.value^2);
            
            obj = complex( unc( a ) ,...
                           unc( b ) ); 
                       
            obj_X = real(obj);
            obj_Y = imag(obj);
            
            % Updating dep property
            [C_X,IA_X,IC_X]=unique([obj1_1.dep,obj1_2.dep,obj2_1.dep,obj2_2.dep]);
            obj_X.dep=C_X;
            
            [C_Y,IA_Y,IC_Y]=unique([obj1_1.dep,obj1_2.dep,obj2_1.dep,obj2_2.dep]);
            obj_Y.dep=C_Y;

            % End of dep property update %

             % Updating grad property
             % Real part
             a = obj2_1.grad * obj1_1.value * (obj2_1.value^2 + obj2_2.value^2)-...
                 (obj1_1.value * obj2_1.value + obj1_2.value * obj2_2.value)*2*obj2_1.value*obj2_1.grad;
             b = obj2_2.grad * obj1_2.value * (obj2_1.value^2 + obj2_2.value^2)-...
                 (obj1_1.value * obj2_1.value + obj1_2.value * obj2_2.value)*2*obj2_2.value*obj2_2.grad;
             
             obj_X.grad=zeros(1,length(IA_X));  
             Grad=[ obj1_1.grad * obj2_1.value/(obj2_1.value^2 + obj2_2.value^2)  ,...
                 obj1_2.grad * obj2_2.value/(obj2_1.value^2 + obj2_2.value^2)...
                     , a/(obj2_1.value^2 + obj2_2.value^2)^2 , b/(obj2_1.value^2 + obj2_2.value^2)^2 ];

             for i=1:length(IC_X)
                  obj_X.grad(IC_X(i))=obj_X.grad(IC_X(i))+Grad(i);
             end
             
             % Imaginary part
             a = -obj2_2.grad * obj1_1.value * (obj2_2.value^2 + obj2_1.value^2)-...
                 (-obj1_1.value * obj2_2.value + obj1_2.value * obj2_1.value)*2*obj2_2.value*obj2_2.grad;
             b = obj2_1.grad * obj1_2.value * (obj2_2.value^2 + obj2_1.value^2)-...
                 (-obj1_1.value * obj2_2.value + obj1_2.value * obj2_1.value)*2*obj2_1.value*obj2_1.grad;
             
             obj_Y.grad=zeros(1,length(IA_Y));  
             Grad=[ -obj1_1.grad * obj2_2.value /(obj2_1.value^2 + obj2_2.value^2)...
                    ,obj1_2.grad * obj2_1.value/(obj2_1.value^2 + obj2_2.value^2)...
                     , b/(obj2_1.value^2 + obj2_2.value^2)^2 , a/(obj2_1.value^2 + obj2_2.value^2)^2 ];

             for i=1:length(IC_Y)
                  obj_Y.grad(IC_Y(i))=obj_Y.grad(IC_Y(i))+Grad(i);
             end
  
            % Updating std_unc property
            % Real part
             s_X = eval_std_unc(obj_X);
             obj_X.std_unc = s_X.std_unc;
             
             % Imaginary part
             s_Y = eval_std_unc(obj_Y);
             obj_Y.std_unc = s_Y.std_unc;
              
             % Final result
             obj = complex(obj_X , obj_Y);
                  
                  
                  
                  
         else
                %two objects are scalar
               if any(imag(obj1)==0) & any(imag(obj2)==0) & isscalar(obj1) & isscalar(obj2) 
                      
                  obj = obj1*inv(obj2);
                  
               %obj2 is scalar and obj1 is an array
               elseif isscalar(obj2) && (ismatrix(obj1)|| isvector(obj1))
                  obj = unc(zeros(size(obj1)));
                  for i=1:size(obj1,1)
                      for j=1:size(obj1,2)
                          obj(i,j)= obj1(i,j)*inv(obj2);                    
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
                    
     end
              
  
          
          
          %% 
          
          function obj=times(obj1,obj2)
              
%               This method implements the overloading of the element-wise 
%               multiplication operator called for the syntax ' .*'
                           
              
              % check if the objects are complex and they are not unc obj
              
              if (~isreal(obj1)) && (~isa(obj1,'unc'))                   
                 obj1 = complex(unc(real(obj1)),unc(imag(obj1)));
                 
              elseif (~isreal(obj2)) && (~isa(obj2,'unc')) 
                 obj2 = complex(unc(real(obj2)),unc(imag(obj2)));
                 
              %check if one of then is no an unc obj  
              elseif (~isa(obj1,'unc'))|| (~isa(obj2,'unc'))  
                obj1=unc(obj1);
                obj2=unc(obj2);   
                
              end         
              
              issc=false; % issc: stands for 'is scalar'
              
              % if there is a scalar, then it should be obj1
              if length(obj2(:))==1 
                  obh=obj1;
                  obj1=obj2;
                  obj2=obh;
              end
              
              if length(obj1)==1 
                  issc=true;
              end              
              
%               
              %%%case with both elements have the same size and one of them is
              %%%an scalar 
              if min(size(obj1)==size(obj2)) || issc
                  obj = unc(zeros(size(obj1)));
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
                  
                     obj = unc(zeros(length(obj1),length(obj2)));
                      for i=1:length(obj1)
                        for j=1:length(obj2)
                           obj(i,j)= multscal(obj1(i),obj2(j));
                        end
                      end
              elseif (isvector(obj1)&& isvector(obj2)) && (iscolumn(obj2)&& isrow(obj1)) 
                     obj = unc(zeros(length(obj2),length(obj1)));
                      for i=1:length(obj2)
                        for j=1:length(obj1)
                           obj(i,j)= multscal(obj2(i),obj1(j));
                        end
                      end
                      
              %%% case matrix-vector multiplication 
              elseif (isvector(obj1)&& ismatrix(obj2)) && (iscolumn(obj1) && length(obj1)==size(obj2,1))
                  obj = unc(zeros(size(obj2))); 
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= multscal(obj1(i),obj2(i,j));
                      end 
                  end 
              elseif (isvector(obj1)&& ismatrix(obj2)) && (isrow(obj1) && length(obj1)==size(obj2,2))
                  obj = unc(zeros(size(obj2))); 
                  for i=1:size(obj2,1)
                      for j=1:size(obj2,2)
                          obj(i,j)= multscal(obj1(j),obj2(i,j));
                      end
                  end
                 
              elseif (isvector(obj2)&& ismatrix(obj1)) && (iscolumn(obj2) && length(obj2)==size(obj1,1))
                  obj = unc(zeros(size(obj1)));
                      for i=1:size(obj1,1)
                          for j=1:size(obj1,2)
                              obj(i,j)= multscal(obj2(i),obj1(i,j));
                          end
                      end 
              elseif (isvector(obj2)&& ismatrix(obj1)) && (isrow(obj2) && length(obj2)==size(obj1,2))
                  obj = unc(zeros(size(obj1)));
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
              
              if imag(obj)~=0
                  
                 obj_X = real(obj);
                 obj_Y = imag(obj);
                 
                 for i=1:size(obj,1)
                     for j=1:size(obj,2)
                     
                        value(i,j) = obj_X(i,j).value + 1i * obj_Y(i,j).value;
                        
                     end
                 end
                 
                  
              else
                  
                    value=reshape([obj.value],size(obj,1),size(obj,2));
                    
              end
              
          end
          

          %% 
          
          function std_unc=gmu(obj)
              
              % gmu(obj): Get the matrix of std_unc proprties of the
              % array obj of uncertainty objects
              
               if imag(obj)~=0
                  
                 obj_X = real(obj);
                 obj_Y = imag(obj);
                 
                 for i=1:size(obj,1)
                     for j=1:size(obj,2)
                         std_unc(i,j) = sqrt( (obj_X(i,j).std_unc)^2 + (obj_Y(i,j).std_unc)^2 );
                     end
                 end
                 
                  
               else
                        std_unc=reshape([obj.std_unc],size(obj,1),size(obj,2)); 
               end
               
          end

          
           %%
          
          function names=gmn(obj)
              
              % gmn(obj): Get the matrix of name proprties of the
              % array obj of uncertainty objects 
              
              names=cell(size(obj));
              for i=1:size(obj,1)
                    for j=1:size(obj,2)
                        if isempty(obj(i,j).name)
                            names(i,j)={'not named'}; 
                        else
                            names(i,j)= obj(i,j).name;
                        end
                    end
              end
          end
          
          
          %%

           function obj=eval_std_unc(obj1)   
               
               % eval_std_unc(obj1): Evaluate the standard uncertainty of the uncertainty
               % object obj1
               
                  obj=obj1;
                  for i=1:size(obj1,1)
                      for k=1:size(obj1,2)
                             
                        % ------------------------------------------------  
                        sum1=0; % initialization  
                        sum2=0; % initialization
                        
                        for j=1:length(obj1(i).grad)
                            for p=1:length(obj1(i).grad)
                                if j~=p
                                        CX=covr( obj1(i,k).dep(j) , obj1(i,k).dep(p) ); 
                                        sum1 = sum1 + obj1(i,k).grad(j) * obj1(i,k).grad(p) * CX ;             
                                 else
                                         sum2 = sum2 + (obj1(i,k).grad(j).*[obj1(i,k).dep(j).std_unc]).^2;     
                                end
                                
                            end
                        end
                        obj(i,k).std_unc = sqrt(sum1 + sum2);
                        
                        % ------------------------------------------------ 
                      end
                  end       
           end   


          %%
          
          function  [obj] = eval_obj(str,obj1 )
              
              % Used in the overloading of common mathematical functions
              
              obj1=unc(obj1);
              
              switch str
                  
                  case 'sin'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                    obj1_X = real(obj1(i,j));
                                    obj1_Y = imag(obj1(i,j));
                            
                                    obj(i,j) = complex( cosh(obj1_Y) * sin(obj1_X) , sinh(obj1_Y) * cos(obj1_X) );
                              end
                          end
                          
                          return
                          
                      else
                          
                            obj=unc(sin(gmv(obj1)),zeros(size(obj1)));
                            dobj = cos(gmv(obj1));
                            
                      end
                      
                      
                  case 'cos'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                  
                                        obj1_X = obj1(i,j).r;
                                        obj1_Y = obj1(i,j).img;
                            
                                        a= complex( -obj1_Y , obj1_X );
                                        b= complex(  obj1_Y , -obj1_X );
                                        c=exp( a ) + exp( b );
                            
                                        obj(i,j) = complex( c.r/2 , c.img/2 );
                              end
                          end
                          
                          return
                          
                      else
                            obj = unc(cos(gmv(obj1)),zeros(size(obj1)));
                            dobj = -sin(gmv(obj1));
                      end
                          
                         
                      
                  case 'tan'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                    obj1_X = real(obj1(i,j));
                                    obj1_Y = imag(obj1(i,j));
                            
                                    obj(i,j) = complex( sin(obj1_X) * cos(obj1_X) / ( cos(obj1_X)^2 + sinh(obj1_Y)^2) , ...
                                                   sinh(obj1_Y) * cosh(obj1_Y) / ( cos(obj1_X)^2 + sinh(obj1_Y)^2) );
                              end
                          end
                          
                          return
                          
                      else
                          
                             obj=unc(tan(gmv(obj1)),zeros(size(obj1)));
                             dobj = 1 + tan(gmv(obj1)).^2;
                             
                      end
                      
                      
                  case 'acos'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                    obj(i,j) = -1i * log( obj1(i,j) + sqrt( obj1(i,j).^2 - 1 ));
                              end
                          end
                          
                          return
                          
                      else
                      
                            obj=unc(acos(gmv(obj1)),zeros(size(obj1)));
                            dobj = (-1./sqrt(1-gmv(obj1).^2));
                            
                      end
                      
                  case 'asin'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj(i,j) = -1i * log( 1i*obj1(i,j) + sqrt( 1 - obj1(i,j).^2));
                                    
                              end
                          end
                          
                          return
                          
                      else
                      
                            obj=unc(asin(gmv(obj1)),zeros(size(obj1)));
                            dobj = (1./sqrt(1-gmv(obj1).^2));
                            
                      end


                  case 'atan'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj(i,j) = (1i/2) * log( (1-1i*obj1(i,j))/(1+1i*obj1(i,j)));
                                    
                              end
                          end
                          
                          return
                          
                      else
                      
                                obj=unc(atan(gmv(obj1)),zeros(size(obj1)));
                                dobj = (1./(1+gmv(obj1).^2));
                                
                      end
                      
                      
     
                  case 'sinh'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj1_X = real(obj1(i,j));
                                    obj1_Y = imag(obj1(i,j));
                          
                                    a= complex( -obj1_X , -obj1_Y ); % a=-obj1
                            
                                    Z = exp(obj1(i,j)) - exp(a) ;
                                    
                                    obj(i,j) = complex( real(Z)/2 , imag(Z)/2 );
                                    
                              end
                          end
                          
                          return
                          
                      else
                            
                             obj=unc(sinh(gmv(obj1)),zeros(size(obj1)));
                             dobj = cosh(gmv(obj1));
                             
                      end
                      
                                
                  case 'cosh'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj1_X = real(obj1(i,j));
                                    obj1_Y = imag(obj1(i,j));
                          
                                    a= complex( -obj1_X , -obj1_Y ); % a=-obj1
                            
                                    Z = exp(obj1(i,j)) + exp(a) ;
                                    obj(i,j) = complex( real(Z)/2 , imag(Z)/2 );
                                    
                              end
                          end
                          
                          return
                          
                      else
                             obj=unc(cosh(gmv(obj1)),zeros(size(obj1)));
                             dobj = sinh(gmv(obj1));
                             
                      end
                      
                      
                  case 'tanh'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj1_X = real(obj1(i,j));
                                    obj1_Y = imag(obj1(i,j));
                          
                                    obj(i,j) = complex( sinh(2*obj1_X) / ( cosh(2*obj1_X) + cos(2*obj1_Y)) ,...
                                                     sin(2*obj1_Y) / ( cosh(2*obj1_X) + cos(2*obj1_Y)) );
                                                 
                              end
                          end
                          
                          return
                          
                      else
                          
                            obj=unc(tanh(gmv(obj1)),zeros(size(obj1)));
                            dobj = 1-tanh(gmv(obj1)).^2;
                            
                      end
                      
                  case 'acosh'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj(i,j) = log( obj1(i,j) + sqrt( obj1(i,j)*obj1(i,j)-1));
                                    
                              end
                          end
                          
                          return
                          
                      else
                          
                            obj=unc(acosh(gmv(obj1)),zeros(size(obj1)));
                            dobj = 1./sqrt(gmv(obj1).^2 - 1);
                            
                      end
                      
                      
                  case 'asinh'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj(i,j) = log( obj1(i,j) + sqrt( obj1(i,j)*obj1(i,j)+1));
                                    
                              end
                          end
                          
                          return
                          
                      else
                            obj=unc(asinh(gmv(obj1)),zeros(size(obj1)));
                            dobj = 1./sqrt(gmv(obj1).^2 + 1);
                            
                      end
                      
                      
                  case 'atanh'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj(i,j) = (1/2) * log( (1+obj1(i,j))/(1-obj1(i,j)));
                                    
                              end
                          end
                          
                          return
                          
                      else
                          
                            obj=unc(atanh(gmv(obj1)),zeros(size(obj1)));
                            dobj = 1./(1-gmv(obj1).^2);
                            
                      end
                      
                      
                  case 'cot'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj1_X = real(obj1(i,j));
                                    obj1_Y = imag(obj1(i,j));
                            
                                    obj(i,j) = complex( sin(obj1_X) * cos(obj1_X) / ( sin(obj1_X)^2 + sinh(obj1_Y)^2) , ...
                                                -sinh(obj1_Y) * cosh(obj1_Y) / ( sin(obj1_X)^2 + sinh(obj1_Y)^2) );
                                            
                              end
                          end
                          
                          return
                          
                      else
                          
                             obj=unc(cot(gmv(obj1)),zeros(size(obj1)));
                             dobj = (-1 - cot(gmv(obj1)).^2);
                      
                      end
                      
                      
                  case 'acot'
                      
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                          
                                    obj(i,j) = (-1i/2) * log( (1i+obj1(i,j))/(-1i+obj1(i,j)));
                                    
                              end
                          end
                          
                          return
                          
                      else
                      
                            obj=unc(acot(gmv(obj1)),zeros(size(obj1)));
                            dobj = (1./(1+gmv(obj1).^2));
                            
                      end

                  case 'coth'
                      
                      if imag(obj1)~=0
                          
                           for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                      
                                    obj1_X = real(obj1(i,j));
                                    obj1_Y = imag(obj1(i,j));
                          
                                    obj(i,j) = complex( sinh(2*obj1_X) / ( cosh(2*obj1_X) - cos(2*obj1_Y)) ,...
                                                    -sin(2*obj1_Y) / ( cosh(2*obj1_X) - cos(2*obj1_Y)) );
                                                
                              end
                           end
                           
                           return
                          
                      else
                          
                            obj=unc(coth(gmv(obj1)),zeros(size(obj1)));
                            dobj = (1- coth(gmv(obj1)).^2);
                            
                      end
                      
                      
                  case 'exp'
                       
                      if imag(obj1)~=0
                          
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                           
                                    obj(i,j) = exp( obj1(i,j).r) * cos( obj1(i,j).img) + 1i * exp( obj1(i,j).r) * sin( obj1(i,j).img);
                                    
                              end
                          end
                        
                          return  
                        
                      else
                          
                        obj=unc(exp(gmv(obj1)),zeros(size(obj1)));
                        dobj = exp(gmv(obj1));
                        
                      end
                      
                          
                  case 'log'
                      
                      if imag(obj1)~=0
                          
                           for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                      
                                    obj(i,j)= log( sqrt( real(obj1(i,j))^2 + imag(obj1(i,j))^2 ))+ 1i * atan2( imag(obj1(i,j)) , real(obj1(i,j)));
                                    
                              end
                           end
                           
                           return         
          
                      else
                          
                            obj=unc(log(gmv(obj1)),zeros(size(obj1)));
                            dobj = 1./gmv(obj1);
                            
                      end
                      
                  case 'log10'
                      
                      if imag(obj1)~=0
                          
                           for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                         
                                 obj(i,j)= ( log( sqrt( real(obj1(i,j))^2 + imag(obj1(i,j))^2 )) + 1i * atan2( imag(obj1(i,j)) , real(obj1(i,j))) )/log(10);
                                 
                              end
                           end
                           
                           return         
          
                      else
                      
                            obj=unc(log10(gmv(obj1)),zeros(size(obj1)));
                            dobj = 1./gmv(obj1)./log(10);
                            
                      end

                      
              end
                      

             
              for i=1:size(obj1,1)
                  for j=1:size(obj,2)
                      
                          
                            obj(i,j).grad = dobj(i,j).*[obj1(i,j).grad];
                            obj(i,j).dep  = obj1(i,j).dep;
                            
                  
                  end
              end
              
             
                obj= eval_std_unc(obj);
                
             

      
          end

 
          %% sinus
          function obj=sin(obj1)
              obj = eval_obj('sin',obj1);
                            
          end
          
          function obj=sind(obj1)
              
              obj11 = obj1/180*pi;
              obj = eval_obj('sin',obj11);
                            
          end
                  
          %% cosine
          function obj=cos(obj1)
              obj = eval_obj('cos',obj1);
          end
          
          function obj=cosd(obj1)
              
              obj11 = obj1/180*pi;
              obj = eval_obj('cos',obj11);
                            
          end
          
          %% tan
          function obj = tan(obj1)
              obj = eval_obj('tan',obj1);

          end
          
          function obj=tand(obj1)
              
              obj11 = obj1/180*pi;
              obj = eval_obj('tan',obj11);
                            
          end
                               
          
          %% Inverse cosine: acos
          function obj=acos(obj1)
              
              obj = eval_obj('acos',obj1);
          end
          
          function obj=acosd(obj1)
     
              obj11 = eval_obj('acos',obj1);
              obj = obj11/pi*180;              
          end
          
          %% Inverse sine: asin
          function obj=asin(obj1)
             
              obj = eval_obj('asin',obj1);
                                     
          end
          
          function obj=asind(obj1)
     
              obj11 = eval_obj('asin',obj1);
              obj = obj11/pi*180;              
          end
          
          %% Inverse tangent: atan
          function obj=atan(obj1)
              obj = eval_obj('atan',obj1);
 
          end
          
          function obj=atand(obj1)
     
              obj11 = eval_obj('atan',obj1);
              obj = obj11/pi*180;              
          end
          
          
          %% Hyperbolic cosine: cosh
          function obj = cosh(obj1)
              
              obj = eval_obj('cosh',obj1);
 
          end
          
          %% Hyperbolic sine: sinh
          function obj = sinh(obj1)
              
              obj = eval_obj('sinh',obj1);
 
              
          end
          
          %% Hyperbolic tangent:tanh
          function obj = tanh(obj1)
              
              obj = eval_obj('tanh',obj1);
 
          end
          
          %% Hyperbolic cosine inverse: acosh
          function obj = acosh(obj1)
              
              obj = eval_obj('acosh',obj1);

          end
          
          %% Hyperbolic sine inverse: asinh
          function obj = asinh(obj1)
                     
              obj = eval_obj('asinh',obj1);

          end
          
          %% Hyperbolic tangent inverse: atanh
          function obj = atanh(obj1)
             
              obj = eval_obj('atanh',obj1);

          end
          
          %% Cotangent : cot
          function obj = cot(obj1)
              
              obj = eval_obj('cot',obj1);

          end
          
          %% Inverse cotangent: acot
          function obj = acot(obj1)
             
              obj = eval_obj('acot',obj1);

          end
          
          %% Hyperbolic cotangent: coth
          function obj = coth(obj1)
              
              obj = eval_obj('coth',obj1);

          end
          
          %% Natural logarithm
          function obj = log(obj1)
              
              obj = eval_obj('log',obj1);

          end
          
          %% Decimal logarithm 
          
          function obj = log10(obj1)
              
              obj = eval_obj('log10',obj1);

          end
              
          
          %% exponential
          
          function obj = exp(obj1)
              
              obj = eval_obj('exp',obj1);

          end
          
                    
          %% 
          
           function obj = atan2(obj1,obj2)
               
               % Overloading the Four-quadrant inverse tangent function
               % Declaring both inputs as unc objects
              obj1 = unc(obj1);
              obj2 = unc(obj2);
                           
              obj=unc( atan2( gmv(obj1) , gmv(obj2) ), zeros(size(obj1)) );
              
              
              
              % Partial derivatives
              dobj_2 = gmv(obj1)./ (gmv(obj1).^2 + gmv(obj2).^2); 
              dobj_1 = -gmv(obj2)./ (gmv(obj1).^2 + gmv(obj2).^2);
              
              
 
              for i=1:size(obj1,1)
                  for j=1:size(obj,2)
                                            
                            [C,~,IC]=unique([obj1(i,j).dep,obj2(i,j).dep]);
                            obj(i,j).dep=C;
                            
                            obj(i,j).grad = zeros(1,length(C));
                            
                            Grad = [ dobj_1(i,j).* obj1(i,j).grad  , dobj_2(i,j).* obj2(i,j).grad ] ;
  
                            for k=1:length(IC)
                                
                                        obj.grad(IC(k))=obj.grad(IC(k))+Grad(k);
                                        
                            end
                  
                  
                  end
              end
              
           
                
              obj = eval_std_unc(obj);
               
              

           end
           
           function obj = atan2d(obj1,obj2)
                obj11 = atan2(obj1,obj2);
                obj = obj11/pi*180;
           end
               
          
                   
          %% 
          
          function obj = max(obj1)
              
              % max(obj1): Return the uncertainty object that has the
              % maximum nominal value of all the elements of obj1
              
              obj= unc(0,0);
              idx=0;
              for i=1:length(obj1)
                 if obj1(i).value==max([obj1.value])
                     idx = idx+1;
                     obj(idx)=obj1(i); 
                 end                 
                     
              end
              obj= unique(obj);            
          end

          
          %% 
          
          function obj = min(obj1)
              
              % min(obj1): Return the uncertainty object that has the
              % minimum nominal value of all the elements of obj1
              
              obj= unc(0,0);
              idx=0;
              for i=1:length(obj1)
                 if obj1(i).value==min([obj1.value])
                     idx = idx+1;
                     obj(idx)=obj1(i);                     
                 end
              end
              obj= unique(obj);            
              
          end
          
          
          %%        
          
          function obj=cumsum(obj1)
              
              % cumsum(obj1): Return the cumulative sum of the elements of obj1  

              obj=unc(0.*gmv(obj1)); % Preallocation. obj must be of the same size as obj1.
              
              W=obj1';
              Q=W(:);
              Q=Q';
              
              J=size(obj1,2);
              
              for i=1:size(obj1,1) 
                  for j=1:size(obj1,2)
                          
                          z=unc(0,0);
                        
                          for k=1:(j+(i-1)*J)
                            z = Q(k) + z;
                          end
                          obj(i,j)=z;
                       
                         
                  end
              end  
              
              
          end 
          
                    
          %% 
          
          function obj=sum(obj1)
              
              % sum(obj1): Sum of elements of obj1
              
                  obj=unc(0,0);
                  for i=1:size(obj1,1)
                    for j=1:size(obj1,2)
                          obj = obj +obj1(i,j);
                    end
                  end
             
          end 
 
          
          %% 
          
           function obj = power(obj1,obj2)
               
               % Overloading of Element-wise powers
              
              obj1=unc(obj1);
              obj2=unc(obj2);
              
             %case obj1 is complex and obj2 are real        
                 if any(any(imag(obj1)~=0)) && all(all(imag(obj2)==0))  
                        % obtain the absolute value and angle of obj1
                        obj1_r = abs(obj1);
                        obj1_p = angle(obj1);
                        % Moivre's formula for complex numbers         
                        obj=exp(obj2 .* log(obj1_r)) .* ( cos(obj2.*obj1_p) + 1i * sin(obj2.*obj1_p));
                   
                 %case obj1 and obj2 are real 
                 elseif all(all(imag(obj1)==0)) && all(all(imag(obj2)==0))
                    
                             %check the dim
                      %%% one of the objects is scalar 
                      if length(obj1)==1
                          obj(size(obj2,1),size(obj2,2))=unc;
                          for i = 1:size(obj2,1)
                              for j = 1:size(obj2,2)
                                obj(i,j) = power_check(obj1,obj2(i,j));   
                              end
                          end
                      elseif length(obj2)==1
                          obj(size(obj1,1),size(obj1,2))=unc;

                          for i = 1:size(obj1,1)
                              for j = 1:size(obj1,2)
                                  obj(i,j) = power_check(obj1(i,j),obj2); 
                              end
                          end
                           %%% vector-vector/ matrix-matrix addition. Objects have the same dim   
                      elseif size(obj1)==size(obj2)
                          obj(size(obj2,1),size(obj2,2))=unc;
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                  obj(i,j)= power_check(obj1(i,j),obj2(i,j));                        
                              end                    
                          end
                          %%% vector-matrix addition 
                      elseif (iscolumn(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,1))
                          obj(size(obj2,1),size(obj2,2))=unc;
                          for i=1:size(obj2,1)
                              for j=1:size(obj2,2)
                                  obj(i,j)= power_check(obj1(i),obj2(i,j));                        
                              end                    
                          end
                      elseif (isrow(obj1) && ismatrix(obj2) && length(obj1)==size(obj2,2))
                          obj(size(obj2,1),size(obj2,2))=unc;
                          for i=1:size(obj2,1)
                              for j=1:size(obj2,2)
                                  obj(i,j)= power_check(obj1(j),obj2(i,j));                        
                              end                    
                          end
                      elseif (iscolumn(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,1))
                          obj(size(obj1,1),size(obj1,2))=unc;
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                  obj(i,j)= power_check(obj2(i),obj1(i,j));                        
                              end                    
                          end
                      elseif (isrow(obj2) && ismatrix(obj1) && length(obj2)==size(obj1,2))
                          obj(size(obj1,1),size(obj1,2))=unc;
                          for i=1:size(obj1,1)
                              for j=1:size(obj1,2)
                                  obj(i,j)= power_check(obj2(j),obj1(i,j));                        
                              end                    
                          end
                          %%% vector-vector addition
                       elseif (isrow(obj1) && iscolumn(obj2))
                          obj(length(obj2),length(obj1))=unc;
                          for i=1:length(obj2)
                              for j=1:length(obj1)
                                  obj(i,j)= power_check(obj2(i),obj1(j));                        
                              end                    
                          end
                       elseif (isrow(obj2) && iscolumn(obj1))
                          obj(length(obj1),length(obj2))=unc;
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
               if (obj1.value > 0 )
                   obj = exp(obj2 .* log(obj1));
               elseif (obj1.value < 0 )
                   obj =((-1).^gmv(obj2)) .* exp(obj2 .* log(abs(obj1)));
               elseif (obj1.value == 0 )  
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
                obj=unc(obj1);
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
              
              L = max(size(obj1),size(obj2));
              if ismatrix(obj1) || ismatrix(obj2)
                   for i=1:L(1)
                       if size(obj1)==size(obj2)
                           obj(i,:)= power(obj1(i,:),obj2(i,:));
                           if ~isreal(obj(i).value)|| ~isreal(obj(i).std_unc)
                               error('Error using realpow. Realpow produced complex result.')
                           end
                       elseif isscalar(obj1) && ~isscalar(obj2)
                           obj(i,:)=power(obj1,obj2(i,:));
                           if ~isreal(obj(i).value)|| ~isreal(obj(i).std_unc)
                               error('Error using realpow. Realpow produced complex result.')
                           end
                       elseif ~isscalar(obj1) && isscalar(obj2)
                           
                           obj(i,:)= power(obj1(i,:),obj2);
                           if ~isreal(obj(i).value)|| ~isreal(obj(i).std_unc)
                               error('Error using realpow. Realpow produced complex result.')
                           end
                       else
                           error('Error using realpow. Matrix dimensions must agree.')
                       end
                   end
                       
              end
          end
                   
          %%
          
          function disp_contribution(obj1)
              
              % disp_contribution(obj1): Display uncertainty contributions
              % of the uncertainty object obj1
              
              contrib = [obj1.grad].*[obj1.dep.std_unc];
              [ValueArray,IndexArray] = sort(abs(contrib));
              IndexArray=IndexArray(ValueArray>max(ValueArray)*1e-16); %take only the indices whose values associated to it are greater than the maximum value times 1e-16. 
              IndexArray=IndexArray(end:-1:1); % Start with largest contributor
              
              % -------------- Creation of the Display frame ----------------------- %
              fprintf ('\nUncertainty Contribution:\n')
              
              MaxNameLength = 0;
              for Index = IndexArray
                  MaxNameLength = max(MaxNameLength,length(char(obj1.dep(1,Index).name)));
              end
              if MaxNameLength>13
                  blanks = repmat (' ',1,MaxNameLength-13);
                  StrLength = fprintf (['Variable Name',blanks,'  | Contribution\n']);
              else
                   StrLength = fprintf ('Variable Name | Contribution\n');
              end
              fprintf ([repmat('-',1,StrLength),'\n']);
              for Index = IndexArray
                  digit = floor(log10(abs(contrib(Index))));
                  disp_std_unc = abs(round(contrib(Index) / 10^(digit-2))*10^(digit-2));
                  NameStr = char(obj1.dep(1,Index).name);
                  if MaxNameLength>12
                    FillerStr = repmat ('.',1,MaxNameLength-length(NameStr));
                  else
                    FillerStr = repmat ('.',1,12-length(NameStr));  
                  end
                  fprintf ('%s %s | %s\n',NameStr,FillerStr,num2str(disp_std_unc));
              end
              % -------------- End of creation of the Display frame ----------------------- %
              
          end
         
          
          %% 
          
           function obj = gt(obj1,obj2)
           
              % This method implements the overloading of the greater 
              % than operator called for the syntax 'obj1 > obj2'
                       
           obj=unc(0,0);
           
           if isa(obj1,'unc') && isa(obj2,'unc')
               obj1=unc(obj1);
               obj2=unc(obj2);
               
               h=obj1-obj2;
               
               obj.value=double(h.value>0); % obj.value=1 if h.value>0 
                                            % obj.value=0 if h.value<0
               
               obj.std_unc=normcdf(h.value,0,h.std_unc); % obj.std_unc= area under the standard
                                                         % normal distribution curve from
                                                         % -infinity to (h.value - 0)/h.std_unc
               
           
               % obj.value represents the logical decision (in double)
               % obj.std_unc represents the probability of obj1 beeing larger than obj2
            
           else 
               % Comparison with a constant
               if (isa(obj1,'unc') && isa(obj2,'double')) || (isa(obj1,'double') && isa(obj2,'unc'))
                   
                   if isa(obj1,'unc')
                       
                       obj.value=double(obj1.value>obj2);
                       obj.std_unc=1-normcdf(obj2,obj1.value,obj1.std_unc);
                       
                   else
                       if isa(obj2,'unc')
                          
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
               
             obj=unc(0,0);
           
             if isa(obj1,'unc') && isa(obj2,'unc')
               obj1=unc(obj1);
               obj2=unc(obj2);
             
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
                if (isa(obj1,'unc') && isa(obj2,'double')) || (isa(obj1,'double') && isa(obj2,'unc'))
                   
                   if isa(obj1,'unc')
                       
                       obj.value=double(obj1.value<obj2);
                       obj.std_unc=normcdf(obj2,obj1.value,obj1.std_unc);
                       
                   else
                       if isa(obj2,'unc')
                          
                          obj.value=double(obj2.value>obj1);
                          obj.std_unc=1-normcdf(obj1,obj2.value,obj2.std_unc);
                         
                       end
                
                   end
                end
                                                       
             end  
           end
           
           %% norm 
           function obj = norm(obj1,arg)
               obj = unc;
               obj_out = unc;
               if nargin == 1 
                                 
                 for i = 1: size(obj1,1)
                     for j = 1:size(obj1,2)
                     obj_out = obj_out + abs(obj1(i,j))^2;
                     end
                 end
                 
                 obj = sqrt(obj_out);
                 
               else
                  for i = 1: size(obj1,1)
                     for j = 1:size(obj1,2)
                     obj_out = obj_out + abs(obj1(i,j))^arg;
                     end
                 end
                 obj = (obj_out)^(1/arg);
                   
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
                    obj1 = [obj1;unc(zeros((N1-N),1))];
                    N = N1; % define the new dimension of N
                end 

                xp = obj1(1:2:end);
                xpp = obj1(2:2:end);
                if N>=8
                     Xp = myfft(xp);
                     Xpp = myfft(xpp);
                     obj = unc(zeros(N,1));
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
            X = unc(X);
            Y = unc(Y);
            Xp = unc(Xp); 

            if isvector(Y) && length(Y)== N
                Yp = unc(zeros(size(Xp)));
                ind_0 = zeros(size(Xp));
                ind_1 = zeros(size(Xp));
                for i = 1:length(Xp)

                    if isempty(find(X(:)== max(X(Xp(i).value> [X(:).value])), 1)) || isempty(find(X(:)== min(X(Xp(i).value< [X(:).value])), 1))
                       Yp(i) = NaN; 
                    else
                     % find the index of the lower and upper points 
                    ind_0(i) = find(X(:)== max(X(Xp(i).value> [X(:).value])));
                    ind_1(i) = find(X(:)== min(X(Xp(i).value< [X(:).value])));   

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
           
        %% Get Jacobian           
        function Jac = get_Jac(obj1, obj2)
             % get jacobian for obj1 (vector) with respect. to obj2 (vector)
             % obj2 must comprise primary RVs

             lo1=length(obj1);
             lo2=length(obj2);
             Jac=zeros(lo1,lo2);
             for i=1:lo1
                 for j=1:lo2
                     I=find(obj1(i).dep==obj2(j));
                     if ~isempty(I)
                        Jac(i,j)=obj1(i).grad(I);
                     end
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
        H_cal = unc(zeros(numb_cal,n+1)); % information matrix 
        Y_cal = unc;
        h = unc(zeros(length(Y),n+1));
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
    
   %% overloading fft method 
   
   
          
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
                           
                     if imag(obj(i,j))~=0 
                           
                            obj_X(i,j) = real(obj(i,j));
                            obj_Y(i,j)= imag(obj(i,j));
                           
                            if obj_Y(i,j).value >= 0
                           
                               fprintf(['    \t',obj_X.GenerateDispStr(obj_X(i,j).value,obj_X(i,j).std_unc) ,' + '...
                                         obj_Y.GenerateDispStr(obj_Y(i,j).value,obj_Y(i,j).std_unc) , ' * i ']);
                           
                            else
                              fprintf(['    \t',obj_X.GenerateDispStr(obj_X(i,j).value,obj_X(i,j).std_unc) ,' - '...
                                        obj_Y.GenerateDispStr(abs(obj_Y(i,j).value),obj_Y(i,j).std_unc) , ' * i ']);
                            end
                                
                   
                       
                     else
                              if imag(obj(i,j))==0
                                  
                            
                                         fprintf(['    \t',obj.GenerateDispStr(obj(i,j).value,obj(i,j).std_unc)]);
                               
                              end  
                            
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
    
       for i=1:r1
           for j=1:c1
               
               if (x(i,j).img==0) && (y(i,j).img==0)
       
                Z(i,j)=unc(); 
                
                Z(i,j).OutDataTypeStr='complex';

                Z(i,j).r = x(i,j); % real part

                Z(i,j).img= y(i,j);  % imaginary part 
   
                Z(i,j).value= x(i,j).value + 1i * y(i,j).value; % Mean value assignment               
                Z(i,j).std_unc= sqrt( (x(i,j).std_unc)^2 + (y(i,j).std_unc)^2 ); % Standard Uncertainty assignment 
                
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
     
     for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                if Z(i,j).r ~=0
                    X(i,j)= Z(i,j).r;
                else
                    X(i,j)=0;
                end
                
         end
     end
             
 end

 
%% 
          
 function Y=imag(Z)
     
     % imag(Z): Return the imaginary part of Z
     
     % --- Handling the error occurring while converting from unc to double--- %
     % This error occurs when the array Z contains both real and complex
     % unc objects and the first element of Z is a real unc object.
     
     c=0;
     
     if (length(Z)>1) && ( Z(1,1).img == 0 ) 
         
         for i=1:size(Z,1)
           for j=1:size(Z,2)
             
               if Z(i,j).img ~=0
                   
                   k=i;
                   l=j;
                   temp=Z(1,1);
                   Z(1,1)=Z(i,j);
                   Z(i,j)=temp;
                   j=size(Z,2); % force to halt the loop
                   i=size(Z,1); % force to halt the loop
                   c=1;
                          
               end
               
           end
         end
         
         
         
     end
     
    % --- End of handling --- %
    
    
        for i=1:size(Z,1)
         for j=1:size(Z,2)
             
             if Z(i,j).img ~=0
                 
                 Y(i,j)= Z(i,j).img;
                 
             else
                 
                 Y(i,j)=0;
                 
             end
       
         end
        end
       
        
       % This will be executed in case the error handling was performed 
       if (c==1)
           
          temp=Y(1,1);
          Y(1,1)=Y(k,l);
          Y(k,l)=temp;
           
       end
        

 end
   
 
%% 
          
 function obj=conj(Z)
       
     % This method takes the complex conjugate of each
     % entry of the matrix Z
     
       for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                 obj(i,j)= Z(i,j).real - 1i * Z(i,j).img;
                
         end
       end
       
             
 end  
 
 
 
 %% 
          
 function obj=abs(Z)
   % overloading the function abs() that returns the absolute value of the uncertainty object Z. When
   % Z is complex it returns the complex modulus (magnitude) of the elements of Z.
   obj = unc(zeros(size(Z)));
       for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                 obj(i,j)= sqrt(Z(i,j).real^2 + Z(i,j).img^2);
                
         end
       end
       
             
 end  
 
 function obj=angle(Z)
      % overloading the function angle() that returns the phase angles, in radians, of the uncertainty object Z.
      obj = unc(zeros(size(Z)));
       for i=1:size(Z,1)
         for j=1:size(Z,2)
             
                 obj(i,j)= atan2(Z(i,j).img, Z(i,j).real);
                
         end
       end
       
             
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
           
          
% %            function DispString = GenerateDispStr(value,std_unc)
% %                
% %                % This method is used by the disp() method to customize 
% %                % the display of uncertainty objects
% %                
% %                 global Number_of_Digits_to_Display;
% %                 d=Number_of_Digits_to_Display-1; 
% %                 
% %                 
% %                 % 1) Clip the standard uncertainty to max 3 digits
% %                 if (std_unc==0) || (log10(abs(value)/abs(std_unc))>12)
% %                     digit=-12;
% %                 else
% %                     digit = floor(log10(abs(std_unc))); % take the integer part no matter is the value of the decimal part
% %                 end
% %                 
% %                 
% %                 
% %                 disp_std_unc = round(std_unc / 10^(digit-d)); % round the value into the closest integer.
% %                 
% %                 
% %                 % 2) Clip the value to show the last d+1 values 
% %                 disp_value = round(value / 10^(digit-d))*10^(digit-d);
% %                 
% %                 
% %                 if disp_std_unc==0
% %                         eval_str = '[num2str(disp_value),''(0)'']';
% %                         DispString = eval (eval_str);
% %                 else
% %                     if digit<0
% %                         eval_str = sprintf ('[sprintf(''%%0.0%if'',disp_value),''('',num2str(disp_std_unc),'')'']',abs(digit)+d);
% %                         DispString = eval (eval_str);
% %                     else
% %                         if digit<d
% %                             disp_std_unc = disp_std_unc*10^(digit-d); %#ok<NASGU>
% %                              eval_str = sprintf ('[sprintf(''%%0.0%if'',disp_value),''('',sprintf(''%%0.0%if'',disp_std_unc),'')'']',d-digit,d-digit);
% %                             DispString = eval (eval_str); 
% %                         else
% %                             disp_std_unc = disp_std_unc*10^(digit-d);
% %                             DispString = [num2str(disp_value),'(',num2str(disp_std_unc),')'];
% %                         end
% %                     end
% %                 end
% %            end
              
        
           
           function [obj_names,obj_array] = get_ws_obj()
               
               % Get unc class objects existing in Workspace
               % obj_names: cell array containing their names
               % obj_array: cell array containing the objects themselves
           
               var_mat=evalin('base','whos'); % the workspace of a function is separate 
                                              % from the base workspace. For
                                              % this reason evalin() is used.

               % First we count the number of class objects

                c=1; % counter
                for i=1:length(var_mat)
                   if strcmp(var_mat(i).class,'unc') 

                       c=c+1; % increment counter 
                   end
                end

                % Second, preallocation

                obj_names=cell(1,c-1);
                obj_array=cell(1,c-1);

                % Finally, store all class objects names in the cell obj_array

                c=1; % counter
                for i=1:length(var_mat)
                   if strcmp(var_mat(i).class,'unc') 
                       obj_names{1,c}=[var_mat(i).name];
                       obj_array{1,c}=evalin('base',obj_names{1,c});
                       c=c+1; % increment counter 
                   end
                end
            

           end
            
                     
                    
           function data = read_data(arg1,arg2,arg3,arg4,arg5)
              
 
              [~,~,ext] = fileparts(arg1); % get file extension
              
              % -------------- Working with Microsoft Excel files  -------------- %
              
              if strcmp(ext , '.xls') || strcmp(ext , '.xlsx')
              
                      if (nargin==5)  % read_data( file_name , sheet_name , mean_val , std_unc , names )

                                    [~,var_names] = xlsread(arg1, arg2 , arg5 );

                                    var_means = xlsread(arg1, arg2 , arg3 );

                                    var_stds = xlsread(arg1, arg2 , arg4 );

                                    data=unc(var_means , var_stds , var_names );

                      else

                          if (nargin==4) 
                              
                              % read_data( file_name , sheet_name , range , 'r/c')
                                  if strcmp(arg4,'r') % data organized in rows
                                      data=xlsread(arg1,arg2,arg3);
                                      var_means=zeros(1,size(data,1));
                                      var_stds=zeros(1,size(data,1));
                                      for i=1:size(data,1)
                                          var_means(i)=mean(data(i,:));
                                          var_stds(i)=std(data(i,:));
                                      end
                                      data=unc(var_means , var_stds );
                                  else
                                      if strcmp(arg4,'c') % data organized in columns
                                          data=xlsread(arg1,arg2,arg3);
                                          var_means=zeros(1,size(data,2));
                                          var_stds=zeros(1,size(data,2));
                                          for i=1:size(data,2)
                                              var_means(i)=mean(data(:,i));
                                              var_stds(i)=std(data(:,i));
                                          end
                                          data=unc(var_means , var_stds );
                                          
                                      else % read_data( file_name , sheet_name , mean_val , std_unc )
                              
                                            var_means = xlsread(arg1, arg2 , arg3 );

                                            var_stds = xlsread(arg1, arg2 , arg4 );

                                            data=unc(var_means , var_stds );
                                      end
                                  end
                          else
                              error('Incorrect number of input arguments!');
                          end
                          
                      end
              else
                      
              % -------------- Working with CSV files  -------------- %
              
                      if strcmp(ext , '.dat') || strcmp(ext , '.csv')    
                          
                          if (nargin==2) % read_data( file_name , 'r/c' )
                              
                              if strcmp(arg2,'r') % data organized in rows
                                      data=csvread(arg1);
                                      var_means=zeros(1,size(data,1));
                                      var_stds=zeros(1,size(data,1));
                                      for i=1:size(data,1)
                                          var_means(i)=mean(data(i,:));
                                          var_stds(i)=std(data(i,:));
                                      end
                                      data=unc(var_means , var_stds );
                                  else
                                      if strcmp(arg2,'c') % data organized in columns
                                          data=csvread(arg1);
                                          var_means=zeros(1,size(data,2));
                                          var_stds=zeros(1,size(data,2));
                                          for i=1:size(data,2)
                                              var_means(i)=mean(data(:,i));
                                              var_stds(i)=std(data(:,i));
                                          end
                                          data=unc(var_means , var_stds );
                                      else
                                          error('Incorrect input arguments!');
                                      end
                              end
                              
                           else
                              if (nargin==5)  % read_data( file_name , R , C , range , 'r/c' )
                                              % reads only the range specified
                                              % by range = [R C R2 C2] where (R,C) is the upper-left corner of
                                              % the data to be read and (R2,C2) is the lower-right corner.
                                  
                                  if strcmp(arg5,'r') % data organized in rows
                                          data=csvread(arg1,arg2,arg3,arg4);
                                          var_means=zeros(1,size(data,1));
                                          var_stds=zeros(1,size(data,1));
                                          for i=1:size(data,1)
                                              var_means(i)=mean(data(i,:));
                                              var_stds(i)=std(data(i,:));
                                          end
                                          data=unc(var_means , var_stds );
                                      else
                                          if strcmp(arg5,'c') % data organized in columns
                                              data=csvread(arg1,arg2,arg3,arg4);
                                              var_means=zeros(1,size(data,2));
                                              var_stds=zeros(1,size(data,2));
                                              for i=1:size(data,2)
                                                  var_means(i)=mean(data(:,i));
                                                  var_stds(i)=std(data(:,i));
                                              end
                                              data=unc(var_means , var_stds );
                                          else
                                              error('Incorrect input arguments!');
                                          end
                                  end
                                                         
                              else
                                     
                              error('Incorrect number of input arguments!');
                              
                              end
                          end
                           
                          
                      else
                          error('Incorrect file type!');
                      end
                  
                      
               end
                                    
                         
                              
            end
           
           
           
           
       end
end

         
          
