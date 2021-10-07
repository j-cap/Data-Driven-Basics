function obj = G_RV(mu,C)
            % Multivariate normally distributed samples
            % examples:
            % scalar:
            % x=G_RV(1,1);
            % vector (row or column):
            % x=G_RV([10,1],[1 0.1;0.1 1]): vector of jointly Gaussian RVs 
            % with covariance C =[1 0.1;0.1 1] and mean mu = [10,1]
                global MCSAMPLES
                obj=unc_t(mu);

                if isscalar(C)
                    C=diag(ones(1,length(mu))*C);
                elseif ismatrix(C) && iscolumn(C) 
                    C=C.';
                end

                R=mvnrnd(mu,C,MCSAMPLES);

                for i=1:length(mu)
                        obj(i).values=R(:,i)';
                        
                end

          end