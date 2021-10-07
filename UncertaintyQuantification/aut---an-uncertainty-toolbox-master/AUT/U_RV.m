function obj = U_RV(mu,intv)
            % uniformly distributed samples over interval intv around mu
            % examples:
            % scalar uniform between 0 and 2:
            % x=U_RV(1,1); 
            % vector (row or column):
            % x=U_RV([10,1]]): vector of indep. unif distr. RVs 

                obj=unc_t(mu);

                global MCSAMPLES

                if isscalar(intv)
                    intv=intv*ones(1,length(mu));
                end

                for i=1:length(mu)
                        obj(i).values=intv(i)*(2*rand(1,MCSAMPLES)-1)+mu(i);
                
                end

          end