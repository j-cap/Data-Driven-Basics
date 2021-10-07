function obj = EXP_RV(mu)
            % Exponentially distributed samples
            % examples:
            % scalar:
            % x=EXP_RV(1);
            % vector (row or column):
            % x=EXP_RV([10,1]]): vector of expontially distr. RVs 
                global MCSAMPLES
                obj=unc_t(mu);

                for i=1:length(mu)
                        obj(i).values=exprnd(mu(i),1,MCSAMPLES);
                        
                end
  
          end