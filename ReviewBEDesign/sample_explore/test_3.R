# TODO: Add comment
# 
# Author: Ted under credit of Young Lee
###############################################################################



Empherical_power_size_loc<-function(x){
	lambda_locs=lambda_alg_Wang(x);
	return( (Empherical_power_core(lambda_loc = lambda_locs, sample_size = as.integer(x*pilot_size) ) )$value );
}
Empherical_power_size_loc(50)




















U
gamma
sigma_square_true

Empherical_power_core<-function(lambda_loc=0, delta_nomials=delta_nomial, 
		sigma_nomial=gamma*sigma_square_true, sample_size=300){
	#	coef_loc = sqrt(n_11^{-1}+n_12^{-1});
	#	lambda is one alg-form
	sample_size=200
	
	n_local=max(as.integer(sample_size), 24);
	n1_local=max((n_local*pilot_dtr), 12);
	t_df=n_local-2;
	
	coef_loc=1/n1_local + 1/(n_local-n1_local);
	return(integrate(function(x1){
						a=pnorm( (U-delta_nomials)/sqrt(sigma_nomial*coef_loc/2)-
										qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) )-
										lambda_loc*sqrt(x1/( t_df*gamma*coef_loc)) );
						
						b=pnorm( (L-delta_nomials)/sqrt(sigma_nomial*coef_loc/2)+
										qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) )+
										lambda_loc*sqrt(x1/( t_df*gamma*coef_loc)) );
						
						return( (a-b)*(exp(-x1/2)*x1^(t_df/2-1)*1/( 2^(t_df/2)*gamma(t_df/2))) );
					}, lower=0, upper=Inf) );
}




Empherical_power_core(sample_size = 10)







test_free<-function(n){
	print(1-1/sqrt(n));
	return(pchisq(q=0.25*n^2, df=n));
}

test_free(100)




























