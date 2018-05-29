# TODO: Add comment
# 
# Author: Ted under credit of Young Lee
###############################################################################

for(prepare_const in 1:1){
	
	GMR=0.95;	
	delta_true=0;	delta_nomial=delta_true
	#+0.05*rnorm(1);
	
	Period_dif=7e-1;				#	futility-para_1
	alpha=0.10;						#	default CI-para
	theta=0.80;						#	tolerance para
	U=-log(theta);	L=log(theta);	#	tolerance border
	
	sigma_square_true=0.25;			#	major para
	gamma=1.1;						#	rate of est sigma_s to real sigma_s
	sigma_nomial=gamma*sigma_square_true;
	
	beta=0.80;						#	power
	lambda=-0.1;					#	assist para
	pilot_power=0.5;	
	pilot_size=0.3;					#	pilot-magnitude
	pilot_dtr=0.5;
	pivotal_dtr=0.5;
};rm(prepare_const);
for(easy_func in c(1)){
	CV_W<-function(s_square_local=0){sqrt(exp(s_square_local)-1);}
	
	CV_W_Inverse<-function(cv_local=0.5){log(cv_local^2+1);}
	
	
	
};rm(easy_func);


###############################################################################
for(Wang_func_component in c(1)){
	
	Empherical_power_core<-function(lamb=0, delta=delta_nomial, sigma=gamma*sigma_square_true,
			sample_size=12, lambda_alg="TRUE"){
		#	coef_loc = sqrt(n_11^{-1}+n_12^{-1});
		#	lambda is one alg-form
		n_local=max(as.integer(sample_size), 24);
		n1_local=max((n_local*pilot_dtr), 12);
		t_df=n_local-2;
		
		coef_loc=1/n1_local + 1/(n_local-n1_local);
		
		if(lambda_alg=="TRUE"){
			return(integrate(function(x1){
								a=pnorm( (U-delta)/sqrt(sigma*coef_loc/2)-
												qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) )-
												lamb*sqrt(x1/( t_df*gamma*coef_loc)) );
								b=pnorm( (L-delta)/sqrt(sigma*coef_loc/2)+
												qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) )+
												lamb*sqrt(x1/( t_df*gamma*coef_loc)) );
								return( (a-b)*(exp(-x1/2)*x1^(t_df/2-1)*1/( 2^(t_df/2)*gamma(t_df/2))) );
							}, lower=sample_size/2, upper=sample_size*2) );
		}
		if(lambda_alg=="FALSE"){
			return(integrate(function(x1){
								a=pnorm( (U-delta-lamb)/sqrt(sigma*coef_loc/2)-
												qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) ) );
								b=pnorm( (L-delta-lamb)/sqrt(sigma*coef_loc/2)+
												qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) ));
								return( (a-b)*(exp(-x1/2)*x1^(t_df/2-1)*1/( 2^(t_df/2)*gamma(t_df/2))) );
							}, lower=sample_size/2, upper=sample_size*2) );
		}
	}
	
	P_uniroot<-function(f=function(x){return(x);}, beta_assign=beta){
		jump_s=1;
		N_L=2*sigma_square_true*(qnorm(1-alpha/2)+qnorm(0.5+beta/2) )^2/U^2;
		
		temp_loc=f(N_L);
		if(temp_loc>=beta_assign){
			print("not available;too weak");
			return(N_L);
		}
		N_U=400.0;
		counter=0;
		while(jump_s>0 & counter<200){
			N_L=N_L+2;
			temp_loc=f(N_L);
			if(temp_loc>beta_loc){
				jump_s=-1;
			}
			if(N_L>=N_U){N_L=N_U;jump_s=-1;}
			counter=counter+1;
		}
		return(N_L);
	}
	
	lambda_alg_Wang<-function(x){
		test_my<-function(y){gamma(y/2+0.5)/gamma(y/2);}
		
		n_star=as.integer(x*pilot_size);
		n1=as.integer(x*pivotal_dtr);n2=x-n1;
		n1_star=as.integer(n_star*pilot_dtr);n2_star=n_star-n1_star;
		
		coef=(test_my(x-2)/test_my( n_star-2 ))*sqrt((n_star-2)/(x-2));
		lambda_algf= coef*qt(p=(2-alpha)/2, df=x-2) * sqrt( 1.0/n1+1.0/n2 ) - qt(p=(2-alpha)/2.0, df=n_star-2)*sqrt( 1.0/n1_star+1.0/n2_star);
		return(lambda_algf);
	}
	lambda_alg_Wang(180)
	
	test_stat<-function(n1_local=10, n_local=50, delta_local=0.0, Period_local=1e-1, s_square_local=5e-1, oringinal="Yes"){
		n2_local=n_local-n1_local;
		if(n2_local>4){
			seq1_local=rnorm(n=n1_local, mean=(delta_local+Period_local)/2.0, sd=sqrt(s_square_local/2.0));
			seq2_local=rnorm(n=n2_local, mean=(delta_local-Period_local)/2.0, sd=sqrt(s_square_local/2.0));
			#	random-samples
			
			sigma_square_hat_local=( var(seq1_local)*(n1_local-1)+var(seq2_local)*(n2_local-1))/(n_local-2);
			delta_hat_local=mean(seq1_local)+mean(seq2_local);
			#	test-statistics
			if(oringinal=="Yes"){
				return(c(sigma_square_hat_local, delta_hat_local, seq1_local, seq2_local));
			}
			else{
				return(c(sigma_square_hat_local, delta_hat_local));
			}
		}
		else{
			print("Error in n2_star!!!Sample size-distri needs to be checked");
			rm(n2_local);
			return(NA);
		}
	}
	
};rm(Wang_func_component);

Wang_Simulation<-function(method="Wang", iter=5000, operation="not_strict"){
	iter_inside=as.integer(iter);
	true_counter=0;
	if(method=="Common_use"){
		for(index in 1:iter_inside){
			X=test_stat(n1_local=n1_star, n_local=n_star, delta_local=delta_true, 
					Period_local=Period_dif, s_square_local=sigma_square_true, oringinal = "No");
			delta_hat=as.numeric(X[2]);sigma_s_hat=as.numeric(X[1]);
			d_L=delta_hat-sqrt(sigma_square_true*(1.0/(n1_star)+1.0/(n2_star) )/2)*qnorm(p=(2-alpha)/2);
			#d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1_star+1.0/n2_star))*qt((2-alpha)/2, df=n_star-2);
			d_U=delta_hat*2-d_L;
			
			if( ( (d_L<U | d_U>L) & operation=="not_strict") | ( (d_U<U & d_L>L) & operation=="strict" ) | (operation=="Pass") ){
				X=test_stat(n1_local=n1, n_local=n, delta_local=delta_true, 
						Period_local=Period_dif, s_square_local=sigma_square_true, oringinal = "No");
				delta_hat=as.numeric(X[2]);sigma_s_hat=as.numeric(X[1]);
				
				#d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1+1.0/n2) )*qt(p=(2-alpha)/2, df=n-2);
				d_L=delta_hat-sqrt(sigma_square_true*(1.0/(n1)+1.0/(n2) )/2)*qnorm(p=(2-alpha)/2);
				d_U=delta_hat*2-d_L;
				if(d_L>L & d_U<U){
					true_counter=true_counter+1;
				}
			}
		}
		return(true_counter/iter_inside);
	}
	
	if(method=="Wang"){
		return(Wang_method(iter=iter, operations=operation));
	}
	if(method=="Wang_pool"){
		return(Wang_method_pool(iter=iter, operations=operation));
	}
	if(method=="Wang_old"){
		return(Wang_method_old(iter=iter, operations=operation));
	}
	if(method=="Dedicated"){
		for(index in 1:iter_inside){
			X=test_stat(n1_local=n1_star, n_local=n_star, delta_local=delta_true, 
					Period_local=Period_dif, s_square_local=sigma_square_true, oringinal = "No");
			delta_hat=as.numeric(X[2]);sigma_s_hat=as.numeric(X[1]);
			d_L=delta_hat-sqrt(sigma_square_true*(1.0/(n1_star)+1.0/(n2_star) )/2)*qnorm(p=(2-alpha)/2);
			#d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1_star+1.0/n2_star))*qt((2-alpha)/2, df=n_star-2);
			d_U=delta_hat*2-d_L;
			
			if( ( (d_L<U | d_U>L) & operation=="not_strict") | ( (d_U<U-lambda & d_L>L+lambda) & operation=="strict" ) | (operation=="Pass") ){
				X=test_stat(n1_local=n1, n_local=n, delta_local=delta_true, 
						Period_local=Period_dif, s_square_local=sigma_square_true, oringinal = "No");
				delta_hat=as.numeric(X[2]);sigma_s_hat=as.numeric(X[1]);
				
				#d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1+1.0/n2) )*qt(p=(2-alpha)/2, df=n-2);
				d_L=delta_hat-sqrt(sigma_square_true*(1.0/(n1)+1.0/(n2) )/2)*qnorm(p=(2-alpha)/2);
				d_U=delta_hat*2-d_L;
				if(d_L>L & d_U<U){
					true_counter=true_counter+1;
				}
			}
		}
		return(true_counter/iter_inside);
	}
	
}

Wang_Simulation()



Wang_method<-function(iter=5000, operations="not_strict"){
	iter_ahead=as.integer(iter);
	true_counter=0;
	
	lambda_alg_Wang_loc<-function(x){
		test_my<-function(y){gamma(y/2+0.5)/gamma(y/2);}
		
		n_star=as.integer(x*pilot_size);
		n1=as.integer(n*pivotal_dtr);			n2=n-n1;
		n1_star=as.integer(n_star*pilot_dtr);	n2_star=n_star-n1_star;
		
		coef=(test_my(x-2)/test_my( n_star-2 ))*sqrt((n_star-2)/(x-2));
		lambda_algf= coef*qt(p=(2-alpha)/2, df=x-2) * sqrt(1.0/n1+1.0/n2) - qt(p=(2-alpha)/2.0, df=n_star-2)*sqrt(1.0/n1_star+1.0/n2_star);
		return(lambda_algf);
	}
	
	Empherical_power_size_loc<-function(x){
		lambda_locs=lambda_alg_Wang_loc(x);
		return( (Empherical_power_core(lamb = lambs, sample_size = as.integer(x*pilot_size) ) )$value );
	}
	
	n0=P_uniroot(f=Empherical_power_size_loc, beta_assign = beta);
	print(n0);
	
	n_Wang=n0;	n_Wang=as.integer(n_Wang/2.0)*2+2;
	n1_Wang=as.integer(n_Wang/2.0);	n2_Wang=n_Wang-n1_Wang;
	
	n_Wang_star=as.integer(max(n_Wang*pilot_size/2.0, 5.0))*2+2;
	n1_Wang_star=as.integer(n_Wang_star/2.0);
	n2_Wang_star=n_Wang_star-n1_Wang_star;
	
	test_my<-function(x){
		gamma(x/2+0.5)/gamma(x/2);
	}
	for(index in 1:iter_ahead){
		X=test_stat(n1_local=n1_Wang_star, n_local=n_Wang_star, delta_local=delta_true, 
				Period_local=Period_dif, s_square_local=sigma_square_true);
		delta_hat=as.numeric(X[2]);sigma_s_hat=as.numeric(X[1]);
		
		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1_Wang_star+1.0/n2_Wang_star))*qt(p=(2-alpha)/2.0, df=n_Wang_star-2);
		d_U=2*delta_hat-d_L;
		
		coef=(test_my(n_Wang-2)/test_my(n_Wang_star-2))*sqrt((n_Wang_star-2)/(n_Wang-2));
		#coef=1;	
		lambda=(coef*qt(p=(2-alpha)/2, df=n_Wang-2)*sqrt(1.0/n1_Wang+1.0/n2_Wang)-qt(p=(2-alpha)/2.0, df=n_Wang_star-2)*sqrt(1.0/n1_Wang_star+1.0/n2_Wang_star) )*sqrt(sigma_s_hat);
		
		if(lambda>0){print(index);}
		if(index==2){print(lambda);}
		if( ( (d_L<U-lambda | d_U>L+lambda) & operations=="not_strict") | ( (d_U<U-lambda | d_L>L+lambda) & operations=="strict" ) ){
			X=test_stat(n1_local=n1_Wang, n_local=n_Wang, delta_local=delta_true,
					Period_local=Period_dif, s_square_local=sigma_square_true);
			delta_hat=as.numeric(X[2]);	sigma_s_hat=as.numeric(X[1]);
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1_Wang+1.0/n2_Wang))*qt((2-alpha)/2, df=n_Wang-2);
			d_U=delta_hat*2-d_L;
			if(d_L>L & d_U<U){
				true_counter=true_counter+1;
			}
		}
	}
	return(true_counter*1.0/iter_ahead);
}

Wang_method_old<-function(iter=5000, operations="not_strict"){
	iter_ahead=as.integer(iter);
	true_counter=0;
	
	lambda_alg_Wang_loc<-function(x){
		test_my<-function(y){gamma(y/2+0.5)/gamma(y/2);}
		
		n_star=as.integer(x*pilot_size);
		n1=as.integer(n*pivotal_dtr);			n2=n-n1;
		n1_star=as.integer(n_star*pilot_dtr);	n2_star=n_star-n1_star;
		
		#coef=(test_my(x-2)/test_my( n_star-2 ))*sqrt((n_star-2)/(x-2));
		coef=1;
		lambda_algf= coef*qt(p=(2-alpha)/2, df=x-2) * sqrt(1.0/n1+1.0/n2) - qt(p=(2-alpha)/2.0, df=n_star-2)*sqrt(1.0/n1_star+1.0/n2_star);
		return(lambda_algf);
	}
	
	Empherical_power_size_loc<-function(x){
		lambda_locs=lambda_alg_Wang_loc(x);
		return( (Empherical_power_core(lamb = lambda_locs, sample_size = as.integer(x*pilot_size) ) )$value );
	}
	
	n0=P_uniroot(f=Empherical_power_size_loc, beta_loc = beta);
	print(n0);
	
	n_Wang=n0;
	n_Wang=as.integer(n_Wang/2.0)*2+2;
	n_Wang=max(n_Wang, 16);
	n1_Wang=as.integer(n_Wang/2.0);
	n2_Wang=n_Wang-n1_Wang;
	
	n_Wang_star=as.integer(max(n_Wang*pilot_size/2.0, 5.0))*2+2;
	n1_Wang_star=as.integer(n_Wang_star/2.0);
	n2_Wang_star=n_Wang_star-n1_Wang_star;
	
	test_my<-function(x){
		gamma(x/2+0.5)/gamma(x/2);
	}
	for(index in 1:iter_ahead){
		X=test_stat(n1_local=n1_Wang_star, n_local=n_Wang_star, delta_local=delta_true, 
				Period_local=Period_dif, s_square_local=sigma_square_true);
		delta_hat=as.numeric(X[2]);sigma_s_hat=as.numeric(X[1]);
		
		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1_Wang_star+1.0/n2_Wang_star))*qt(p=(2-alpha)/2.0, df=n_Wang_star-2);
		d_U=2*delta_hat-d_L;
		
		#coef=(test_my(n_Wang-2)/test_my(n_Wang_star-2))*sqrt((n_Wang_star-2)/(n_Wang-2));
		coef=1;	
		lambda=(coef*qt(p=(2-alpha)/2, df=n_Wang-2)*sqrt(1.0/n1_Wang+1.0/n2_Wang)-qt(p=(2-alpha)/2.0, df=n_Wang_star-2)*sqrt(1.0/n1_Wang_star+1.0/n2_Wang_star) )*sqrt(sigma_s_hat);
		#rm(X, delta_hat, sigma_s_hat);
		if(lambda>0){print(index);}
		if( ( (d_L<U-lambda | d_U>L+lambda) & operations=="not_strict") | ( (d_U<U-lambda & d_L>L+lambda) & operations=="strict" ) ){
			X=test_stat(n1_local=n1_Wang, n_local=n_Wang, delta_local=delta_true,
					Period_local=Period_dif, s_square_local=sigma_square_true);
			delta_hat=as.numeric(X[2]);sigma_s_hat=as.numeric(X[1]);
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1_Wang+1.0/n2_Wang))*qt((2-alpha)/2, df=n_Wang-2);
			d_U=delta_hat+sqrt(sigma_s_hat*(1.0/n1_Wang+1.0/n2_Wang))*qt((2-alpha)/2, df=n_Wang-2);
			rm(X);
			if(d_L>L & d_U<U){
				true_counter=true_counter+1;
			}
		}
	}
	
	rm(index, n_Wang, n_Wang_star, n1_Wang, n1_Wang_star, n2_Wang, n2_Wang_star)
	seq1=NA;seq2=NA;
	return(true_counter*1.0/iter_ahead);
}

Wang_method_pool<-function(iter=5000, operations="not_strict"){
	iter_ahead=as.integer(iter);
	true_counter=0;
	
	lambda_alg_Wang_loc<-function(x){
		test_my<-function(y){gamma(y/2+0.5)/gamma(y/2);}
		
		n_star=as.integer(x*pilot_size);
		n1=as.integer(n*pivotal_dtr);			n2=n-n1;
		n1_star=as.integer(n_star*pilot_dtr);	n2_star=n_star-n1_star;
		
		coef=(test_my(x-2)/test_my( n_star-2 ))*sqrt((n_star-2)/(x-2));
		lambda_algf= coef*qt(p=(2-alpha)/2, df=x-2) * sqrt(1.0/n1+1.0/n2) - qt(p=(2-alpha)/2.0, df=n_star-2)*sqrt(1.0/n1_star+1.0/n2_star);
		return(lambda_algf);
	}
	
	Empherical_power_size_loc<-function(x){
		lambda_locs=lambda_alg_Wang_loc(x);
		return( (Empherical_power_core(lamb = lambda_locs, sample_size = as.integer(x*pilot_size) ) )$value );
	}
	
	n0=P_uniroot(f=Empherical_power_size_loc, beta_loc = beta);
	print(n0);
	print(Empherical_power_size_loc(n0));
	
	n_Wang=n0;	n_Wang=as.integer(n_Wang/2.0)*2+2;
	n_Wang=max(n_Wang, 16);
	n1_Wang=as.integer(n_Wang/2.0);
	n2_Wang=n_Wang-n1_Wang;
	
	n_Wang_star=as.integer(max(n_Wang*pilot_size/2.0, 5.0))*2+2;
	n1_Wang_star=as.integer(n_Wang_star/2.0);
	n2_Wang_star=n_Wang_star-n1_Wang_star;
	
	test_my<-function(x){
		gamma(x/2+0.5)/gamma(x/2);
	}
	
	for(index in 1:iter_ahead){
		seq1=rnorm(n=n1_Wang_star, mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0));
		seq2=rnorm(n=n2_Wang_star, mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0));
		
		sigma_s_hat=( var(seq1)*(n1_Wang_star-1)+var(seq2)*(n2_Wang_star-1))/(n_Wang_star-2);
		
		delta_hat=mean(seq1)+mean(seq2);
		#	test-statistics
		
		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/n1_Wang_star + 1.0/n2_Wang_star))*qt(p=(2-alpha)/2.0, df=n_Wang_star-2);
		d_U=2*delta_hat-d_L;
		
		coef=test_my(n_Wang-2)/test_my(n_Wang_star-2);
		#coef=1;
		lambda=(coef*qt(p=(2-alpha)/2, df=n_Wang-2)*sqrt(1.0/n1_Wang+1.0/n2_Wang)-qt(p=(2-alpha)/2.0, df=n_Wang_star-2)*sqrt(1.0/n1_Wang_star+1.0/n2_Wang_star) )*sqrt(sigma_s_hat);
		#rm(X, delta_hat, sigma_s_hat);
		if(lambda>=0){print(index);}
		
		if( ( (d_L<U-lambda | d_U>L+lambda) & operations=="not_strict") | ( (d_U<U-lambda & d_L>L+lambda) & operations=="strict" ) ){
			seq1=c(seq1, rnorm(n=n1_Wang, mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)) );
			seq2=c(seq2, rnorm(n=n2_Wang, mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)) );
			
			#	random-samples
			
			sigma_s_hat=( var(seq1)*(n1_Wang_star+n1_Wang-1)+var(seq2)*(n2_Wang_star+n2_Wang-1))/(n_Wang+n_Wang_star-2);
			delta_hat=mean(seq1)+mean(seq2);
			
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/(n1_Wang+n1_Wang_star)+1.0/(n2_Wang_star+n2_Wang) ))*qt((2-alpha)/2, df=n_Wang+n_Wang_star-2);
			d_U=2*delta_hat-d_L;
			if(d_L>L & d_U<U){
				true_counter=true_counter+1;
			}
		}
	}
	return(true_counter*1.0/iter_ahead);
}

pilot_size=0.4
delta_nomial
delta_true=U

Wang_Simulation(method="Common_use", iter=100000, operation = "strict")
Wang_Simulation(method="Wang", iter=40000, operation = "strict")
Wang_Simulation(method="Wang_pool", iter=40000, operation = "strict")
Wang_Simulation(method="Wang_old", iter=40000, operation = "strict")
Wang_Simulation(method="Dedicated", iter=40000, operation = "strict")

Wang_Simulation(method="Common_use", iter=40000, operation = "not_strict")
Wang_Simulation(method="Wang", iter=40000, operation = "not_strict")
Wang_Simulation(method="Wang_pool", iter=40000, operation = "not_strict")
Wang_Simulation(method="Wang_old", iter=40000, operation = "not_strict")

#	Common_use		

beta=0.8
delta_true=log(GMR)
delta_true=0

change_CV=0.5
CV_W_Inverse(change_CV)
for(modify_procedure_1 in c(1)){
	#sigma_square_true=CV_W_Inverse(change_CV);
	
	n0=12
	for(indexs in 1:10){
		record=n0;
		n0=as.integer(2*sigma_nomial*( qnorm(p=1-alpha/2)+qnorm(p=0.5+beta/2) )^2 * 1/(U^2) );
		if(abs(record-n0)<=1){break;}
	};if(indexs==10){print("potential under estimmate");}
	rm(indexs);
	#n0=144
	n=n0;	n=as.integer(n/2.0)*2+2;
	n1=as.integer(n/2.0);	n2=n-n1;
	
	n_star=as.integer(max(n*pilot_size/2.0, 6.0))*2;
	n1_star=as.integer(n_star/2.0);	n2_star=n_star-n1_star;
	
	print(n);
	print(n_star);
};rm(modify_procedure_1);

pilot_size=0.3



beta=0.8
delta_true=log(GMR)
change_CV=0.5
CV_W_Inverse(change_CV)
for(modify_procedure in c(1)){
	sigma_square_true=CV_W_Inverse(change_CV);
	
	Empherical_power_size_D<-function(x){
		return( (Empherical_power_core(lamb = lambda, sample_size = as.integer(x*pilot_size) ) )$value - 0.80);
	}
	
	
	n0=P_uniroot(f=Empherical_power_size_D, beta_loc = beta);
	#rm(lambda_alg_Wang, Empherical_power_size)
	
	n0=60;
	n=n0;	n=as.integer(n/2.0)*2+2;
	n1=as.integer(n/2.0);	n2=n-n1;
	
	n_star=as.integer(max(n*pilot_size/2.0, 6.0))*2;
	n1_star=as.integer(n_star/2.0);	n2_star=n_star-n1_star;
	
	print(n);
	print(n_star);
};rm(modify_procedure);





