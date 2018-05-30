# TODO: Add comment
# 
# Author: Ted under credit of Young Lee
###############################################################################

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


for(Wang_func_component in c(1)){
	
	Empherical_power_core<-function(lambda_loc=0, delta_nomials=delta_nomial, 
			sigma_nomial=gamma*sigma_square_true, sample_size=12, lambda_alg="TRUE"){
		#	coef_loc = sqrt(n_11^{-1}+n_12^{-1});
		#	lambda is one alg-form
		n_local=max(as.integer(sample_size), 24);
		n1_local=max((n_local*pilot_dtr), 12);
		t_df=n_local-2;
		
		coef_loc=1/n1_local + 1/(n_local-n1_local);
		
		if(lambda_alg=="TRUE"){
			return(integrate(function(x1){
								a=pnorm( (U-delta_nomials)/sqrt(sigma_nomial*coef_loc/2)-
												qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) )-
												lambda_loc*sqrt(x1/( t_df*gamma*coef_loc)) );
								b=pnorm( (L-delta_nomials)/sqrt(sigma_nomial*coef_loc/2)+
												qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) )+
												lambda_loc*sqrt(x1/( t_df*gamma*coef_loc)) );
								return( (a-b)*(exp(-x1/2)*x1^(t_df/2-1)*1/( 2^(t_df/2)*gamma(t_df/2))) );
							}, lower=sample_size/2, upper=sample_size*2) );
		}
		if(lambda_alg=="FALSE"){
			return(integrate(function(x1){
								a=pnorm( (U-delta_nomials-lambda_loc)/sqrt(sigma_nomial*coef_loc/2)-
												qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) ) );
								b=pnorm( (L-delta_nomials-lambda_loc)/sqrt(sigma_nomial*coef_loc/2)+
												qt(p=1-alpha/2, df=t_df)*sqrt(x1/(gamma*t_df) ));
								return( (a-b)*(exp(-x1/2)*x1^(t_df/2-1)*1/( 2^(t_df/2)*gamma(t_df/2))) );
							}, lower=sample_size/2, upper=sample_size*2) );
		}
	}
	
	X1=Empherical_power_core(lambda_loc=-0.2, sample_size = 60, lambda_alg = "FALSE")
	
	P_uniroot<-function(f=function(x){return(x);}, beta_loc=beta){
		jump_s=1;
		N_L=2*sigma_square_true*(qnorm(1-alpha/2)+qnorm(0.5+beta/2) )^2/U^2;
		
		temp_loc=f(N_L);
		if(temp_loc>=beta_loc){
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
	
	
};rm(Wang_func_component);




###	Wang_Method
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
		return( (Empherical_power_core(lambda_loc = lambda_locs, sample_size = as.integer(x*pilot_size) ) )$value );
	}
	
	n0=P_uniroot(f=Empherical_power_size_loc, beta_loc = beta);
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
		return( (Empherical_power_core(lambda_loc = lambda_locs, sample_size = as.integer(x*pilot_size) ) )$value );
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
		return( (Empherical_power_core(lambda_loc = lambda_locs, sample_size = as.integer(x*pilot_size) ) )$value );
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
		return( (Empherical_power_core(lambda_loc = lambda, sample_size = as.integer(x*pilot_size) ) )$value - 0.80);
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




test_mine<-function(n){
	alpha_lo=0.05;
	return(qnorm(1-alpha_lo)-qt(1-alpha_lo, df=n));
}





power_local<-function(choice="N", sigma=0.5, mean_nomial=1, N1_loc=5, N2_loc=5, alpha_loc=0.05, lambda_loc=0){
	lambda_loc=min(lambda_loc, 0);
	if(choice=="N"){
		return( pnorm( (U-mean_nomial-lambda_loc)/(sqrt(sigma*(1/N1_loc+1/N2_loc) )) - qnorm(1-alpha_loc/2))-pnorm( (L-mean_nomial+lambda_loc)/(sigma*sqrt(1/N1_loc+1/N2_loc))+qnorm(1-alpha_loc/2)));
	}
	
	if(choice=="T"){
		return( pnorm( (U-mean_nomial-lambda_loc)/(sqrt(sigma*(1/N1_loc+1/N2_loc) )) - qt(p=1-alpha_loc/2, df=N1_loc+N2_loc-2) )-pnorm( (L-mean_nomial+lambda_loc)/(sigma*sqrt(1/N1_loc+1/N2_loc))+qt(p=1-alpha_loc/2, df=N1_loc+N2_loc-2)));
	}
}


Fuglsang_Simulation<-function(Order=c(1000, "Yes", "Yes", "No", "No"), N1_dtr=c(1, 0.3), N2_dtr=0.5){
	iter_inside=max(as.integer(Order[1]), 1000);
	pilot_permission=Order[2];
	pilot_veto=Order[3];
	GMR_obs=Order[4];
	CV_obs=Order[5];
	
	CV_sigma=gamma*sigma_square_true;
	delta_est=log(GMR);
	
	true_counter=0;
	pilot_indicator=-1;
	N=N1_dtr[1]*12;
	N1=as.integer(N*N1_dtr[2]);N2=N-N1;
	for(index in 1:iter_inside){
		X=test_stat(n1_local=N1, n_local = N, delta_local = delta_true,
				Period_local = Period_dif, s_square_local = sigma_square_true);
		
		delta_hat=X[2];	sigma_s_hat=X[1]
		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha)/2.0, df=N-2);
		d_U=2*delta_hat-d_L;
		
		if(d_L > L & d_U<U){pilot_indicator=10;}
		if(d_L > U | d_U<L){pilot_indicator=-10;}
		
		if( (pilot_veto=="No" & pilot_indicator<0) | (pilot_indicator==0) | (pilot_indicator>0 & pilot_permission=="No") ){
			if(CV_obs=="Yes"){CV_sigma=2*sigma_s_hat;}
			if(GMR_obs=="Yes"){delta_est=delta_hat;}
			
			for(prepare_again in c(1)){
				Empherical_power_size_loc<-function(x){
					return(Empherical_power_core(lambda_loc = 0, sample_size = as.integer(x*pilot_size))$value);
				}
				N=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
			};rm(prepare_again);
			
			N1=as.integer(N*N2_dtr);
			N2=N-N1;
			X=test_stat(n1_local=N1, n_local = N, delta_local = delta_true,
					Period_local = Period_dif, s_square_local = sigma_square_true);
			
			delta_hat=X[2];sigma_s_hat=X[1]
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha)/2.0, df=N-2);
			d_U=2*delta_hat-d_L;
			if(d_L>L & d_U<U){
				true_counter=true_counter+1;
			}
		}
		if(pilot_indicator>0 & pilot_permission=="Yes"){
			true_counter=true_counter+1;
		}
	}
	return(true_counter/iter_inside);
}

delta_true=log(GMR)
delta_true=U


#	pilot_permission=Order[2];	pilot_veto=Order[3];	GMR_obs=Order[4];	CV_obs=Order[5];
Fuglsang_Simulation(Order=c(10000, "No", "Yes", "Yes", "Yes"), N2_dtr = 0.3, N1_dtr =c(5,0.3));










Potvin_Simulation<-function(method="A", iter=5000, N1_dtr=0.3, N2_dtr=0.3){
	if(method=="A"){
		return(Potvin_method_A(iter=iter, N1_dtr = N1_dtr, N2_dtr = N2_dtr));
	}
	
	if(method=="B"){
		return(Potvin_method_B(iter=iter, N1_dtr = N1_dtr, N2_dtr = N2_dtr));
	}
	
	if(method=="C"){
		return(Potvin_method_C(iter=iter, N1_dtr = N1_dtr, N2_dtr = N2_dtr));
	}
	
	if(method=="D"){
		return(Potvin_method_D(iter=iter, N1_dtr = N1_dtr, N2_dtr = N2_dtr));
	}
	
	else{
		print("method alg error!");
		return(NA);
	}
	
}

Potvin_method_A<-function(iter=5000, N1_dtr=0.3, N2_dtr=0.3){
	iter_local=as.integer(iter);
	alpha_local=0.05;
	true_counter=0;
	
	for(prepare_again in c(1)){
		Empherical_power_size_loc<-function(x){
			return(Empherical_power_core(sample_size = as.integer(x*pilot_size))$value);
		}
		N_fix=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
	};rm(prepare_again);
	
	print(N_fix);
	
	N=N_fix;N1=as.integer(N*N1_dtr);N2=N-N1;
	for(index in 1:iter_local){
		seq1_local=c( rnorm(n=N1, mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0) ) );
		seq2_local=c( rnorm(n=N2, mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0) ) );
		delta_hat=mean(seq1_local)+mean(seq2_local);
		sigma_s_hat=( var(seq1_local)*(length(seq1_local)-1)+var(seq2_local)*(length(seq2_local)-1))/(length(seq1_local)+length(seq2_local)-2);
		
		if(power_local(sigma=sigma_s_hat, mean_nomial=log(GMR), N1_loc=N1, N2_loc=N2, alpha_loc=alpha_local)>=0.8){
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
			d_U=delta_hat*2-d_L;
			if(d_L>L & d_U<U){true_counter=true_counter+1;}
		}
		else{
			N1_old=N1;N2_old=N2;
			for(prepare_again in c(1)){
				Empherical_power_size_loc<-function(x){
					return(Empherical_power_core(sample_size = as.integer(x*pilot_size), sigma_nomial = sigma_s_hat )$value);
				}
				N=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
			};rm(prepare_again);
			#N=as.integer(2*sigma_s_hat*(qnorm(0.5+beta/2)+qnorm(1-alpha_local/2))^2*1.0/U^2);
			N1=as.integer(N*N2_dtr);N2=N-N1;
			
			seq1_local=c(seq1_local, rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
			seq2_local=c(seq2_local, rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
			delta_hat=mean(seq1_local)+mean(seq2_local);
			sigma_s_hat=( var(seq1_local)*(length(seq1_local)-1)+var(seq2_local)*(length(seq2_local)-1))/(length(seq1_local)+length(seq2_local)-2);
			
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
			d_U=delta_hat*2-d_L;
			if(d_L>L*1.0 & d_U<U*1.0){true_counter=true_counter+1;}
		}
		rm(seq1_local, seq2_local);
	}
	return(true_counter/iter_local);
}
Potvin_method_B<-function(iter=5000, N1_dtr=0.3, N2_dtr=0.3){
	iter_local=as.integer(iter);
	alpha_local=0.05;
	true_counter=0;
	
	alpha_local=0.0294;
	
	for(prepare_again in c(1)){
		Empherical_power_size_loc<-function(x){
			return(Empherical_power_core(sample_size = as.integer(x*pilot_size))$value);
		}
		n0=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
	};rm(prepare_again);
	
	print(N_fix);
	
	N_fix=n0;N=N_fix;N1=as.integer(N*N1_dtr);N2=N-N1;	
	for(index in 1:iter_local){
		seq1=c( rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		seq2=c( rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		delta_hat=mean(seq1)+mean(seq2);
		sigma_s_hat=( var(seq1)*(length(seq1)-1)+var(seq2)*(length(seq2)-1))/(length(seq1)+length(seq2)-2);
		
		
		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
		d_U=delta_hat*2-d_L;
		
		if(d_L>L*1 & d_U<U){true_counter=true_counter+1;}
		else{
			if(power_local(sigma=sigma_s_hat, mean_nomial=log(GMR), N1_loc=N1, N2_loc=N2, alpha_loc=alpha_local)<0.8){
				for(prepare_again in c(1)){
					Empherical_power_size_loc<-function(x){
						return(Empherical_power_core(sample_size = as.integer(x*pilot_size), sigma_nomial = sigma_s_hat )$value);
					}
					N=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
				};rm(prepare_again);
				#N=as.integer(2*sigma_s_hat*(qnorm(0.5+beta/2)+qnorm(1-alpha_local/2))^2*1.0/U^2);
				N1=as.integer(N*N2_dtr);N2=N-N1;
				
				seq1=c(seq1, rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
				seq2=c(seq2, rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
				delta_hat=mean(seq1)+mean(seq2);
				sigma_s_hat=( var(seq1)*(length(seq1)-1)+var(seq2)*(length(seq2)-1))/(length(seq1)+length(seq2)-2);
				
				d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
				d_U=delta_hat*2-d_L;
				if(is.na(delta_hat)){print(delta_hat);}
				if(d_L>L & d_U<U*1){true_counter=true_counter+1;}
				rm(seq1, seq2);
			}
		}
	}
	return(true_counter/iter_local);
}
Potvin_method_C<-function(iter=5000, N1_dtr=0.3, N2_dtr=0.3){
	iter_local=as.integer(iter);
	alpha_local=0.05;
	true_counter=0;
	
	alpha_local=0.05;
	for(prepare_again in c(1)){
		Empherical_power_size_loc<-function(x){
			return(Empherical_power_core(sample_size = as.integer(x*pilot_size))$value);
		}
		n0=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
	};rm(prepare_again);
	
	print(N_fix);
	
	N_fix=n0;N=N_fix;N1=as.integer(N*N1_dtr);N2=N-N1;
	for(index in 1:iter_local){
		seq1=c( rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		seq2=c( rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		delta_hat=mean(seq1)+mean(seq2);
		sigma_s_hat=( var(seq1)*(length(seq1)-1)+var(seq2)*(length(seq2)-1))/(length(seq1)+length(seq2)-2);
		
		if(power_local(sigma=sigma_s_hat, mean_nomial=log(GMR), N1_loc=N1, N2_loc=N2, alpha_loc=alpha_local)>=0.8){
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
			d_U=delta_hat*2-d_L;
			if(d_L>L & d_U<U){true_counter=true_counter+1;}
		}
		else{
			alpha_local=0.0294;
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
			d_U=delta_hat*2-d_L;
			
			if(d_L>L & d_U<U){true_counter=true_counter+1;}
			else{
				for(prepare_again in c(1)){
					Empherical_power_size_loc<-function(x){
						return(Empherical_power_core(sample_size = as.integer(x*pilot_size), sigma_nomial = sigma_s_hat )$value);
					}
					N=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
				};rm(prepare_again);
				#N=as.integer(2*sigma_s_hat*(qnorm(0.5+beta/2)+qnorm(1-alpha_local/2))^2*1.0/U^2);
				N1=as.integer(N*N2_dtr);N2=N-N1;
				
				seq1=c(seq1, rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
				seq2=c(seq2, rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
				delta_hat=mean(seq1)+mean(seq2);
				sigma_s_hat=( var(seq1)*(length(seq1)-1)+var(seq2)*(length(seq2)-1))/(length(seq1)+length(seq2)-2);
				
				d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
				d_U=delta_hat*2-d_L;
				if(d_L>L & d_U<U){true_counter=true_counter+1;}
				rm(seq1, seq2);
			}
		}
	}
	
	return(true_counter/iter_local);
}
Potvin_method_D<-function(iter=5000, N1_dtr=0.3, N2_dtr=0.3){
	iter_local=as.integer(iter);
	alpha_local=0.05;
	true_counter=0;
	
	alpha_local=0.05;
	for(prepare_again in c(1)){
		Empherical_power_size_loc<-function(x){
			return(Empherical_power_core(sample_size = as.integer(x*pilot_size))$value);
		}
		n0=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
	};rm(prepare_again);
	
	print(N_fix);
	
	N_fix=n0;N=N_fix;N1=as.integer(N*N1_dtr);N2=N-N1;
	for(index in 1:iter_local){
		seq1=c( rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		seq2=c( rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		delta_hat=mean(seq1)+mean(seq2);
		sigma_s_hat=( var(seq1)*(length(seq1)-1)+var(seq2)*(length(seq2)-1))/(length(seq1)+length(seq2)-2);
		
		if(power_local(sigma=sigma_s_hat, mean_nomial=log(GMR), N1_loc=N1, N2_loc=N2, alpha_loc=alpha_local)>=0.8){
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
			d_U=delta_hat*2-d_L;
			if(d_L>L*1 & d_U<U){true_counter=true_counter+1;}
		}
		else{
			alpha_local=0.028;
			d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
			d_U=delta_hat*2-d_L;
			
			if(d_L>L & d_U<U*1){true_counter=true_counter+1;}
			else{
				for(prepare_again in c(1)){
					Empherical_power_size_loc<-function(x){
						return(Empherical_power_core(sample_size = x, sigma_nomial = sigma_s_hat )$value);
					}
					N=P_uniroot(f=Empherical_power_size_loc, beta_loc = 0.8);
				};rm(prepare_again);
				#N=as.integer(2*sigma_s_hat*(qnorm(0.5+beta/2)+qnorm(1-alpha_local/2))^2*1.0/U^2);
				N1=as.integer(N*N2_dtr);N2=N-N1;
				
				seq1=c(seq1, rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
				seq2=c(seq2, rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
				delta_hat=mean(seq1)+mean(seq2);
				sigma_s_hat=( var(seq1)*(length(seq1)-1)+var(seq2)*(length(seq2)-1))/(length(seq1)+length(seq2)-2);
				
				d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=(2-alpha_local)/2.0, df=N-2);
				d_U=delta_hat*2-d_L;
				if(d_L>L & d_U<U){true_counter=true_counter+1;}
				rm(seq1, seq2);
			}
		}
	}
	return(true_counter/iter_local);
}



pilot_power=0.8
pilot_power=0.3
delta_true=log(GMR)
delta_true=U
Potvin_Simulation(method="A", iter=10000)





test_tear<-function(x){
	print(x);
	sigma_nomial=0.111;
	
	print(sigma_nomial);
	return(x);
}
test_tear(11);print(sigma_nomial);



Psi_Simulation<-function(iter=5000, Psi=0.5, Max=200, N=12){
	
	iter_inside=as.integer(iter);
	
	#necessary yakso
	#here d_ij var=sigma^2
	alpha_Psi=0.10;
	sigma_nomial_Psi=2*sigma_nomial;
	delta_nomial_Psi=delta_nomial;
	
	sample_size_Psi<-function(n, Psi_l=0.5)
	{
		n1=as.integer(n*pivotal_dtr);			n2=n-n1;
		
		n_star=as.integer(n*pilot_size);	
		n1_star=as.integer(n_star*pilot_dtr);	n2_star=n_star-n1_star;
		
		time1=4/(1/n1_star+1/n2_star);
		time2=4/(1/(n1+n1_star)+1/(n2+n2_star) );
		time3=4/(1/n1+1/n2);
		
#		G1=qt(p=1-alpha/2, df=n_star-2)/sqrt(n_star-2);
#		a_U1=(U-delta_nomial)/sqrt( sigma_nomial/time1 );
#		a_L1=(L-delta_nomial)/sqrt( sigma_nomial/time1 );
#		r_U1=max(a_U1, (Psi_l*U-delta_nomial)/sqrt( sigma_nomial/time1 ) );
#		r_L1=min(a_L1, (Psi_l*L-delta_nomial)/sqrt( sigma_nomial/time1 ) );
		
		return(integrate( function(x){
							G1=qt(p=1-alpha/2, df=n_star-2)/sqrt(n_star-2);
							a_U1=(U-delta_nomial)/sqrt( sigma_nomial/time1 )-G1*sqrt(x);
							a_L1=(L-delta_nomial)/sqrt( sigma_nomial/time1 )+G1*sqrt(x);
							r_U1=max(a_U1, (Psi_l*U-delta_nomial)/sqrt( sigma_nomial/time1 )+G1*sqrt(x) );
							r_L1=min(a_L1, (Psi_l*L-delta_nomial)/sqrt( sigma_nomial/time1 )-G1*sqrt(x) );
							
							XX=pnorm(a_U1)-pnorm(a_L1)+
									integrate(function(x1){
												G2=qt(p=1-alpha/2, df=n+n_star-2);
												a_U=(U-delta_nomial)/sqrt( sigma_nomial/time2 )-G2;
												a_L=(L-delta_nomial)/sqrt( sigma_nomial/time2 )+G2;
												return(dnorm(x1)*(pnorm(a_U*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1) - 
																	pnorm(a_L*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1) ));
											},
											lower=r_L1, upper=a_L1 )$value+
									integrate(function(x1){
												G2=qt(p=1-alpha/2, df=n+n_star-2);
												a_U=(U-delta_nomial)/sqrt( sigma_nomial/time2 )-G2;
												a_L=(L-delta_nomial)/sqrt( sigma_nomial/time2 )+G2;
												return(dnorm(x1)*(pnorm(a_U*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1) - 
																	pnorm(a_L*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1) ));
											},
											lower=a_U1, upper=r_U1 )$value;
							return(XX*dchisq(x, df=n_star-2));},
						lower=n_star/2, upper=n_star*2)$value );
	}
	
	N0=1200;
	while(sample_size_Psi(N0)<=0.8 & N0<200){
		N0=N0+2;
	}
	print(N0);
	
	#N1=as.integer(sigma_nomial*(qnorm(1-alpha_local/2)+qnorm(0.5+beta*pilot_power/2))^2*1/U^2 )+40;
	
	N0=N;
	N1=as.integer(N0*pilot_size);
	N11=as.integer(N1*pilot_dtr);			N12=N1-N11;
	
	N2=N0;	N21=as.integer(N2*pivotal_dtr);	N22=N2-N21;
	
	TPNF<-matrix(rep(0,4), nrow=2, ncol=2, dimnames=list(c("True","False"), c("Positive","Negative")));
	Acc<-matrix(rep(0,6), nrow=3, ncol=2, dimnames=list(c("Pilot_True","Pilot_False", "Not_Sure"), c("Pivotal_Positive","Pivotal_Negative")));
	
	for(index in 1:iter_inside){
		pilot_indicator_1=0;	pilot_indicator_2=0;
		TF_stat=0;
		
		#stage_1
		seq1=rnorm(n=N11, mean=Period_dif+delta_nomial_Psi, sd=sqrt(2*sigma_square_true));
		seq2=rnorm(n=N12, mean=Period_dif-delta_nomial_Psi, sd=sqrt(2*sigma_square_true));
		
		delta_hat=(mean(seq1)-mean(seq2) )/2;
		sigma_s_hat=(sum((seq1-mean(seq1))^2)+sum((seq2-mean(seq2))^2) )/(N1-2);
		xi=sqrt(sigma_s_hat*(1/N11+1/N12)/4 );
		Z=(delta_hat-delta_nomial_Psi)/xi;
		
		a_U1=(U-delta_nomial_Psi)/xi - qt(p=1-alpha_Psi/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi);
		a_L1=(L-delta_nomial_Psi)/xi + qt(p=1-alpha_Psi/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi);
		r_U1=max(a_U1, (Psi*U - delta_nomial_Psi)/xi + qt(p=1-alpha_Psi/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi) );
		r_L1=min(a_L1, (Psi*L - delta_nomial_Psi)/xi - qt(p=1-alpha_Psi/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi) );
		
		if(Z<a_U1 & Z>a_L1){
			pilot_indicator_1=1;
			if(TF_stat==0){
				if(abs(delta_nomial_Psi)<U){
					TPNF["True", "Positive"]=1+TPNF["True", "Positive"];
				}
				else{
					TPNF["False", "Positive"]=1+TPNF["False", "Positive"];
				}
				TF_stat=1;
			}
		}
		
		if( ((r_L1< Z & Z<=a_L1) | (a_U1<= Z & Z<r_U1) )){
			pilot_indicator_1=3;
		}
		
		if( (Z<=r_L1 | Z>=r_U1) ){
			pilot_indicator_1=2;
			if(TF_stat==0){
				if(abs(delta_nomial_Psi)<U){
					TPNF["True", "Negative"]=1+TPNF["True", "Negative"];
				}
				else{
					TPNF["False", "Negative"]=1+TPNF["False", "Negative"];
				}
				TF_stat=1;
			}
		}
		#	pilot_indicator means that
		#	0:?		1:pass		2:fail
		#stage_2
		seq1=c(seq1, rnorm(n=N21, mean=Period_dif+delta_nomial, sd=sqrt(2*sigma_square_true)) );
		seq2=c(seq2, rnorm(n=N22, mean=Period_dif-delta_nomial, sd=sqrt(2*sigma_square_true)) );
		
		delta_hat=(mean(seq1)-mean(seq2))/2;
		sigma_s_hat=(sum((seq1-mean(seq1))^2)+sum((seq2-mean(seq2))^2) )/(N1+N2-2) ;
		xi=sqrt(sigma_s_hat*(1/(N11+N21)+1/(N12+N22) )/4);
		Z=(delta_hat-delta_nomial_Psi)/xi;
		
		a_U2=(U-delta_nomial_Psi)/xi - qt(p=1-alpha_Psi/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi);
		a_L2=(L-delta_nomial_Psi)/xi + qt(p=1-alpha_Psi/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi);
		if(Z<a_U2 & Z>a_L2){
			if(TF_stat==0){
				if(abs(delta_nomial_Psi)<U){
					TPNF["True", "Positive"]=1+TPNF["True", "Positive"];
				}
				else{
					TPNF["False", "Positive"]=1+TPNF["False", "Positive"];
				}
				TF_stat=1;
			}
			pilot_indicator_2=1;
		}
		
		else{
			if(TF_stat==0){
				if(abs(delta_nomial_Psi)<U){
					TPNF["True", "Negative"]=1+TPNF["True", "Negative"];
				}
				else{
					TPNF["False", "Negative"]=1+TPNF["False", "Negative"];
				}
				TF_stat=1;
			}
			pilot_indicator_2=2;
		}
		
		Acc[pilot_indicator_1, pilot_indicator_2]=1+Acc[pilot_indicator_1, pilot_indicator_2];	
	}
	print(Acc);
	return(TPNF);
}

#	remember to add limitations 
Psi_Simulation(N=70)
delta_nomial=U+0.01
delta_nomial=0.01



beta=0.6
delta_true=log(GMR)
change_CV=0.32431
CV_W_Inverse(0.32431)
sigma_square_true
for(prepare_again in c(1)){
	sigma_square_true=CV_W_Inverse(change_CV);
	Empherical_power_size_loc<-function(x){
		return(Empherical_power_core(lambda_loc = 0, sample_size = x)$value);
	}
	n0=P_uniroot(f=Empherical_power_size_loc);
	
	n=n0;	n=as.integer(n/2.0)*2+2;
	n=max(n, 16);
	n1=as.integer(n/2.0);
	n2=n-n1;
	
	n_star=as.integer(max(n*pilot_size/2.0, 5.0))*2+2;
	n1_star=as.integer(n_star/2.0);
	n2_star=n_star-n1_star;
};rm(prepare_again);
n
n_star












death_list=as.character(ls() );death_list_1=death_list[1]
for(name in death_list[2:length(death_list)]){
	death_list_1<-paste(death_list_1, name, sep=",")
}
death_list_1;rm(death_list, name, death_list_1)


