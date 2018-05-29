# TODO: Add comment
# 
# Author: Ted under credit of Young Lee
###############################################################################


for(prepare_const in 1:1){
	
	GMR=0.95;	delta_true=log(GMR);
	delta_nomial=delta_true+0.01*rnorm(1);
	#	original assumption
	
	Period_dif=7e-1;				#	futility-para_1
	
	sigma_square_true=0.25;			#	major para
	
	alpha=0.10;						#	control TIE
	beta=0.80;						#	power/recall rate
	
	theta=0.80;						#	tolerance para
	U=-log(theta);
	L=log(theta);
	
	gamma=1.1;						#	rate of est sigma_s to real sigma_s
	
	size_rate=0.5;					#	pilot-magnitude
	lambda=0;
	pilot_press=0.25;	
	
	#n=100;		n1=45;		n2=n-n1;
	#	for pivotal size/re-est sample size
	
	#n_star=40;	n1_star=25;	n2_star=n_star-n1_star;
	#	pilot size
	
};rm(prepare_const);



CV_W<-function(s_square_local=0){
	sqrt(exp(s_square_local)-1);}
CV_W_Inverse<-function(cv_local=0.5){
	log(cv_local^2+1);
}

test_stat<-function(n1_local=10, n_local=50, delta_local=0.0, 
		Period_local=1e-1, s_square_local=5e-1, oringinal="Yes")
{
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



Empherical_power_core<-function(lambda_loc=lambda, delta_nomials=delta_nomial, 
		sigma_nomial=gamma*sigma_square_true, sample_size=12, sample_dtr=0.5){
	#	coef_loc = sqrt(n_11^{-1}+n_12^{-1});
	n_local=max(as.integer(sample_size), 12);
	n1_local=max((n_local*sample_dtr), 4);
	t_df=n_local-2;
	
	coef_loc=sqrt(1/n1_local+1/(n_local-n1_local));
	return(integrate(function(x){
						a=pnorm((U-lambda_loc-delta_nomials)/sqrt(sigma_nomial*coef_loc^2*1/2)-qt(p=1-alpha/2, df=t_df)*sqrt(2*x*x/sigma_nomial) );
						b=pnorm((L+lambda_loc-delta_nomials)/sqrt(sigma_nomial*coef_loc^2*1/2)-qt(p=1-alpha/2, df=t_df)*sqrt(2*x*x/sigma_nomial) );
						return(a-b)*(exp(-x/2)*x^(t_df/2-1)*1/(2^(t_df/2)*gamma(t_df/2)));
					}, lower=0, upper=Inf) );
}

P_uniroot<-function(f=function(x){return(x);}, beta_loc=beta){
	jump_s=1;
	temp_loc=f(300);
	if(temp_loc<=beta_loc){
		print("not available;it exceeds Sup");
		return(300);
	}
	N_L=12.0;	N_U=300.0;
	counter=0;
	while(jump_s>0 & counter<20){
		temp=(N_U+N_L)/2;
		temp_loc=f(temp);
		if(temp_loc>beta_loc){
			N_U=temp;
		}
		else{
			N_L=temp;
		}
		if(N_U-N_L<=1){
			jump_s=-1;
			counter=22;
		}
		counter=counter+1;
	}
	return(as.integer(N_U));
}




size_rate=0.3
delta_nomial
delta_true

Wang_Simulation(method="Common_use", iter=400000, operation = "strict")
Wang_Simulation(method="Wang", iter=400000, operation = "strict")
Wang_Simulation(method="Wang_pool", iter=400000, operation = "strict")
Wang_Simulation(method="Wang_old", iter=400000, operation = "strict")

Wang_Simulation(method="Common_use", iter=400000, operation = "not_strict")
Wang_Simulation(method="Wang", iter=400000, operation = "not_strict")
Wang_Simulation(method="Wang_pool", iter=400000, operation = "not_strict")
Wang_Simulation(method="Wang_old", iter=400000, operation = "not_strict")

#	Common_use		
beta=0.8
delta_true=log(GMR)
change_CV=0.32431
CV_W_Inverse(0.32431)
sigma_square_true
for(prepare_again in c(1)){
	#sigma_square_true=CV_W_Inverse(change_CV);
	Empherical_power_size_loc<-function(x){
		return(Empherical_power_core(lambda_loc = 0, sample_size = x)$value);
	}
	n0=P_uniroot(f=Empherical_power_size_loc, beta_loc = beta);
	
	n=n0;	n=as.integer(n/2.0)*2+2;
	n=max(n, 16);
	n1=as.integer(n/2.0);
	n2=n-n1;
	
	n_star=as.integer(max(n*size_rate/2.0, 5.0))*2+2;
	n1_star=as.integer(n_star/2.0);
	n2_star=n_star-n1_star;
};rm(prepare_again);
n
n_star
#	0.32431	0.75406	0.7555525	0.910695	0.7550125
#	0.25	0.22606	0.753435	0.61554		0.48647		
delta_nomial

?setwd()


#x=c("delta_true", "delta_nomial", "gamma", "sigma_square_true", "size_rate", "pilot_press", "beta", "Real_Percent", "Fake_Percent");

x=c("Method", "state", "delta_true", "delta_nomial", "gamma", "sigma_square_true", "size_rate", "pilot_press", "beta", "Real_Percent", "Fake_Percent", "N", "N_star");

write(x, append=TRUE, file="D:/JAVA_Palace/learning_pilot_trial_simulation/Verify_onWang_20180308/news.txt", sep=",")



for(prepare_4 in 1:5){

	for(prepare_tri in 2:10){
	
		#Wang_Simulation(method="Wang", iter=400000, operation = "strict")
		#Wang_Simulation(method="Wang_pool", iter=400000, operation = "strict")
		#Wang_Simulation(method="Wang_old", iter=400000, operation = "strict")
		beta=0.1*prepare_tri;
		delta_true=log(GMR);
		delta_nomial=delta_true+0.02*rnorm(1)
	
		for(prepare_again in c(1)){
			#sigma_square_true=CV_W_Inverse(change_CV);
			Empherical_power_size_loc<-function(x){
				return(Empherical_power_core(lambda_loc = 0, sample_size = x)$value);
			}
			n0=P_uniroot(f=Empherical_power_size_loc, beta_loc = beta);
		
			n=n0;	n=as.integer(n/2.0)*2+2;
			n=max(n, 16);
			n1=as.integer(n/2.0);
			n2=n-n1;
		
			n_star=as.integer(max(n*size_rate/2.0, 5.0))*2+2;
			n1_star=as.integer(n_star/2.0);
			n2_star=n_star-n1_star;
		};rm(prepare_again);
	
		delta_true=log(GMR);
		Real_Percent=Wang_Simulation(method="Common_use", iter=100000, operation = "strict");
		delta_true=U;
		Fake_Percent=Wang_Simulation(method="Common_use", iter=100000, operation = "strict");
		x=c("Common_use", "strict", delta_true, delta_nomial, gamma, sigma_square_true, size_rate, pilot_press, beta, Real_Percent, Fake_Percent, n, n_star);
		write(x, append=TRUE, file="D:/JAVA_Palace/learning_pilot_trial_simulation/Verify_onWang_20180308/news.txt", sep="\t")
	
	
		delta_true=log(GMR);
		Real_Percent=Wang_Simulation(method="Wang_old", iter=100000, operation = "strict");
		delta_true=U;
		Fake_Percent=Wang_Simulation(method="Wang_old", iter=100000, operation = "strict");
		x=c("Wang_old", "strict", delta_true, delta_nomial, gamma, sigma_square_true, size_rate, pilot_press, beta, Real_Percent, Fake_Percent, n, n_star);
		write(x, append=TRUE, file="D:/JAVA_Palace/learning_pilot_trial_simulation/Verify_onWang_20180308/news.txt", sep="\t")
	
		delta_true=log(GMR);
		Real_Percent=Wang_Simulation(method="Wang_pool", iter=100000, operation = "strict");
		delta_true=U;
		Fake_Percent=Wang_Simulation(method="Wang_pool", iter=100000, operation = "strict");
		x=c("Wang_pool", "strict", delta_true, delta_nomial, gamma, sigma_square_true, size_rate, pilot_press, beta, Real_Percent, Fake_Percent, n, n_star);
		write(x, append=TRUE, file="D:/JAVA_Palace/learning_pilot_trial_simulation/Verify_onWang_20180308/news.txt", sep="\t")
	
		delta_true=log(GMR);
		Real_Percent=Wang_Simulation(method="Wang", iter=100000, operation = "strict");
		delta_true=U;
		Fake_Percent=Wang_Simulation(method="Wang", iter=100000, operation = "strict");
		x=c("Wang", "strict", delta_true, delta_nomial, gamma, sigma_square_true, size_rate, pilot_press, beta, Real_Percent, Fake_Percent, n, n_star);
		write(x, append=TRUE, file="D:/JAVA_Palace/learning_pilot_trial_simulation/Verify_onWang_20180308/news.txt", sep="\t")
	
	#x=c("Wang", "state", "delta_true", "delta_nomial", "gamma", "sigma_square_true", "size_rate", "pilot_press", "beta", "Real_Percent", "Fake_Percent");

	#x=c("Wang", "strict", delta_true, delta_nomial, gamma, sigma_square_true, size_rate, pilot_press, beta, Real_Percent, Fake_Percent);
	#x=100000;
	#write(x, append=TRUE, file="D:/JAVA_Palace/learning_pilot_trial_simulation/Verify_onWang_20180308/news.txt", sep="\t")
	}
}
































n0=12
for(indexs in 1:10){
	record=n0;	n0=as.integer(2*gamma*sigma_square_true*(qt(p=1-alpha/2, df=n0-2)+qt(p=0.5+beta/2, df=n0-2))^2*1/U^2);
	if(abs(record-n0)<=1){break;}
};if(indexs==10){print("potential under estimmate");}
rm(indexs);






