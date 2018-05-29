# TODO: Add comment
# 
# Author: Ted under credit of Young Lee
###############################################################################


Psi_Simulation_Const<-function(iter=5000, Psi=0.5, Max=200, N=12, beta1=0.5, beta2=0.8, alpha_AL=0.05){
	
	iter_inside=as.integer(iter);
	sigma_nomial_Psi=2*sigma_nomial;
	delta_nomial_Psi=delta_nomial;
	N0=12;
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
	
	while(sample_size_Psi(N0, Psi_l = Psi)<=0.8 & N0<100){N0=N0+2;}
	print(N0);
	N0=max(N0, N);
	print(N0);
	
	#N1=as.integer(sigma_nomial*(qnorm(1-alpha_local/2)+qnorm(0.5+beta*pilot_power/2))^2*1/U^2 )+40;
	
	N1=as.integer(N0*pilot_size);
	N11=as.integer(N1*pilot_dtr);			N12=N1-N11;
	
	N2=N0;	N21=as.integer(N2*pivotal_dtr);	N22=N2-N21;
	
	TPNF<-matrix(rep(0,4), nrow=2, ncol=2, dimnames=list(c("True","False"), c("Positive","Negative")));
	Acc<-matrix(rep(0,6), nrow=3, ncol=2, dimnames=list(c("Pilot_True","Pilot_False", "Not_Sure"), c("Pivotal_Positive","Pivotal_Negative")));
	Brown_Acc<-matrix(rep(0,6), nrow=3, ncol=2, dimnames=list(c("B_Pilot_True","B_Pilot_False", "B_Not_Sure"), c("B_Pivotal_Positive","B_Pivotal_Negative")));
	
	for(index in 1:iter_inside){
		pilot_indicator_1=0;	pilot_indicator_2=0;	pilot_indicator_3=0;
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
		
		t=N1/(N1+N0);
		CP=pnorm(delta_hat*sqrt(4/( (1/N11+1/N12)*(t-t*t)*sigma_nomial_Psi ) ) - qnorm(1-alpha_AL/2)/sqrt(1-t) );
		
		if(Z<a_U1 & Z>a_L1){
			pilot_indicator_1=1;
			if(CP<beta2){
				pilot_indicator_3=1;}
			if(CP>=beta2){
				pilot_indicator_3=3;}
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
			if(CP>=beta1){
				pilot_indicator_3=2;}
			if(CP<1-beta2){
				pilot_indicator_3=1;}
			if(pilot_indicator_3==0){
				pilot_indicator_3=3;}
		}
		
		if( (Z<=r_L1 | Z>=r_U1) ){
			pilot_indicator_1=2;
			if(pilot_indicator_3==0 & CP>beta1){
				pilot_indicator_3=2;}
			if(pilot_indicator_3==0 & CP<=beta1){
				pilot_indicator_3=3;}
			
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
		Brown_Acc[pilot_indicator_3, pilot_indicator_2]=1+Acc[pilot_indicator_3, pilot_indicator_2];	
	}
	print(Acc);
	print(Brown_Acc);
	return(TPNF);
}
Psi_Simulation_Const(Psi=0.9, N=12, beta1=0.5, beta2=0.8)


