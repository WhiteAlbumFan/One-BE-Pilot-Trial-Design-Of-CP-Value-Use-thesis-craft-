# TODO: Add comment
# 
# Author: Ted under credit of Young Lee
############################################################################################################
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
	CV_W<-function(s_square_local=0){sqrt(exp(s_square_local)-1);};
	CV_W_Inverse<-function(cv_local=0.5){log(cv_local^2+1);};
	# coefficient of variable
};rm(easy_func);


############################################################################################################
test_tears<-function(n=12){
	U*U/(4*2*sigma_nomial*qt(1-alpha/2, df=n-2))*(n);
};

for(time in 12:100){
	print(test_tears(time) );
};rm(time);
# this test_tears function is meant to show the values of one specification during BE design
CV_W(0.05)





Psi_Simulation<-function(iter=5000, Psi_initial="Assign", Psi_assign=0.5, Max=100, N=NA){
	
	iter_inside=as.integer(iter);		N0=12;
	sigma_nomial_Psi=2*sigma_nomial;	sigma_s_true_Psi=2*sigma_square_true;
	sample_size_Psi<-function(n=60){
		n1=as.integer(n*pivotal_dtr);		n2=n-n1;
		n_star=max(as.integer(n*pilot_size), 12);
		n1_star=as.integer(n_star*pilot_dtr);	n2_star=n_star-n1_star;
		
		
		time1=4/(1/n1_star+1/n2_star);		time2=4/( 1/(n1+n1_star)+1/(n2+n2_star) );
		time3=4/(1/n1+1/n2);				C1=2*sigma_nomial_Psi*qt(1-alpha/2, df=n_star-2)^2*1/(U*U);
		#C1=1;
		#C1=U*U/(4*sigma_nomial_Psi*qt(1-alpha/2, df=n_star-2))
		
		if(Psi_initial=="Sqrt_Const"){Psi_l=1-sqrt(C1/time1);};
		if(Psi_initial=="Sqrt"){Psi_l=1-sqrt(1/time1);};
		if(Psi_initial=="Assign"){Psi_l=Psi_assign;};
		
		return(
				integrate(function(x){
							#print(x);
							G1=qt(p=1-alpha/2, df=n_star-2)/sqrt(n_star-2);
							a_U1=(U-delta_nomial)/sqrt( sigma_nomial_Psi/time1) -G1*sqrt(x);					
							a_L1=(L-delta_nomial)/sqrt( sigma_nomial_Psi/time1) +G1*sqrt(x);
							#print(a_L1);
							r_U1=max(a_U1, (Psi_l*U-delta_nomial)/sqrt( sigma_nomial_Psi/time1 )+G1*sqrt(x) );
							#print(r_U1);
							r_L1=min(a_L1, (Psi_l*L-delta_nomial)/sqrt( sigma_nomial_Psi/time1 )-G1*sqrt(x) );
							XX=pnorm(a_U1)-pnorm(a_L1)+
									integrate(function(x1){
												G2=qt(p=1-alpha/2, df=n+n_star-2);
												a_U=(U-delta_nomial)/sqrt( sigma_nomial_Psi/time2 )-G2;
												a_L=(L-delta_nomial)/sqrt( sigma_nomial_Psi/time2 )+G2;
												return(dnorm(x1)*max(pnorm(a_U*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1)- 
																		pnorm(a_L*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1), 0) );
											},lower=r_L1, upper=a_L1 )$value+
									integrate(function(x1){
												G2=qt(p=1-alpha/2, df=n+n_star-2);
												a_U=(U-delta_nomial)/sqrt( sigma_nomial_Psi/time2 )-G2;
												a_L=(L-delta_nomial)/sqrt( sigma_nomial_Psi/time2 )+G2;
												return(dnorm(x1)*max(pnorm(a_U*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1)-
																		pnorm(a_L*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1), 0)  );
											},lower=a_U1, upper=r_U1 )$value;
							return(XX*dchisq(x, df=n_star-2));
						},lower=(n_star-2)/2, upper=2*(n_star-2) )$value);
	}
	
	if(is.na(N)){
		K=sample_size_Psi(n=N0);if(is.na(K)){K=0;}
		while(K<=beta & N0<Max ){
			N0=N0+2;K=sample_size_Psi(n=N0);
			if(is.na(K)){K=0;}
		};rm(K);
	}
	else{N0=N;};
	print("SampleSizeE.S.T.");print(N0);  print("Traditional:");
	print( as.integer(2*sigma_nomial_Psi*(qnorm(1-alpha/2)+qnorm(0.5+beta*pilot_power/2))^2*1/U^2 ));
	
	for(prepare_Psi in c(1) ){
		N1=max(as.integer(N0*pilot_size), 6);
		N11=as.integer(N1*pilot_dtr);	N12=N1-N11;
		
		N2=N0;N21=as.integer(N2*pivotal_dtr);	N22=N2-N21;
		time1=4/(1/N11+1/N12);			C1=2*sigma_nomial_Psi*qt(1-alpha/2, df=N1-2)^2*1/(U*U);
		
		if(Psi_initial=="Sqrt_Const"){Psi=1-sqrt(C1/time1);};
		if(Psi_initial=="Sqrt"){Psi=1-sqrt(1/time1);};
		if(Psi_initial=="Assign"){Psi=Psi_assign;};	rm(time1, C1);
		
		TPNF<-matrix(rep(0,4), nrow=2, ncol=2, dimnames=list(c("True","False"),
		                                                     c("Positive","Negative")));
		Acc<-matrix(rep(0,6), nrow=3, ncol=2, dimnames=list(c("Pilot_True","Pilot_False", "Not_Sure"),
		                                                    c("Pivotal_Positive","Pivotal_Negative")));
		Fact="True";if(abs(delta_true)>=U){Fact="False";};
	}
	
	for(index in 1:iter_inside){
		TF_stat=0;				pivotal_indicator=0;
		pilot_indicator_1=0;
		
		#stage_1
		seq1=rnorm(n=N11, mean=Period_dif+delta_true, sd=sqrt(sigma_s_true_Psi));
		seq2=rnorm(n=N12, mean=Period_dif-delta_true, sd=sqrt(sigma_s_true_Psi));
		
		delta_hat=(mean(seq1)-mean(seq2) )/2;
		sigma_s_hat=(sum((seq1-mean(seq1))^2)+sum((seq2-mean(seq2))^2) )/(N1-2);
		xi=sqrt(sigma_nomial_Psi*(1/N11+1/N12)/4 );
		Z=(delta_hat-delta_nomial)/xi;
		
		a_U1=(U-delta_nomial)/xi - qt(p=1-alpha/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi);
		a_L1=(L-delta_nomial)/xi + qt(p=1-alpha/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi);
		r_U1=max(a_U1, (Psi*U - delta_nomial)/xi + qt(p=1-alpha/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi) );
		r_L1=min(a_L1, (Psi*L - delta_nomial)/xi - qt(p=1-alpha/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi) );
		

		for(CP_Late in c(1)){
			if(Z<a_U1 & Z>a_L1){
				pilot_indicator_1=1;
			};
			if( ((r_L1< Z & Z<=a_L1) | (a_U1<=Z & Z<r_U1) )){
				pilot_indicator_1=3;
			};
			if( (Z<=r_L1 | Z>=r_U1) ){
				pilot_indicator_1=2;
			};
		};rm(CP_Late);
		

		
		#stage_2
		seq1=c(seq1, rnorm(n=N21, mean=Period_dif+delta_true, sd=sqrt(sigma_s_true_Psi)) );
		seq2=c(seq2, rnorm(n=N22, mean=Period_dif-delta_true, sd=sqrt(sigma_s_true_Psi)) );
		
		delta_hat=(mean(seq1)-mean(seq2))/2;
		sigma_s_hat=(sum((seq1-mean(seq1))^2)+sum((seq2-mean(seq2))^2) )/(N1+N2-2) ;
		xi=sqrt(sigma_nomial_Psi*(1/(N11+N21)+1/(N12+N22) )/4);
		Z=(delta_hat-delta_nomial)/xi;
		
		a_U2=(U-delta_nomial)/xi - qt(p=1-alpha/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi);
		a_L2=(L-delta_nomial)/xi + qt(p=1-alpha/2, df=N1-2)*sqrt(sigma_s_hat/sigma_nomial_Psi);
		
		if(Z<a_U2 & Z>a_L2){pivotal_indicator=1;}
		else{pivotal_indicator=2;};
		
		Acc[pilot_indicator_1, pivotal_indicator]=1+Acc[pilot_indicator_1, pivotal_indicator];	
	}
	
	for(results in c(1)){
		print(Acc);	
		TPNF[Fact, "Positive"]=Acc[1,1]+Acc[1,2]+Acc[3,1];
		TPNF[Fact, "Negative"]=iter_inside-TPNF[Fact, "Positive"];
		print(TPNF);
	};rm(results);
	
	return(TPNF);
}


#	remember to add limitations 
Psi_Simulation(iter=5000, Psi_initial = "Assign", Psi_assign=0.5)
beta=0.5

pilot_size=0.25

60*pilot_size*(40+181+477)/5000+60*(4300/5000)	#(4+477+4044)/5000;	(477+36)/(5000-258-4044);	1-(258+4044)/5000
60*4224/5000+70*pilot_size*(1-4224/5000)		#(519+42+3668)/5000;(519+16)/598	;561/5000	(561+37)/5000

delta_true=1*U
delta_true=0.01

sigma_square_true=0.1;sigma_nomial=gamma*sigma_square_true
#U^2/(4*sigma_nomial*qnorm(1-alpha/2)^2)







Psi_Simulation_Sqrt_B<-function(iter=5000, Max=200, N=12){
	
	iter_inside=as.integer(iter);
	
	#necessary yakso
	#here d_ij var=sigma^2
	alpha_Psi=0.10;
	sigma_nomial_Psi=2*sigma_nomial;
	delta_nomial_Psi=delta_nomial;
	N0=12;
	sample_size_Psi<-function(n){
		
		Psi_l=1-1/sqrt(n);
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
							print(Psi_l);
							XX=pnorm(a_U1)-pnorm(a_L1)+
									integrate(function(x1){
												G2=qt(p=1-alpha/2, df=n+n_star-2);
												a_U=(U-delta_nomial)/sqrt( sigma_nomial/time2 )-G2;
												a_L=(L-delta_nomial)/sqrt( sigma_nomial/time2 )+G2;
												return(dnorm(x1)*max(pnorm(a_U*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1) - 
																		pnorm(a_L*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1), 0) );
											},
											lower=r_L1, upper=a_L1 )$value+
									integrate(function(x1){
												G2=qt(p=1-alpha/2, df=n+n_star-2);
												a_U=(U-delta_nomial)/sqrt( sigma_nomial/time2 )-G2;
												a_L=(L-delta_nomial)/sqrt( sigma_nomial/time2 )+G2;
												return(dnorm(x1)*max(pnorm(a_U*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1) - 
																		pnorm(a_L*(n+n_star)/n*sqrt(time3/time2)-sqrt(n_star^2*time2/(n^2*time1) )*x1), 0)  );
											},
											lower=a_U1, upper=r_U1 )$value;
							return(XX*dchisq(x, df=n_star-2));},
						lower=(n_star-2)^2*(1-Psi_l)*0.0601, upper=(n_star-2)^2*(1)*0.0601)$value );
	}
	
	while(sample_size_Psi(N0)<=0.8 & N0<200){
		N0=N0+2;
	}
	print(N0);
	N0=max(N0, N);
	print(N0);
	
	#N1=as.integer(sigma_nomial*(qnorm(1-alpha_local/2)+qnorm(0.5+beta*pilot_power/2))^2*1/U^2 )+40;
	Psi=1-1/sqrt(N0);
	N1=as.integer(N0*pilot_size);
	N11=as.integer(N1*pilot_dtr);			N12=N1-N11;
	
	N2=N0;	N21=as.integer(N2*pivotal_dtr);	N22=N2-N21;
	
	TPNF<-matrix(rep(0,4), nrow=2, ncol=2, dimnames=list(c("True","False"), c("Positive","Negative")));
	Acc<-matrix(rep(0,6), nrow=3, ncol=2, dimnames=list(c("Pilot_True","Pilot_False", "Not_Sure"), c("Pivotal_Positive","Pivotal_Negative")));
	
	for(index in 1:iter_inside){
		pilot_indicator_1=0;	pivotal_indicator=0;
		TF_stat=0;
		
		#stage_1
		seq1=rnorm(n=N11, mean=Period_dif+delta_true, sd=sqrt(2*sigma_square_true));
		seq2=rnorm(n=N12, mean=Period_dif-delta_true, sd=sqrt(2*sigma_square_true));
		
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
		seq1=c(seq1, rnorm(n=N21, mean=Period_dif+delta_true, sd=sqrt(2*sigma_square_true)) );
		seq2=c(seq2, rnorm(n=N22, mean=Period_dif-delta_true, sd=sqrt(2*sigma_square_true)) );
		
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
			pivotal_indicator=1;
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
			pivotal_indicator=2;
		}
		
		Acc[pilot_indicator_1, pivotal_indicator]=1+Acc[pilot_indicator_1, pivotal_indicator];	
	}
	print(Acc);
	return(TPNF);
}


















clear_list=as.character(ls() );clear_list_1=clear_list[1]
for(name in clear_list[2:length(clear_list)]){
  clear_list_1<-paste(clear_list_1, name, sep=",")
}
clear_list_1;rm(clear_list, name, clear_list_1)







