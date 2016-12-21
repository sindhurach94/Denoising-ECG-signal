clc;
x=0.01:0.01:10;
 li=30/72;  

  a_pwav=0.25;
      d_pwav=0.09;
      t_pwav=0.16;  
     
      a_qwav=0.025;
      d_qwav=0.066;
      t_qwav=0.166;
      
      a_qrswav=1.6;
      d_qrswav=0.11;
      
      a_swav=0.25;
      d_swav=0.066;
      t_swav=0.09;
      
      a_twav=0.35;
      d_twav=0.142;
      t_twav=0.2;
      
      a_uwav=0.035;
      d_uwav=0.0476;
      t_uwav=0.433;
      
       pwav=p_wav(x,a_pwav,d_pwav,t_pwav,li);
       qwav=q_wav(x,a_qwav,d_qwav,t_qwav,li);
       qrswav=qrs_wav(x,a_qrswav,d_qrswav,li);
       swav=s_wav(x,a_swav,d_swav,t_swav,li);
       twav=t_wav(x,a_twav,d_twav,t_twav,li);
       uwav=u_wav(x,a_uwav,d_uwav,t_uwav,li);
       
 ecg=pwav+qrswav+twav+swav+qwav+uwav;
 namp=[0.6184 0.3475 0.1956 0.1100 0.0618 0.0348 0.0196 0.0110];
 for k1=1:8,
 noise1=namp(k1)*cos(2*pi*10*x);
% noise = wgn(100,1,0);
input1=noise1+ecg;

snr1(k1)= mean( ecg .^ 2 ) / mean( noise1 .^ 2 );
snr_db1(k1) = 10 * log10( snr1(k1) );

desired=ecg;
filterOrder=5;
nCoeff=filterOrder+1;
nIterations=length(desired);
errorVector1=zeros(nIterations,1);
outputVector=zeros(nIterations,1);
coeffVector=zeros(nCoeff,(nIterations+1));
S_d=eye(nCoeff);
p_d=zeros(nCoeff,1);
xk=zeros(nCoeff-1,1);
lambda=0.91398;
final=zeros(10,1);

% dwt sureshrink algorithm

[c,d]=dwt(input1,'haar');
[b,a] =butter(4, 0.5, 'high') ;
dii=filter(b,a,[c,d]);
di1=downsample(dii,1,0);
 [b,a] =butter(4, 0.5) ;
 aii=filter(b,a,[c,d]);
 ai1=downsample(aii,1,0);
 
 [b,a] =butter(4, 0.5, 'high') ;
dii=filter(b,a,ai1);
di2=downsample(dii,1,0);
 [b,a] =butter(4, 0.5) ;
 aii=filter(b,a,ai1);
 ai2=downsample(aii,1,0);
 
 p=2*log(5);
 lambda=sqrt(p);
 
 for it=1:500,
     if di1(it)>lambda
         di1(it)=di1(it)-lambda;
     elseif di1(it)<lambda
         di1(it)=di1(it)+lambda;
     else
         di1(it)=0;
     end
 end
 
 
 for it=1:250
       if di2(it)>lambda
         di2(it)=di2(it)-lambda;
     elseif di2(it)<lambda
         di2(it)=di2(it)+lambda;
     else
         di2(it)=0;
     end
 end
 dij=upsample(di2,1,0);
 [bn,an] =butter(4, 0.5, 'high') ;
 dii2=filter(bn,an,dij);
 
 aij=upsample(ai2,1,0);
 [bn,an] =butter(4, 0.5) ;
 aii2=filter(bn,an,aij);
 
 a1i2=[aii2 dii2];
  a1ij=upsample(a1i2,1,0);
 [bn2,an2] =butter(4, 0.5) ;
 aii3=filter(bn2,an2,a1ij);
 
  d1ij2=upsample(di1,1,0);
 [bn2,an2] =butter(4, 0.5, 'high') ;
 dii3=filter(bn2,an2,d1ij2);
 output1=1.5*[aii3 dii3];
  y=idwt(aii3(1:1:1000),dii3,'haar');
 [v2 z1]=alignsignals(y,desired,20);
 errorVector1=z1-v2(1:1:1012);

[k j]=alignsignals(y,input1,20);
err1(k1)=0;
for i=1:1:1000,
    err1(k1)=err1(k1)+((j(i)-k(i)).^2);
end
err_sureshrink1(k1)=err1(k1);

prd_sureshrink1(k1)=0;
sum11=0;
sum21=0;
for i=1:1:nIterations-5,
    sum11=sum11+(k(i)-j(i)).^2;
    sum21=sum21+(j(i)).^2;
end
prd_sureshrink1(k1)=100*sqrt(sum11/sum21);

snre_sureshrink1(k1)=0;
sum11=0;
sum21=0;
for i=1:1:nIterations-5,
    sum11=sum11+(y(i).^2);
    sum21=sum21+((j(i)-k(i)).^2);
end
snre_sureshrink1(k1)=10*log10(sum11/sum21);

%dwt neighblock
lambdac1=4.505;
z1=(log10(5)/2)+(2*max(1,log10(5)/4));
sum1=0;
i=1;
j=1;
s=0;
[c,d]=dwt(input1,'haar');
while i<=495
   ci= c(i:1:i+5).^2;
    for a=1:5,
        s=s+ci(a);
    end
    while j<=i+5
        c1(j)=(max(0,(s.^2-lambdac1*z1)/s.^2))*c(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end
i=1;
j=1;
s=0;
while i<=495
    di=d(i:1:i+5).^2;
    for a=1:5,
        s=s+di(a);
    end
    while j<=i+5
        d1(j)=(max(0,(s.^2-lambdac1*z1)/s.^2))*d(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end

output1=1.095*idwt(c1,d1,'haar');

errorVector1=desired(1:1:996)-output1;


err1(k1)=0;
for i=1:1:nIterations-5,
    err1(k1)=err1(k1)+((input1(i)-output1(i)).^2);
end
err_neighblock1(k1)=err1(k1);

prd_neighblock1(k1)=0;
sum11=0;
sum21=0;
for i=1:1:nIterations-5,
    sum11=sum11+(input1(i)-output1(i)).^2;
    sum21=sum21+(input1(i)).^2;
end
prd_neighblock1(k1)=100*sqrt(sum11/sum21);
snre_neighblock1(k1)=0;
sum11=0;
sum21=0;
for i=1:1:nIterations-5,
    sum11=sum11+(output1(i).^2);
    sum21=sum21+((input1(i)-output1(i)).^2);
end
snre_neighblock1(k1)=10*log10(sum11/sum21);

 end

%%
clc;
x=0.01:0.01:10;
 li=30/72;  

  a_pwav=0.25;
      d_pwav=0.09;
      t_pwav=0.16;  
     
      a_qwav=0.025;
      d_qwav=0.066;
      t_qwav=0.166;
      
      a_qrswav=1.6;
      d_qrswav=0.11;
      
      a_swav=0.25;
      d_swav=0.066;
      t_swav=0.09;
      
      a_twav=0.35;
      d_twav=0.142;
      t_twav=0.2;
      
      a_uwav=0.035;
      d_uwav=0.0476;
      t_uwav=0.433;
      
       pwav=p_wav(x,a_pwav,d_pwav,t_pwav,li);
       qwav=q_wav(x,a_qwav,d_qwav,t_qwav,li);
       qrswav=qrs_wav(x,a_qrswav,d_qrswav,li);
       swav=s_wav(x,a_swav,d_swav,t_swav,li);
       twav=t_wav(x,a_twav,d_twav,t_twav,li);
       uwav=u_wav(x,a_uwav,d_uwav,t_uwav,li);
       
 ecg=pwav+qrswav+twav+swav+qwav+uwav;
 namp=[0.6184 0.3475 0.1956 0.1100 0.0618 0.0348 0.0196 0.0110];
 for k1=1:8,
 noise2=namp(k1)*cos(2*pi*50*x);
% noise = wgn(100,1,0);
input2=noise2+ecg;

snr2(k1)= mean( ecg .^ 2 ) / mean( noise2 .^ 2 );
snr_db2(k1) = 10 * log10( snr2(k1) );

desired=ecg;
filterOrder=5;
nCoeff=filterOrder+1;
nIterations=length(desired);
errorVector2=zeros(nIterations,1);
outputVector=zeros(nIterations,1);
coeffVector=zeros(nCoeff,(nIterations+1));
S_d=eye(nCoeff);
p_d=zeros(nCoeff,1);
xk=zeros(nCoeff-1,1);
lambda2=0.91398;
final=zeros(10,1);

% dwt sureshrink algorithm

[c2,d2]=dwt(input2,'haar');
[b2,a2] =butter(4, 0.5, 'high') ;
dii2=filter(b2,a2,[c2,d2]);
di12=downsample(dii2,1,0);
 [b2,a2] =butter(4, 0.5) ;
 aii2=filter(b2,a2,[c2,d2]);
 ai12=downsample(aii2,1,0);
 
 [b2,a2] =butter(4, 0.5, 'high') ;
dii2=filter(b2,a2,ai12);
di22=downsample(dii2,1,0);
 [b2,a2] =butter(4, 0.5) ;
 aii2=filter(b2,a2,ai12);
 ai22=downsample(aii2,1,0);
 
 p2=2*log(5);
 lambda2=sqrt(p2);
 
 for it=1:500,
     if di12(it)>lambda2
         di12(it)=di12(it)-lambda2;
     elseif di12(it)<lambda2
         di12(it)=di12(it)+lambda2;
     else
         di12(it)=0;
     end
 end
 
 
 for it=1:250
       if di22(it)>lambda2
         di22(it)=di22(it)-lambda2;
     elseif di22(it)<lambda2
         di22(it)=di22(it)+lambda2;
     else
         di22(it)=0;
     end
 end
 dij2=upsample(di22,1,0);
 [bn2,an2] =butter(4, 0.5, 'high') ;
 dii2=filter(bn2,an2,dij2);
 
 aij2=upsample(ai22,1,0);
 [bn2,an2] =butter(4, 0.5) ;
 aii2=filter(bn2,an2,aij2);
 
 a1i2=[aii2 dii2];
  a1ij2=upsample(a1i2,1,0);
 [bn2,an2] =butter(4, 0.5) ;
 aii32=filter(bn2,an2,a1ij2);
 
  d1ij2=upsample(di12,1,0);
 [bn2,an2] =butter(4, 0.5, 'high') ;
 dii32=filter(bn2,an2,d1ij2);
 output2=1.5*[aii32 dii32];
  y2=idwt(aii32(1:1:1000),dii32,'haar');
 [v2 z2]=alignsignals(y2,desired,20);
 errorVector2=z2-v2(1:1:1012);

[k j]=alignsignals(y2,input2,20);
err2(k1)=0;
for i=1:1:1000,
    err2(k1)=err2(k1)+((j(i)-k(i)).^2);
end
err_sureshrink2(k1)=err2(k1);

prd_sureshrink2(k1)=0;
sum12=0;
sum22=0;
for i=1:1:nIterations-5,
    sum12=sum12+(k(i)-j(i)).^2;
    sum22=sum22+(j(i)).^2;
end
prd_sureshrink2(k1)=100*sqrt(sum12/sum22);

snre_sureshrink2(k1)=0;
sum12=0;
sum22=0;
for i=1:1:nIterations-5,
    sum12=sum12+(y2(i).^2);
    sum22=sum22+((j(i)-k(i)).^2);
end
snre_sureshrink2(k1)=10*log10(sum12/sum22);

%dwt neighblock
lambdac2=4.505;
z2=(log10(5)/2)+(2*max(1,log10(5)/4));
sum2=0;
i=1;
j=1;
s=0;
[c2,d2]=dwt(input2,'haar');
while i<=495
   ci= c2(i:1:i+5).^2;
    for a2=1:5,
        s=s+ci(a2);
    end
    while j<=i+5
        c2(j)=(max(0,(s.^2-lambdac2*z2)/s.^2))*c2(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end
i=1;
j=1;
s=0;
while i<=495
    di=d2(i:1:i+5).^2;
    for a2=1:5,
        s=s+di(a2);
    end
    while j<=i+5
        d1(j)=(max(0,(s.^2-lambdac2*z2)/s.^2))*d2(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end

output2=1.095*idwt(c2,d2,'haar');

% errorVector2=desired(1:1:996)-output2;


err2(k1)=0;
for i=1:1:nIterations-5,
    err2(k1)=err2(k1)+((input2(i)-output2(i)).^2);
end
err_neighblock2(k1)=err2(k1);

prd_neighblock2(k1)=0;
sum12=0;
sum22=0;
for i=1:1:nIterations-5,
    sum12=sum12+(input2(i)-output2(i)).^2;
    sum22=sum22+(input2(i)).^2;
end
prd_neighblock2(k1)=100*sqrt(sum12/sum22);
snre_neighblock2(k1)=0;
sum12=0;
sum22=0;
for i=1:1:nIterations-5,
    sum12=sum12+(output2(i).^2);
    sum22=sum22+((input2(i)-output2(i)).^2);
end
snre_neighblock2(k1)=10*log10(sum12/sum22);

 end
%%
clc;
x=0.01:0.01:10;
 li=30/72;  

  a_pwav=0.25;
      d_pwav=0.09;
      t_pwav=0.16;  
     
      a_qwav=0.025;
      d_qwav=0.066;
      t_qwav=0.166;
      
      a_qrswav=1.6;
      d_qrswav=0.11;
      
      a_swav=0.25;
      d_swav=0.066;
      t_swav=0.09;
      
      a_twav=0.35;
      d_twav=0.142;
      t_twav=0.2;
      
      a_uwav=0.035;
      d_uwav=0.0476;
      t_uwav=0.433;
      
       pwav=p_wav(x,a_pwav,d_pwav,t_pwav,li);
       qwav=q_wav(x,a_qwav,d_qwav,t_qwav,li);
       qrswav=qrs_wav(x,a_qrswav,d_qrswav,li);
       swav=s_wav(x,a_swav,d_swav,t_swav,li);
       twav=t_wav(x,a_twav,d_twav,t_twav,li);
       uwav=u_wav(x,a_uwav,d_uwav,t_uwav,li);
       
 ecg=pwav+qrswav+twav+swav+qwav+uwav;
 namp=[0.6184 0.3475 0.1956 0.1100 0.0618 0.0348 0.0196 0.0110];
 for k1=1:8,
 noise3=namp(k1)*cos(2*pi*80*x);
% noise = wgn(100,1,0);
input3=noise3+ecg;

snr3(k1)= mean( ecg .^ 2 ) / mean( noise3 .^ 2 );
snr_db3(k1) = 10 * log10( snr3(k1) );

desired=ecg;
filterOrder=5;
nCoeff=filterOrder+1;
nIterations=length(desired);
errorVector3=zeros(nIterations,1);
outputVector=zeros(nIterations,1);
coeffVector=zeros(nCoeff,(nIterations+1));
S_d=eye(nCoeff);
p_d=zeros(nCoeff,1);
xk=zeros(nCoeff-1,1);
lambda3=0.91398;
final=zeros(10,1);

% dwt sureshrink algorithm

[c3,d3]=dwt(input3,'haar');
[b3,a3] =butter(4, 0.5, 'high') ;
dii3=filter(b3,a3,[c3,d3]);
di13=downsample(dii3,1,0);
 [b3,a3] =butter(4, 0.5) ;
 aii3=filter(b3,a3,[c3,d3]);
 ai13=downsample(aii3,1,0);
 
 [b3,a3] =butter(4, 0.5, 'high') ;
dii3=filter(b3,a3,ai13);
di23=downsample(dii3,1,0);
 [b3,a3] =butter(4, 0.5) ;
 aii3=filter(b3,a3,ai13);
 ai23=downsample(aii3,1,0);
 
 p3=2*log(5);
 lambda3=sqrt(p3);
 
 for it=1:500,
     if di13(it)>lambda3
         di13(it)=di13(it)-lambda3;
     elseif di13(it)<lambda3
         di13(it)=di13(it)+lambda3;
     else
         di13(it)=0;
     end
 end
 
 
 for it=1:250
       if di23(it)>lambda3
         di23(it)=di23(it)-lambda3;
     elseif di23(it)<lambda3
         di23(it)=di23(it)+lambda3;
     else
         di23(it)=0;
     end
 end
 dij3=upsample(di23,1,0);
 [bn3,an3] =butter(4, 0.5, 'high') ;
 dii3=filter(bn3,an3,dij3);
 
 aij3=upsample(ai23,1,0);
 [bn3,an3] =butter(4, 0.5) ;
 aii3=filter(bn3,an3,aij3);
 
 a1i3=[aii3 dii3];
  a1ij2=upsample(a1i3,1,0);
 [bn3,an3] =butter(4, 0.5) ;
 aii33=filter(bn3,an3,a1ij2);
 
  d1ij3=upsample(di13,1,0);
 [bn3,an3] =butter(4, 0.5, 'high') ;
 dii33=filter(bn3,an3,d1ij3);
 output3=1.5*[aii33 dii33];
  y3=idwt(aii33(1:1:1000),dii33,'haar');
 [v3 z3]=alignsignals(y3,desired,20);
 errorVector3=z3-v3(1:1:1012);

[k j]=alignsignals(y3,input3,20);
err3(k1)=0;
for i=1:1:1000,
    err3(k1)=err3(k1)+((j(i)-k(i)).^2);
end
err_sureshrink3(k1)=err3(k1);

prd_sureshrink3(k1)=0;
sum13=0;
sum23=0;
for i=1:1:nIterations-5,
    sum13=sum13+(k(i)-j(i)).^2;
    sum23=sum23+(j(i)).^2;
end
prd_sureshrink3(k1)=100*sqrt(sum13/sum23);

snre_sureshrink3(k1)=0;
sum13=0;
sum23=0;
for i=1:1:nIterations-5,
    sum13=sum13+(y3(i).^2);
    sum23=sum23+((j(i)-k(i)).^2);
end
snre_sureshrink3(k1)=10*log10(sum13/sum23);

%dwt neighblock
lambdac3=4.505;
z3=(log10(5)/2)+(2*max(1,log10(5)/4));
sum3=0;
i=1;
j=1;
s=0;
[c3,d3]=dwt(input3,'haar');
while i<=495
   ci= c3(i:1:i+5).^2;
    for a3=1:5,
        s=s+ci(a3);
    end
    while j<=i+5
        c1(j)=(max(0,(s.^2-lambdac3*z3)/s.^2))*c3(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end
i=1;
j=1;
s=0;
while i<=495
    di=d3(i:1:i+5).^2;
    for a3=1:5,
        s=s+di(a3);
    end
    while j<=i+5
        d3(j)=(max(0,(s.^2-lambdac3*z3)/s.^2))*d3(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end

output3=1.095*idwt(c3,d3,'haar');

% errorVector3=desired(1:1:996)-output3;


err3(k1)=0;
for i=1:1:nIterations-5,
    err3(k1)=err3(k1)+((input3(i)-output3(i)).^2);
end
err_neighblock3(k1)=err3(k1);

prd_neighblock3(k1)=0;
sum13=0;
sum23=0;
for i=1:1:nIterations-5,
    sum13=sum13+(input3(i)-output3(i)).^2;
    sum23=sum23+(input3(i)).^2;
end
prd_neighblock3(k1)=100*sqrt(sum13/sum23);
snre_neighblock3(k1)=0;
sum13=0;
sum23=0;
for i=1:1:nIterations-5,
    sum13=sum13+(output3(i).^2);
    sum23=sum23+((input3(i)-output3(i)).^2);
end
snre_neighblock3(k1)=10*log10(sum13/sum23);

 end
%%
clc;
x=0.01:0.01:10;
 li=30/72;  

  a_pwav=0.25;
      d_pwav=0.09;
      t_pwav=0.16;  
     
      a_qwav=0.025;
      d_qwav=0.066;
      t_qwav=0.166;
      
      a_qrswav=1.6;
      d_qrswav=0.11;
      
      a_swav=0.25;
      d_swav=0.066;
      t_swav=0.09;
      
      a_twav=0.35;
      d_twav=0.142;
      t_twav=0.2;
      
      a_uwav=0.035;
      d_uwav=0.0476;
      t_uwav=0.433;
      
       pwav=p_wav(x,a_pwav,d_pwav,t_pwav,li);
       qwav=q_wav(x,a_qwav,d_qwav,t_qwav,li);
       qrswav=qrs_wav(x,a_qrswav,d_qrswav,li);
       swav=s_wav(x,a_swav,d_swav,t_swav,li);
       twav=t_wav(x,a_twav,d_twav,t_twav,li);
       uwav=u_wav(x,a_uwav,d_uwav,t_uwav,li);
       
 ecg=pwav+qrswav+twav+swav+qwav+uwav;
 namp=[0.6184 0.3475 0.1956 0.1100 0.0618 0.0348 0.0196 0.0110];
 for k1=1:8,
 noise4=namp(k1)*cos(2*pi*100*x);
% noise = wgn(100,1,0);
input4=noise4+ecg;

snr4(k1)= mean( ecg .^ 2 ) / mean( noise4 .^ 2 );
snr_db4(k1) = 10 * log10( snr4(k1) );

desired=ecg;
filterOrder=5;
nCoeff=filterOrder+1;
nIterations=length(desired);
errorVector4=zeros(nIterations,1);
outputVector=zeros(nIterations,1);
coeffVector=zeros(nCoeff,(nIterations+1));
S_d=eye(nCoeff);
p_d=zeros(nCoeff,1);
xk=zeros(nCoeff-1,1);
lambda4=0.91398;
final=zeros(10,1);

% dwt sureshrink algorithm

[c4,d4]=dwt(input4,'haar');
[b4,a4] =butter(4, 0.5, 'high') ;
dii4=filter(b4,a4,[c4,d4]);
di14=downsample(dii4,1,0);
 [b4,a4] =butter(4, 0.5) ;
 aii4=filter(b4,a4,[c4,d4]);
 ai14=downsample(aii4,1,0);
 
 [b4,a4] =butter(4, 0.5, 'high') ;
dii4=filter(b4,a4,ai14);
di24=downsample(dii4,1,0);
 [b4,a4] =butter(4, 0.5) ;
 aii4=filter(b4,a4,ai14);
 ai24=downsample(aii4,1,0);
 
 p4=2*log(5);
 lambda4=sqrt(p4);
 
 for it=1:500,
     if di14(it)>lambda4
         di14(it)=di14(it)-lambda4;
     elseif di14(it)<lambda4
         di14(it)=di14(it)+lambda4;
     else
         di14(it)=0;
     end
 end
 
 
 for it=1:250
       if di24(it)>lambda4
         di24(it)=di24(it)-lambda4;
     elseif di24(it)<lambda4
         di24(it)=di24(it)+lambda4;
     else
         di24(it)=0;
     end
 end
 dij4=upsample(di24,1,0);
 [bn4,an4] =butter(4, 0.5, 'high') ;
 dii4=filter(bn4,an4,dij4);
 
 aij4=upsample(ai24,1,0);
 [bn4,an4] =butter(4, 0.5) ;
 aii4=filter(bn4,an4,aij4);
 
 a1i4=[aii4 dii4];
  a1ij4=upsample(a1i4,1,0);
 [bn4,an4] =butter(4, 0.5) ;
 aii34=filter(bn4,an4,a1ij4);
 
  d1ij4=upsample(di14,1,0);
 [bn4,an4] =butter(4, 0.5, 'high') ;
 dii34=filter(bn4,an4,d1ij4);
 output4=1.5*[aii34 dii34];
  y4=idwt(aii34(1:1:1000),dii34,'haar');
 [v4 z4]=alignsignals(y4,desired,20);
 errorVector4=z4-v4(1:1:1012);

[k j]=alignsignals(y4,input4,20);
err4(k1)=0;
for i=1:1:1000,
    err4(k1)=err4(k1)+((j(i)-k(i)).^2);
end
err_sureshrink4(k1)=err4(k1);

prd_sureshrink4(k1)=0;
sum14=0;
sum24=0;
for i=1:1:nIterations-5,
    sum14=sum14+(k(i)-j(i)).^2;
    sum24=sum24+(j(i)).^2;
end
prd_sureshrink4(k1)=100*sqrt(sum14/sum24);

snre_sureshrink4(k1)=0;
sum14=0;
sum24=0;
for i=1:1:nIterations-5,
    sum14=sum14+(y4(i).^2);
    sum24=sum24+((j(i)-k(i)).^2);
end
snre_sureshrink4(k1)=10*log10(sum14/sum24);

%dwt neighblock
lambdac4=4.505;
z4=(log10(5)/2)+(2*max(1,log10(5)/4));
sum4=0;
i=1;
j=1;
s=0;
[c4,d4]=dwt(input4,'haar');
while i<=495
   ci= c4(i:1:i+5).^2;
    for a4=1:5,
        s=s+ci(a4);
    end
    while j<=i+5
        c1(j)=(max(0,(s.^2-lambdac4*z4)/s.^2))*c4(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end
i=1;
j=1;
s=0;
while i<=495
    di=d4(i:1:i+5).^2;
    for a4=1:5,
        s=s+di(a4);
    end
    while j<=i+5
        d1(j)=(max(0,(s.^2-lambdac4*z4)/s.^2))*d4(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end

output4=1.095*idwt(c1,d1,'haar');

errorVector4=desired(1:1:996)-output4;


err4(k1)=0;
for i=1:1:nIterations-5,
    err4(k1)=err4(k1)+((input4(i)-output4(i)).^2);
end
err_neighblock4(k1)=err4(k1);

prd_neighblock4(k1)=0;
sum14=0;
sum24=0;
for i=1:1:nIterations-5,
    sum14=sum14+(input4(i)-output4(i)).^2;
    sum24=sum24+(input4(i)).^2;
end
prd_neighblock4(k1)=100*sqrt(sum14/sum24);
snre_neighblock4(k1)=0;
sum14=0;
sum24=0;
for i=1:1:nIterations-5,
    sum14=sum14+(output4(i).^2);
    sum24=sum24+((input4(i)-output4(i)).^2);
end
snre_neighblock4(k1)=10*log10(sum14/sum24);

 end
%%
clc;
x=0.01:0.01:10;
 li=30/72;  

  a_pwav=0.25;
      d_pwav=0.09;
      t_pwav=0.16;  
     
      a_qwav=0.025;
      d_qwav=0.066;
      t_qwav=0.166;
      
      a_qrswav=1.6;
      d_qrswav=0.11;
      
      a_swav=0.25;
      d_swav=0.066;
      t_swav=0.09;
      
      a_twav=0.35;
      d_twav=0.142;
      t_twav=0.2;
      
      a_uwav=0.035;
      d_uwav=0.0476;
      t_uwav=0.433;
      
       pwav=p_wav(x,a_pwav,d_pwav,t_pwav,li);
       qwav=q_wav(x,a_qwav,d_qwav,t_qwav,li);
       qrswav=qrs_wav(x,a_qrswav,d_qrswav,li);
       swav=s_wav(x,a_swav,d_swav,t_swav,li);
       twav=t_wav(x,a_twav,d_twav,t_twav,li);
       uwav=u_wav(x,a_uwav,d_uwav,t_uwav,li);
       
 ecg=pwav+qrswav+twav+swav+qwav+uwav;
 namp=[0.6184 0.3475 0.1956 0.1100 0.0618 0.0348 0.0196 0.0110];
 for k1=1:8,
 noise5=namp(k1)*cos(2*pi*150*x);
% noise = wgn(100,1,0);
input5=noise5+ecg;

snr5(k1)= mean( ecg .^ 2 ) / mean( noise5 .^ 2 );
snr_db5(k1) = 10 * log10( snr5(k1) );

desired=ecg;
filterOrder=5;
nCoeff=filterOrder+1;
nIterations=length(desired);
errorVector5=zeros(nIterations,1);
outputVector=zeros(nIterations,1);
coeffVector=zeros(nCoeff,(nIterations+1));
S_d=eye(nCoeff);
p_d=zeros(nCoeff,1);
xk=zeros(nCoeff-1,1);
lambda5=0.91398;
final=zeros(10,1);

% dwt sureshrink algorithm

[c5,d5]=dwt(input5,'haar');
[b5,a5] =butter(4, 0.5, 'high') ;
dii5=filter(b5,a5,[c5,d5]);
di15=downsample(dii5,1,0);
 [b5,a5] =butter(4, 0.5) ;
 aii5=filter(b5,a5,[c5,d5]);
 ai15=downsample(aii5,1,0);
 
 [b5,a5] =butter(4, 0.5, 'high') ;
dii5=filter(b5,a5,ai15);
di25=downsample(dii5,1,0);
 [b5,a5] =butter(4, 0.5) ;
 aii5=filter(b5,a5,ai15);
 ai25=downsample(aii5,1,0);
 
 p5=2*log(5);
 lambda5=sqrt(p5);
 
 for it=1:500,
     if di15(it)>lambda5
         di15(it)=di15(it)-lambda5;
     elseif di15(it)<lambda5
         di15(it)=di15(it)+lambda5;
     else
         di15(it)=0;
     end
 end
 
 
 for it=1:250
       if di25(it)>lambda5
         di25(it)=di25(it)-lambda5;
     elseif di25(it)<lambda5
         di25(it)=di25(it)+lambda5;
     else
         di25(it)=0;
     end
 end
 dij5=upsample(di25,1,0);
 [bn5,an5] =butter(4, 0.5, 'high') ;
 dii5=filter(bn5,an5,dij5);
 
 aij5=upsample(ai25,1,0);
 [bn5,an5] =butter(4, 0.5) ;
 aii5=filter(bn5,an5,aij5);
 
 a1i5=[aii5 dii5];
  a1ij5=upsample(a1i5,1,0);
 [bn5,an5] =butter(4, 0.5) ;
 aii35=filter(bn5,an5,a1ij5);
 
  d1ij5=upsample(di15,1,0);
 [bn5,an5] =butter(4, 0.5, 'high') ;
 dii35=filter(bn5,an5,d1ij5);
 output5=1.5*[aii35 dii35];
  y5=idwt(aii35(1:1:1000),dii35,'haar');
 [v5 z5]=alignsignals(y5,desired,20);
 errorVector5=z5-v5(1:1:1012);

[k j]=alignsignals(y5,input5,20);
err5(k1)=0;
for i=1:1:1000,
    err5(k1)=err5(k1)+((j(i)-k(i)).^2);
end
err_sureshrink5(k1)=err5(k1);

prd_sureshrink5(k1)=0;
sum15=0;
sum25=0;
for i=1:1:nIterations-5,
    sum15=sum15+(k(i)-j(i)).^2;
    sum25=sum25+(j(i)).^2;
end
prd_sureshrink5(k1)=100*sqrt(sum15/sum25);

snre_sureshrink5(k1)=0;
sum15=0;
sum25=0;
for i=1:1:nIterations-5,
    sum15=sum15+(y5(i).^2);
    sum25=sum25+((j(i)-k(i)).^2);
end
snre_sureshrink5(k1)=10*log10(sum15/sum25);

%dwt neighblock
lambdac5=4.505;
z5=(log10(5)/2)+(2*max(1,log10(5)/4));
sum5=0;
i=1;
j=1;
s=0;
[c5,d5]=dwt(input5,'haar');
while i<=495
   ci= c5(i:1:i+5).^2;
    for a5=1:5,
        s=s+ci(a5);
    end
    while j<=i+5
        c1(j)=(max(0,(s.^2-lambdac5*z5)/s.^2))*c5(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end
i=1;
j=1;
s=0;
while i<=495
    di=d5(i:1:i+5).^2;
    for a5=1:5,
        s=s+di(a5);
    end
    while j<=i+5
        d1(j)=(max(0,(s.^2-lambdac5*z5)/s.^2))*d5(j);
        j=j+1;
    end
    i=i+4;
    s=0;
end

output5=1.095*idwt(c1,d1,'haar');

errorVector5=desired(1:1:996)-output5;


err5(k1)=0;
for i=1:1:nIterations-5,
    err5(k1)=err5(k1)+((input5(i)-output5(i)).^2);
end
err_neighblock5(k1)=err5(k1);

prd_neighblock5(k1)=0;
sum15=0;
sum25=0;
for i=1:1:nIterations-5,
    sum15=sum15+(input5(i)-output5(i)).^2;
    sum25=sum25+(input5(i)).^2;
end
prd_neighblock5(k1)=100*sqrt(sum15/sum25);
snre_neighblock5(k1)=0;
sum15=0;
sum25=0;
for i=1:1:nIterations-5,
    sum15=sum15+(output5(i).^2);
    sum25=sum25+((input5(i)-output5(i)).^2);
end
snre_neighblock5(k1)=10*log10(sum15/sum25);

 end

 figure;
 subplot(2,2,1)
 hold on;
plot(snr_db1,prd_neighblock1,'red');
hold on;
plot(snr_db1,prd_sureshrink1);
hold off;
xlim([10,40]);
ylim([0,25]);
legend('Neighblock','Sureshrink');
xlabel('SNR (dB)','FontSize',12);
ylabel('Estimated prd (%)','FontSize',12);
title(' noise signal of frequency 10 Hz','FontSize',12.5);
subplot(2,2,2)
plot(snr_db2,prd_neighblock2,'red');
hold on;
plot(snr_db2,prd_sureshrink2);
hold off;
xlim([10,40]);
ylim([0,25]);
legend('Neighblock','Sureshrink');
xlabel('SNR (dB)','FontSize',12);
ylabel('Estimated prd (%)','FontSize',12);
title(' noise signal of frequency 50 Hz','FontSize',12.5);
% figure
subplot(2,2,3)
plot(snr_db3,prd_neighblock3,'red');
hold on;
plot(snr_db3,prd_sureshrink3);
hold off;
xlim([10,40]);
ylim([0,25]);
legend('Neighblock','Sureshrink');
xlabel('SNR (dB)','FontSize',12);
ylabel('Estimated prd (%)','FontSize',12);
title('noise signal of frequency 80 Hz','FontSize',12.5);
subplot(2,2,4)
plot(snr_db4,prd_neighblock4,'red');
hold on;
plot(snr_db4,prd_sureshrink4);
hold off;
xlim([10,40]);
ylim([0,25]);
xlabel('SNR (dB)','FontSize',12);
ylabel('Estimated prd (%)','FontSize',12);
title(' noise signal of frequency 100 Hz','FontSize',12.5);
suptitle('Estimated PRD Vs SNR (dB) for Neighblock and Sureshrink algorithm with differnt frequencies of noise signal');
legend('Neighblock','Sureshrink');

 figure;
 subplot(2,2,1)
 hold on;
plot(snr_db1,snre_neighblock1,'red');
hold on;
plot(snr_db1,snre_sureshrink1);
hold off;
xlim([10,40]);
ylim([0,25]);
legend('Neighblock','Sureshrink');
xlabel('SNR (dB)','FontSize',12);
ylabel('Estimated SNR (dB)','FontSize',12);
title('noise signal of frequency 10 Hz','FontSize',12.5);
subplot(2,2,2)
plot(snr_db2,snre_neighblock2,'red');
hold on;
plot(snr_db2,snre_sureshrink2);
hold off;
xlim([10,40]);
ylim([0,25]);
legend('Neighblock','Sureshrink');
xlabel('SNR (dB)','FontSize',12);
ylabel('Estimated SNR (dB)','FontSize',12);
title(' noise signal of frequency 50 Hz','FontSize',12.5);
% figure
subplot(2,2,3)
plot(snr_db3,snre_neighblock3,'red');
hold on;
plot(snr_db3,snre_sureshrink3);
hold off;
xlim([10,40]);
ylim([0,25]);
legend('Neighblock','Sureshrink');
xlabel('SNR (dB)','FontSize',12);
ylabel('Estimated SNR dB)','FontSize',12);
title(' noise signal of frequency 80 Hz','FontSize',12.5);
subplot(2,2,4)
plot(snr_db4,snre_neighblock4,'red');
hold on;
plot(snr_db4,snre_sureshrink4);
hold off;
xlim([10,40]);
ylim([10,30]);

xlabel('SNR (dB)','FontSize',12);
ylabel('Estimated SNR dB)','FontSize',12);
title('noise signal of frequency 100 Hz','FontSize',12.5);

suptitle('Estimated SNR Vs SNR (dB)for Neighblock and Sureshrink algorithm with differnt frequencies of noise signal');
%  figure;
legend('Neighblock','Sureshrink');