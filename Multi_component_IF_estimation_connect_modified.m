clc
clear all;
close all
SampFreq = 128;
delta=4;
num=2;
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('D:\TFSA7\TFSA7');%
t = 0:1/SampFreq:1-1/SampFreq;

Sig7=1*exp(1i*(2*pi*(20*t.^2))+1i*(2*pi*(10*t))); 

Sig1 = 1*exp(1i*(-2*pi*(1*20*t.^2))+1i*(2*pi*(50*t))); %300tªÚ’ﬂ150t

Sig2 = 1*exp(1i*(2*pi*(62*t -20*t.^2)));

Sig3 = exp(1i*(2*pi*(5*t +0*t.^3)));





Sig4 =1*exp(1i*(2*pi*(5*t +15*t.^3)));
Sig5 =1*exp(1i*(2*pi*(15*t +15*t.^3)));

Sig6 =1*exp(1i*(2*pi*(50*t -1*13*t.^3)));
Sig8=exp(1i*(2*pi*(60*t )));

%Sig =1*Sig3 +1*Sig1 +0*Sig1+0*Sig2;%+0*Sig3;
%IF_O(:,3)=52*t.^0;

IF_O(:,2)=56*t.^2+5;
%IF_O(:,3)=45*t.^2+20;
%IF_O(:,4)=8*t.^0;
IF_O(:,1)=-60*t.^2+62;


num=2;
Sig =1*Sig5 +1*Sig6+0*Sig4+1*Sig8;
%Sig =1*Sig5 +1*Sig3;

% =0*Sig1 +1*Sig2+1*Sig7+1*Sig3;

% Sig2 =1*exp(1i*(2*pi*(15*t +15*t.^3)));
%     Sig3 =1*exp(1i*(2*pi*(50*t -1*13*t.^3)));
%     Sig4=exp(1i*(2*pi*(60*t )));
% 
 %Sig=awgn(Sig,5,'measured');




%[tfd,orient]=HTFD_neww(Sig,2,30,64);
 tfd = quadtfd(Sig, length(Sig)/4-1, 1, 'mb',0.1,length(Sig));
tfd(tfd<0)=0;
orient=orient_img_calc(tfd ,1, 1, 1);
orient=180*(orient)/(pi);
%tfd=tfd.^2;
tfsapl(Sig,tfd);
%[IF,out, peaks] = component_linking_neww_direc(tfd,orient,0.1,length(Sig)/2,15,8);
%figure; imagesc(out)
%orient(orient>90)=180-orient(orient>90);

[IF,out, peaks] = component_linking_new(tfd,orient,0.1,length(tfd)/8,10);
[IF]= merge_IFs(IF,orient,10,30);%[IF,out, peaks] = component_linking_neww_direc(tfd,orient,0.2,length(Sig)/4,30,5);
jj=1;

for ii=1:length(IF)
    if length(IF{ii})>length(Sig)/2
        IF1{jj}=IF{ii};
        jj=jj+1;
    end
end
IF=IF1;

figure;imagesc(peaks)
figure; imagesc(out)

IF_est=fill_zeros(IF);
figure; plot(IF_est.')

figure; imagesc(out)
IF_est=fill_zeros(IF);
figure;plot(IF_est.');

[a,b]=size(IF_est);
IF_image=zeros(size(out));
IF_est=round(IF_est);
%[IF_est,~] = RPRG(IF_est,20);
for ii=1:a
    for jj=1:b
        if IF_est(ii,jj)~=0
        IF_image(IF_est(ii,jj),jj)=1;
        end
    end
end
figure;imagesc(IF_image)

% % IF_est=zeros(length(IF),length(Sig));
% % for i=1:length(IF)
% %     A=IF{i};
% %     t_axis=A(:,2);
% %     f_axis=A(:,1);
% %     for tt=1:length(A)
% %         IF_est(i,t_axis(tt))=f_axis(tt);
% %     end
% %     IF_est1(i,:)=medfilt1( IF_est(i,:),round(length(Sig)/2));
% %     IF_est(i,IF_est(i,:)==0)=IF_est1(i,IF_est(i,:)==0);
% % end
% figure; plot(IF_est.')
% 
% [findex]  = extridge_SVD_IF_spec(Sig, num, 3);
% 
% %  [fidexmult,s_sig] = extridge_mult(Sig, 3, 5,3);
% % % %[fidexmult,s_sig] = extridge_mult_new_modified1(Sig, 2, 5);
% % % 
%  [findex,interset] = RPRG(findex,25);
%   figure; plot(findex.')

