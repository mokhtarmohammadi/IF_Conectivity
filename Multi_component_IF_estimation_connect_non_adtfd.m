clc
clear all;
close all
delta=4;
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('D:\TFSA7\TFSA7');%
SampFreq=128*2;
t = 0:1/SampFreq:1-1/SampFreq;

%Sig7=1*exp(1i*(2*pi*(20*t.^2))+1i*(2*pi*(10*t))); 

Sig1 = 1*exp(1i*(-0.5*pi*(1*SampFreq*t.^2))+1i*(pi*(SampFreq*t))); %300t或者150t

Sig2 = 1*exp(1i*(pi*(0*t +0.5*SampFreq*t.^2)));

Sig3 = exp(1i*(pi*(0.5*SampFreq*t)));

Sig4 = 1*exp(1i*(-1*pi*(1*SampFreq*t.^3)/3)+1i*(pi*(0.95*SampFreq*t))); %300t或者150t

Sig5 = 1*exp(1i*(pi*(0.05*SampFreq*t +SampFreq*t.^3/3)));


num=2;
Sig =1*Sig1 +1*Sig2+1*Sig3;
%Sig =1*Sig4 +1*Sig5+1*Sig3;
Sig =1*Sig1 +1*Sig2+1*Sig3;
%Sig =1*Sig3 +1*Sig2+1*Sig4;

  

f = linspace(0,SampFreq/2,length(Sig));


%[tfd,orient]=HTFD_neww(Sig,3,12,64);
     %   tfd = quadtfd(Sig, length(Sig)/4-1, 1, 'specx',length(Sig)/4-1,'hamm',length(Sig));
 %       tfd = quadtfd(Sig, length(Sig)/8-1, 8, 'mb',0.1,length(Sig)/8);
%[tfd,orient]=HTFD_neww(Sig,3,24,64*2);
 %  figure; tfsapl(Sig,tfd)
tic
 tfd = quadtfd(Sig, length(Sig)/4-1, 1, 'mb',0.2,length(Sig));
tfd(tfd<0)=0;
orient=orient_img_calc(tfd ,1, 1, 1);
orient=180*(orient)/(pi);
toc


% COde for real-life signals
tic
[tfd,orient]=HTFD_neww(Sig,2,30,64);
toc
%[IF,out, peaks] = component_linking_neww(tfd,orient,0.2,length(tfd)/2,40,20);
[IF,out, peaks] = component_linking_new(tfd,orient,0.02,length(tfd)/8,10);
[IF]= merge_IFs(IF,orient,20,30,length(Sig)/2);%[IF,out, peaks] = component_linking_neww_direc(tfd,orient,0.2,length(Sig)/4,30,5);

%figure;tfsapl(Sig,peaks)

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

%tfsapl(Sig,IF_image.*tfd)
figure;imagesc(tfd); 
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure;imagesc(peaks); 

set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);
figure;imagesc(IF_image); 

set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);