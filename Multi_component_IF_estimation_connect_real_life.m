clc
clear all;
%close all
delta=4;
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('D:\TFSA7\TFSA7');%



Sig=bat_signal();
SampFreq=143000;
f = linspace(0,SampFreq/2,length(Sig));
t = 0:1/SampFreq:1-1/SampFreq;
 %Sig=decimate(Sig,4);
%Sig=filter([1 -1],1,Sig);
% 
  n=0:length(Sig)-1;
  Sig=Sig.'+0.075*cos(2*pi*(0.1*n+2*0.00012*n.^2)).*gausswin(400).';
%   Sig=hilbert(Sig);
% t=0:1/16:8-1/16;
%[tfd,orient]=HTFD_new_dec(Sig,3,15,64,4);
%[tfd1,orient]=HTFD_neww(Sig,3,12,64);
[tfd,orient]=HTFD_neww(Sig,3,12,128);
%[tfd,orient]=HTFD_neww(Sig,2,40,64*2);
%tfd=min(tfd,tfd1);
%tfsapl(Sig,tfd);

%tfd=min(tfd,tfd1);
tfd(1:5,:)=0;
[IF,out, peaks] = component_linking_new(tfd,orient,0.005,length(Sig)/50,15);
%[IF,out, peaks] = component_linking_neww_direc(tfd,orient,0.01,length(Sig)/4,30,15);
[IF]= merge_IFs(IF,orient,30,50,length(Sig)/8);%[IF,out, peaks] = component_linking_neww_direc(tfd,orient,0.2,length(Sig)/4,30,5);

%figure;tfsapl(Sig,peaks)
%figure; tfsapl(Sig,out.*tfd)

IF_est=fill_zeros(IF);
plot(IF_est.');

[a,b]=size(IF_est);
IF_image=zeros(size(out));
IF_est=round(IF_est);
for ii=1:a
    for jj=1:b
        if IF_est(ii,jj)~=0
        IF_image(IF_est(ii,jj),jj)=1;
        end
    end
end

%tfsapl(Sig,IF_image.*tfd)
figure;imagesc(t,f,tfd); 
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure;imagesc(t,f,peaks); 
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);
figure;imagesc(t,f,IF_image); 
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);