clc
clear all;
close all
SampFreq = 128;
delta=4;
num=2;
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('D:\TFSA7\TFSA7');%
t = 0:1/SampFreq:1-1/SampFreq;

code=3;


if code==1  % well separated signals
   Sig5 =1*exp(1i*(2*pi*(15*t +15*t.^3)));

    Sig6 =1*exp(1i*(2*pi*(5*t +1*15*t.^3)));
    
    SigA=Sig5+Sig6;

    IF_O(:,1)=45*t.^2+15;
    IF_O(:,2)=45*t.^2+5;
    num=2;

elseif code==2 % Two component signal
   Sig5 =1*exp(1i*(2*pi*(15*t +15*t.^3)));

    Sig6 =1*exp(1i*(2*pi*(50*t -1*13*t.^3)));
    
    SigA=Sig5+Sig6;

    IF_O(:,1)=45*t.^2+15;
    IF_O(:,2)=-39*t.^2+50;
    
    num=2;
elseif code==3  % three component LFM signals
    Sig7=1*exp(1i*(2*pi*(20*t.^2))+1i*(2*pi*(10*t))); 

    Sig1 = 1*exp(1i*(-2*pi*(1*20*t.^2))+1i*(2*pi*(50*t))); %300tªÚ’ﬂ150t

    Sig2 = 1*exp(1i*(2*pi*(62*t -20*t.^2)));

    Sig3 = exp(1i*(2*pi*(5*t )));

    SigA =1*Sig1 +1*Sig2+1*Sig7+0*Sig3;

    IF_O(:,1)=-40*t+50;
    IF_O(:,2)=-40*t+62;
    IF_O(:,3)=40*t+10;
    %IF_O(:,4)=5;
    
    num=3;
elseif code==4 % 2 NLFM and 1 tone
    Sig1 =1*exp(1i*(2*pi*(5*t +15*t.^3)));
    Sig2 =1*exp(1i*(2*pi*(15*t +15*t.^3)));
    Sig3 =1*exp(1i*(2*pi*(50*t -1*13*t.^3)));
    Sig4=exp(1i*(2*pi*(60*t )));
    IF_O(:,1)=45*t.^2+15;
    IF_O(:,2)=45*t.^2+5;
    IF_O(:,3)=-39*t.^2+50;
    %IF_O(:,3)=60;
    SigA=1*Sig1+1*Sig2+1*Sig3+0*Sig4;
    num=3;
end
IF_O=IF_O/(SampFreq/2);
 %Sig=awgn(SigA,5,'measured');
Sig=SigA;

f = linspace(0,SampFreq/2,length(Sig));
[tfd,orient]=HTFD_neww(Sig,2,30,84);
imagesc(t,f,tfd); 
% axis([0 1 -SampFreq/2 SampFreq/2]);
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);


[IF,out, peaks] = component_linking_new(tfd,orient,0.02,length(tfd)/8,10);
[IF]= merge_IFs(IF,orient,20,30,length(Sig)/2);
figure;
imagesc(t,f,peaks); 
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

IF_est=fill_zeros(IF);

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

figure;
imagesc(t,f,IF_image); 
% axis([0 1 -SampFreq/2 SampFreq/2]);
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);


figure;% plot(IF_est.')
plot(t(1:4:end),IF_est(:,1:4:end).'/2,'o',t,IF_O'.*64,'g','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
%title('(b)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20)
set(gca,'linewidth',2);
axis([0 1 0 64]);

[IF,out, peaks] = component_linking(tfd,0.1,64);
 %figure; imagesc(out)
figure;
imagesc(t,f,out); 
% axis([0 1 -SampFreq/2 SampFreq/2]);
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

 
 
 
 IF_est=fill_zeros(IF);
 figure;% plot(IF_est.')
plot(t(1:4:end),IF_est(:,1:4:end).'/2,'o',t,IF_O'.*64,'g','linewidth',2);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
%title('(c)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20)
set(gca,'linewidth',2);
axis([0 1 0 64]);

