clc;
clear;
addpath(pwd,'m-files');
close all;lambda=1064;C1=[0.24,0.57,0.25];C2=[1,0.38,0];C3=[0.25,0.41,0.88];C4=[0.89,0.09,0.05];
path1=[pwd,'\lisa_quat_scatpath.txt'];
Cass1=Scatt_Para_der(path1,lambda*10^-3,-1,1,-1);
Calc=Roughness_Surface;
 Cass1.Plot_sactter(10,2,3,C1,1)
 %%
figure('Color',[1 1 1])
 h1=histogram(Cass1.deviation,30);
[PDF1]=PDFcal(h1);
figure('Color',[1 1 1])
area(PDF1(:,1),PDF1(:,2), 'LineWidth', 1.5,'EdgeColor',C1,'FaceColor',C1,'FaceAlpha',0.1)
 grid on
xlabel('{\theta_s}-{\theta_i} [бу]')
ylabel('Probability Distribution')
set(gca,'Fontname','Palatino Linotype')
set(gca,'Fontsize',15)
set(gcf,'Position',[510 500 500 390]);
set(gca,'Position',[.2 .16 .7 .75]);
%%
% figure('Color',[1 1 1])
% h1=histogram(Cass1.f,100);
fmin4=min(Cass1.f);
fmax4=max(Cass1.f);
%%
[ScatPSD_fit1,BandLimScatPSD,RMS]=PsdDet(Calc,5,fmin4,fmax4,lambda);
figure('color',[1 1 1])
loglog(ScatPSD_fit1(:,1),ScatPSD_fit1(:,2),'lineWidth', 1.5)
hold on
area(BandLimScatPSD(:,1),BandLimScatPSD(:,2), 'LineWidth', 1.5,'EdgeColor',C1,'FaceColor',C1,'FaceAlpha',0.1)
hold off
save('D:\baidupan\Ciomp\PHD\CODE\LISA\Total\SijingPSD2.txt','BandLimScatPSD','-ascii') 
%%
function [Cass]=Scatt_Para_der(path1,lambda,a1,a2,a3)
Cass=Scattering_cal;
Cass=Cass.load_data(path1);
Cass=Cass.Surface_para;
Cass=Cass.Space_para;
Cass=Cass.Scatter_para(0.74,-5.42,0.0117,a1,a2,a3);
Cass.deviation=Cass.scatter_para.theta_s-Cass.scatter_para.theta_i;
Cass.f=abs((sind(Cass.scatter_para.theta_s)-sind(Cass.scatter_para.theta_i))/lambda);
end

function [PDF1]=PDFcal(h1)
 PDF1(:,1)=h1.BinEdges(1:end);
 PDF1(:,2)=[h1.BinCounts/sum(h1.BinCounts),0];
end
function [t,y]=Histplot(Cass,step)
figure('Color',[1 1 1])
 h1=histogram(Cass.deviation,step);
 t=h1.BinEdges(1:end);y=[0,h1.BinCounts/size(Cass.deviation,1)];
close gcf
end
%%
function [ScatPSD]=Scatt_Psd_cal(ARS,AOI,theta_smin,theta_smax,step,Sign)
Qc=Opaque_Optical_factor;
Qc.theta_smin=theta_smin;Qc.theta_smax=theta_smax;
Qc.step=step;
Qc=Qc.Kienzle(AOI);
ScatPSD(:,1)=ARS.fx5;
ScatPSD(:,2)=ARS.AOI(:,2)./Qc.Qss';
if Sign>0
A=find(ARS.fx5>0);   
 ScatPSD=sortrows(abs(ScatPSD(A,:)),1);
elseif Sign<0
 A=find(ARS.fx5<0);   
 ScatPSD=sortrows(abs(ScatPSD(A,:)),1);
else
     ScatPSD=sortrows(abs(ScatPSD(A,:)),1);
end
X=find(ScatPSD(:,2)>1E10);
ScatPSD(X,:)=[];
% semilogy(ARS.AOI(:,1),Qc.Qss)
end

function [Num]=index_deter(A,index,lower,upper)
A1=abs(A(:,1)-index);
A2=abs(A(:,1)-index+lower);
A3=abs(A(:,1)-index-upper);
Num(1)=find(A1==min(A1));
Num(2)=find(A2==min(A2));
Num(3)=find(A3==min(A3));
end

function [ScatPSD_fit1,BandLimScatPSD,RMS]=PsdDet(Calc,AOI1,fmin,fmax,lambda)
theta_min=-85;theta_max=85;
step1=0.05;
BRDF_fit1.AOI(:,1)=theta_min:step1:theta_max;
BRDF_fit1.fx5=(sind(BRDF_fit1.AOI(:,1))-sind(AOI1))/(lambda*10^-3);
BRDF_fit1.AOI(:,2)=Harvey(AOI1,BRDF_fit1(:,1).AOI(:,1),0.74,-5.42,0.0117).*cosd(BRDF_fit1.AOI(:,1));
[ScatPSD_fit1]=Scatt_Psd_cal(BRDF_fit1,AOI1,theta_min,theta_max,step1,-1);
PSD1_cross_min=index_deter(ScatPSD_fit1,0.2057,0,0);
PSD1_cross_man=index_deter(ScatPSD_fit1,0.2991,0,0);
BandLimScatPSD=[ScatPSD_fit1(PSD1_cross_min(1):PSD1_cross_man(1),1),ScatPSD_fit1(PSD1_cross_min(1):PSD1_cross_man(1),2)];
RMS1=Calc.PSD2RMS(ScatPSD_fit1(20:end,1),ScatPSD_fit1(20:end,2),3);
RMSSijing=Calc.PSD2RMS(BandLimScatPSD(:,1),BandLimScatPSD(:,2),3);
RMS=[RMS1,RMSSijing];
end
