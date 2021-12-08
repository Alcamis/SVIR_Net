clear
time_start=tic;
clc
close all

a1=0; % strength of the coupling term for S
a2=0.01 % strength of the coupling term for I

t0=0;
tf=500;
t_step=0.1;
t_span=t0:t_step:tf;

filename='adjacency_matrix_Barabasi_network.txt'; % connectivity matrix
L=readmatrix(filename,'Delimiter',' ');

N=size(L,1); % number of communities
[laplacian_L,K]=Laplacian_matrix(L);

% mu:       births occur with the same constant rate mu of deaths

% beta:     the transmission rate of susceptible to infected individuals

% phi:      0 <= phi <= 100, the percentage of the vaccinated population in the susceptible population

% rho:      0 <= rho <= 1, the effectiveness of vaccination in the community
%           rho=0 means that the vaccine is perfectly effective
%           rho=1 means that the vaccine has no effect

% lambda:   the rate of recovered individuals in the community

% delta:    the rate for the loss of immunity

% theta:    the rate by which vaccination loses effect

% Bronze community
mu_b=0;
beta_b=0.02;
phi_b=1e-16;
rho_b=0;
lambda_b=0.01;
delta_b=0.0001;
theta_b=0;

% Silver community
mu_s=0;
beta_s=0.02;
phi_s=40;
rho_s=0.6;
lambda_s=0.01;
delta_s=0.0001;
theta_s=0;

% Gold community
mu_g=0;
beta_g=0.02;
phi_g=80;
rho_g=0.3;
lambda_g=0.01;
delta_g=0.0001;
theta_g=0;

Bronze=[mu_b beta_b phi_b rho_b lambda_b delta_b theta_b];
Silver=[mu_s beta_s phi_s rho_s lambda_s delta_s theta_s];
Gold= [mu_g beta_g phi_g rho_g lambda_g delta_g theta_g];

C=zeros(N,7);

cbn=5;
csn=30;
cgn=N;

C(1:cbn,:)=repmat(Bronze,cbn,1);
C(cbn+1:csn,:)=repmat(Silver,csn-cbn,1);
C(csn+1:cgn,:)=repmat(Gold,cgn-csn,1);

xic=zeros(4*N,1);
for j=1:N
    xic(1+4*(j-1))=1-1e-16; % S
    xic(2+4*(j-1))=0;       % V
    xic(3+4*(j-1))=1e-16;   % I
    xic(4+4*(j-1))=0;       % R
end

options=odeset('NonNegative',[1:4*N]);
[T,Y]=ode45(@(t,y) eqs_SVIR_coupled(a1,a2,N,t,y,C,laplacian_L),t_span,xic,options);







xx=[1:N]';
yy=strsplit(num2str(xx(:)'));

BIs=Y(:,3:4:4*N); % all Is


M=sum(Y,2); % total population: M = Mbronze + Msilver + Mgold

Mbronze=sum(Y(:,1:4*cbn),2); % total population of Bronze communities
Msilver=sum(Y(:,4*cbn+1:4*csn),2); % total population of Silver communities
Mgold=sum(Y(:,4*csn+1:4*cgn),2); % total population of Gold communities




figure(1)
set(gca,'fontsize',15)
colormap jet
imagesc(xx,T,BIs);
h=colorbar;
ylabel(h, '$I$','fontsize',15,'Interpreter','latex')

xticks([1 0 0 0 0 0 0 0 0 0 0 0 0]+[0:5:N])

caxis([0 max(max(BIs))])
xlabel('Community index','fontsize',15)
ylabel('$t$','fontsize',15,'Interpreter','latex')
caxis([1e-16 1])
set(gca,'ColorScale','log')
set(gca,'YDir','normal')

ax=gca; 
ax.YAxis.Exponent=0;
ax.FontSize=15;

xline(cbn,'--w','LineWidth',4)
xline(csn,'--w','LineWidth',4)

text(0.0,1.05,"Bronze",'Units','normalized','HorizontalAlignment','center','fontsize',15)
text(0.24,1.05,"Silver",'Units','normalized','HorizontalAlignment','center','fontsize',15)
text(0.24+.5,1.05,"Gold",'Units','normalized','HorizontalAlignment','center','fontsize',15)

pos = get(gca,'Position');
xoffset = 0.02;
pos(1) = pos(1) + xoffset;
set(gca, 'Position', pos)

set(h,'Ticks', [1e-16 1e-12 1e-8 1e-4 1], 'TickLabels',["10^{-16}" "10^{-12}" "10^{-8}" "10^{-4}" 1])

fileID = fopen(strcat('Is_vs_time_spatiotemporal_plot_N=',string(N),'_a1=',string(a1),'_a2=',string(a2),'_Barabasi_',string(cbn),'bronze_',string(csn-cbn),'silver_',string(cgn-csn),'gold_info.txt'),'w');
fprintf(fileID,'a1=%26.16f\n',a1);
fprintf(fileID,'a2=%26.16f\n',a2);
fprintf(fileID,'t0=%26.16f\n',t0);
fprintf(fileID,'tf=%26.16f\n',tf);
fprintf(fileID,'t_step=%26.16f\n',t_step);
fprintf(fileID,'N=%d\n',N);

fprintf(fileID,'\n');
fprintf(fileID,'mu:       births occur with the same constant rate mu of deaths\n');
fprintf(fileID,'beta:     the transmission rate of susceptible to infected individuals\n');
fprintf(fileID,'phi:      0 <= phi <= 100, the percentage of the vaccinated population in the susceptible population\n');
fprintf(fileID,'rho:      0 <= rho <= 1, the effectiveness of vaccination in the community\n');
fprintf(fileID,'          rho=0 means that the vaccine is perfectly effective\n');
fprintf(fileID,'          rho=1 means that the vaccine has no effect\n');
fprintf(fileID,'lambda:   the rate of recovered individuals in the community\n');
fprintf(fileID,'delta:    the rate for the loss of immunity\n');
fprintf(fileID,'theta:    the rate by which vaccination loses effect\n');

fprintf(fileID,'\n');
fprintf(fileID,'Bronze communities: From 1 to %d\n',cbn);
fprintf(fileID,'mu=%26.16f\n',mu_b);
fprintf(fileID,'beta=%26.16f\n',beta_b);
fprintf(fileID,'phi=%26.16f\n',phi_b);
fprintf(fileID,'rho=%26.16f\n',rho_b);
fprintf(fileID,'lambda=%26.16f\n',lambda_b);
fprintf(fileID,'delta=%26.16f\n',delta_b);
fprintf(fileID,'theta=%26.16f\n',theta_b);

fprintf(fileID,'\n');
fprintf(fileID,'Silver communities: From %d to %d\n',cbn+1,csn);
fprintf(fileID,'mu=%26.16f\n',mu_s);
fprintf(fileID,'beta=%26.16f\n',beta_s);
fprintf(fileID,'phi=%26.16f\n',phi_s);
fprintf(fileID,'rho=%26.16f\n',rho_s);
fprintf(fileID,'lambda=%26.16f\n',lambda_s);
fprintf(fileID,'delta=%26.16f\n',delta_s);
fprintf(fileID,'theta=%26.16f\n',theta_s);

fprintf(fileID,'\n');
fprintf(fileID,'Gold communities: From %d to %d\n',csn+1,cgn);
fprintf(fileID,'mu=%26.16f\n',mu_g);
fprintf(fileID,'beta=%26.16f\n',beta_g);
fprintf(fileID,'phi=%26.16f\n',phi_g);
fprintf(fileID,'rho=%26.16f\n',rho_g);
fprintf(fileID,'lambda=%26.16f\n',lambda_g);
fprintf(fileID,'delta=%26.16f\n',delta_g);
fprintf(fileID,'theta=%26.16f\n',theta_g);

fprintf(fileID,'\n');
fprintf(fileID,'cbn=%d\n',cbn);
fprintf(fileID,'csn=%d\n',csn);
fprintf(fileID,'cgn=%d\n',cgn);

for j=1:N
    fprintf(fileID,'\n');
    fprintf(fileID,'xic(%d)=%26.16f S%d\n',1+4*(j-1),xic(1+4*(j-1)),j); % S
    fprintf(fileID,'xic(%d)=%26.16f V%d\n',2+4*(j-1),xic(2+4*(j-1)),j); % V
    fprintf(fileID,'xic(%d)=%26.16f I%d\n',3+4*(j-1),xic(3+4*(j-1)),j); % I
	fprintf(fileID,'xic(%d)=%26.16f R%d\n',4+4*(j-1),xic(4+4*(j-1)),j); % R
end

saveas(gcf,strcat('Is_vs_time_spatiotemporal_plot_N=',string(N),'_a1=',string(a1),'_a2=',string(a2),'_Barabasi_',string(cbn),'bronze_',string(csn-cbn),'silver_',string(cgn-csn),'gold.eps'),'epsc');
savefig(strcat('Is_vs_time_spatiotemporal_plot_N=',string(N),'_a1=',string(a1),'_a2=',string(a2),'_Barabasi_',string(cbn),'bronze_',string(csn-cbn),'silver_',string(cgn-csn),'gold.fig'));
set(gca,'ColorScale','log')

save(strcat('Is_vs_time_spatiotemporal_plot_N=',string(N),'_a1=',string(a1),'_a2=',string(a2),'_Barabasi_',string(cbn),'bronze_',string(csn-cbn),'silver_',string(cgn-csn),'gold.mat'))
set(gca,'ColorScale','log')





figure(2)
semilogy(T,Y(:,1),'g-' ,T,Y(:,2),'b-' ,T,Y(:,3),'r:',T,Y(:,4),'c-.','LineWidth',1.5)
legend('$S_1(t)$','$V_1(t)$','$I_1(t)$','$R_1(t)$','fontsize',15,'Interpreter','latex','Location','northeastoutside')
xlabel('$t$','fontsize',15,'Interpreter','latex')
ylim([1e-25 1e2])

%                 [left bottom width height]
% legend('Position',[0.78 0.55 0.1 0.2])
legend('boxoff')

ax=gca;
ax.FontSize=15;

saveas(gcf,strcat('Bronze_BN.eps'),'epsc');
savefig(strcat('Bronze_BN.fig'));





figure(3)
semilogy(T,Y(:,1+4*cbn),'g-' ,T,Y(:,2+4*cbn),'b-' ,T,Y(:,3+4*cbn),'r:',T,Y(:,4+4*cbn),'c-.','LineWidth',1.5)
legend('$S_6(t)$','$V_6(t)$','$I_6(t)$','$R_6(t)$','fontsize',15,'Interpreter','latex','Location','northeastoutside')
xlabel('$t$','fontsize',15,'Interpreter','latex')
ylim([1e-25 1e2])

%                 [left bottom width height]
% legend('Position',[0.78 0.55 0.1 0.2])
legend('boxoff')

ax=gca; 
ax.FontSize=15;

saveas(gcf,strcat('Silver_BN.eps'),'epsc');
savefig(strcat('Silver_BN.fig'));





figure(4)
semilogy(T,Y(:,1+4*csn),'g-' ,T,Y(:,2+4*csn),'b-' ,T,Y(:,3+4*csn),'r:',T,Y(:,4+4*csn),'c-.','LineWidth',1.5)
legend('$S_{31}(t)$','$V_{31}(t)$','$I_{31}(t)$','$R_{31}(t)$','fontsize',15,'Interpreter','latex','Location','northeastoutside')
xlabel('$t$','fontsize',15,'Interpreter','latex')
ylim([1e-25 1e2])

%                 [left bottom width height]
% legend('Position',[0.78 0.55 0.1 0.2])
legend('boxoff')

ax=gca; 
ax.FontSize=15;

saveas(gcf,strcat('Gold_BN.eps'),'epsc');
savefig(strcat('Gold_BN.fig'));





figure(5)
plot5=plot(T,M,'g-',T,repmat(N,size(T,1),1),'black-.',T,Mbronze,'g-',T,repmat(cbn,size(T,1),1),'black-.',T,Msilver,'b-',T,repmat(csn-cbn,size(T,1),1),'black-.',T,Mgold,'r',T,repmat(cgn-csn,size(T,1),1),'black-.');
plot5(1).LineWidth=1.5;
plot5(2).LineWidth=0.4;
plot5(3).LineWidth=1.5;
plot5(4).LineWidth=0.4;
plot5(5).LineWidth=1.5;
plot5(6).LineWidth=0.4;
plot5(7).LineWidth=1.5;
plot5(8).LineWidth=0.4;

legend('$M$','','$M_B$','','$M_S$','','$M_G$','','fontsize',15,'Interpreter','latex','Location','northeastoutside')
xlabel('$t$','fontsize',15,'Interpreter','latex')
ylim([0 70])

%                [left bottom width height]
% legend('Position',[0.78 0.55 0.1 0.2])
legend('boxoff')

ax=gca; 
ax.FontSize=15;

saveas(gcf,strcat('M_M_b_M_s_M_g_vs_time_BN.eps'),'epsc');
savefig(strcat('M_M_b_M_s_M_g_vs_time_BN.fig'));

fprintf('Elapsed time (in DD:HH:MM:SS.MS) [%s]\n',datestr(toc(time_start)/(24*60*60),'DD:HH:MM:SS.FFF'));
fprintf(fileID,'\nElapsed time (in DD:HH:MM:SS.MS) [%s]\n',datestr(toc(time_start)/(24*60*60),'DD:HH:MM:SS.FFF'));
fclose(fileID);