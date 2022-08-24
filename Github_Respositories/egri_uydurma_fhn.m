clc; clear;

% fhn osilator egri uydurma polinomal yontemi

syms tt

%step size 
h = 0.1;
N = 350-1;
t(1) = 0;

a_sp1 = 0.07; b_sp1 = 0.08; c_sp1 = 0.5; aa_sp1 = 0;   f_sp1=1;     I_sp1=0.15;
a_sp2 = 0.2; b_sp2 = 0.1; c_sp2 = 2; aa_sp2 = 0;   f_sp2=1;     I_sp2=0.05; 
a_sp3 = 0.35; b_sp3 = 0.4; c_sp3 = 2.5; aa_sp3 = 0;   f_sp3=1;     I_sp3=0.75; 
a_sp7 = 0.75; b_sp7 = 1.2; c_sp7 = 7; aa_sp7 = 0;   f_sp7=1;     I_sp7=0.55; 

v_sp1(1) = -2;      w_sp1(1) = 1.16388;
v_sp2(1) = -2;      w_sp2(1) = 1.16388;
v_sp3(1) = -2;      w_sp3(1) = 1.16388; 
v_sp7(1) = -2;      w_sp7(1) = 1.16388;

wet=1;

[v_sp1, w_sp1] = function_fhn_1n(a_sp1,b_sp1,c_sp1,I_sp1,h,N,t,v_sp1,w_sp1,aa_sp1,f_sp1,wet);
[v_sp2, w_sp2] = function_fhn_1n(a_sp2,b_sp2,c_sp2,I_sp2,h,N,t,v_sp2,w_sp2,aa_sp2,f_sp2,wet);
[v_sp3, w_sp3] = function_fhn_1n(a_sp3,b_sp3,c_sp3,I_sp3,h,N,t,v_sp3,w_sp3,aa_sp3,f_sp3,wet);
[v_sp7, w_sp7] = function_fhn_1n(a_sp7,b_sp7,c_sp7,I_sp7,h,N,t,v_sp7,w_sp7,aa_sp7,f_sp7,wet);

figure(1); clf(1);
plot(v_sp1(2,:),v_sp1(1,:),'k','LineWidth',2)
xlabel('Number of Samples')
ylabel('Membran Potential (mV)')
set(gca,'Fontsize',30)
grid on


x = v_sp1(1,:);
u = w_sp1(1,:);

% fhn_known = load('mat_fhn_sp135.mat');
% x = fhn_known.v_sp7;
% u = fhn_known.w_sp7;

for i=1:1:150  
    
    m = i; %polinomun derecesi
       
    t = 1:1:length(x);
    
    a = polyfit(t,x,m);
    b = polyfit(t,u,m);
    
    v_sym_t = poly2sym(a,tt); % zaman t ye göre v yi yeniden oluþturur.
    w_sym_t = poly2sym(b,tt);
    
    for k = 1:length(t)
    
        fhn_cal_v(k) = polyval(a,t(k));
        fhn_cal_w(k) = polyval(b,t(k));
    
    end
    
    m_s_e_v = immse(x,fhn_cal_v);
    m_a_e_v = mae(x,fhn_cal_v);
    r_m_s_e_v = sqrt(m_s_e_v);
    n_r_m_s_e_v = (r_m_s_e_v/(max(x) - min(x)))*100;
    corr = corr2(x,fhn_cal_v);
    
    m_s_e_w = immse(u,fhn_cal_w);
    m_a_e_w = mae(u,fhn_cal_w);
    r_m_s_e_w = sqrt(m_s_e_w);
    n_r_m_s_e_w = (r_m_s_e_w/(max(u) - min(u)))*100;
    corr = corr2(u,fhn_cal_w);

    NRMSE(i) = n_r_m_s_e_v;

end

[Min_vlue,Min_locs] = min(NRMSE);

fig1 = figure('Position',get(0,'Screensize'));
plot(NRMSE(1:80),'LineStyle','-','Marker','o','Color','k','MarkerSize',20,'linewidth',4)
grid on; hold on;
plot(Min_locs,Min_vlue,'Marker','*','Color','r','MarkerSize',30,'linewidth',6)
xline(Min_locs, 'Color', 'r', 'LineWidth', 4)
textLabel = sprintf('Min of %.4f NRMSE at Polynominal Order=%.0f', Min_vlue, Min_locs);
text(0, 130, textLabel, 'fontSize', 50, 'Color', 'r', 'VerticalAlignment','middle')
ylabel('NRMSE')
xlabel('Polynomial Order');
set(gca,'Fontsize',60);
saveas(fig1, 'fhn.jpg');
%%
m = find(NRMSE==min(NRMSE)); %polinomun derecesi

% m = 5; %polinomun derecesi
% x = v_sp2(1,:);
% u = w_sp2(1,:);

t = 1:1:length(x);

a = polyfit(t,x,m);
b = polyfit(t,u,m);

v_sym_t = poly2sym(a,tt); % zaman t ye göre v yi yeniden oluþturur.
w_sym_t = poly2sym(b,tt);

for k = 1:length(t)

    fhn_cal_v(k) = polyval(a,t(k));
    fhn_cal_w(k) = polyval(b,t(k));

end

m_s_e_v = immse(x,fhn_cal_v)
m_a_e_v = mae(x,fhn_cal_v)
r_m_s_e_v = sqrt(m_s_e_v)
n_r_m_s_e_v = (r_m_s_e_v/(max(x) - min(x)))*100
corr = corr2(x,fhn_cal_v);

m_s_e_w = immse(u,fhn_cal_w)
m_a_e_w = mae(u,fhn_cal_w)
r_m_s_e_w = sqrt(m_s_e_w)
n_r_m_s_e_w = (r_m_s_e_w/(max(u) - min(u)))*100
corr = corr2(u,fhn_cal_w);

fig1 = figure('Position',get(0,'Screensize'));
plot(x,'o','Color','r','linewidth',10);
grid on; hold on;
plot(fhn_cal_v,'-.','Color','k','linewidth',10);
ylabel('v, [V]')
xlabel('Number of Samples');
set(gca,'Fontsize',60);
legend({'Original','Curve Fitting Polynomial Method'},'Location','southwest');
saveas(fig1, 'fhn_sp7_v.jpg');

fig2 = figure('Position',get(0,'Screensize'));
plot(u,'o','Color','r','linewidth',10);
grid on; hold on;
plot(fhn_cal_w,'-.','Color','k','linewidth',10);
ylabel('u, [V]')
xlabel('Number of Samples');
set(gca,'Fontsize',60);
legend({'Original','Curve Fitting Polynomial Method'},'Location','southwest');
saveas(fig2, 'fhn_sp7_w.jpg');

fig3 = figure('Position',get(0,'Screensize'));
plot(x,u,'LineStyle','-','Marker','none','Color','r','MarkerSize',25,'linewidth',4)
grid on; hold on;
plot(fhn_cal_v,fhn_cal_w,'LineStyle','-','Marker','none','Color','k','MarkerSize',25,'linewidth',4)
ylabel('u')
xlabel('v');
set(gca,'Fontsize',60);
legend({'Original','Curve Fitting Polynomial Method'},'Location','southwest');
saveas(fig3, 'fhn_sp7_vw.jpg');