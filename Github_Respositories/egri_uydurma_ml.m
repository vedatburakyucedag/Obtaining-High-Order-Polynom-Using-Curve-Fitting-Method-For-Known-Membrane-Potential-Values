clc; clear;

% fhn osilator egri uydurma polinomal yontemi

syms tt

%step size 
hh = 0.1; 
N = 700-1;
t(1) = 0;

C_sp1 = 8; E_Ca_sp1 = 20;  E_K_sp1 = - 60; E_Leak_sp1 = - 70; g_Ca_sp1 = 20;  g_K_sp1 = 20; g_Leak_sp1 = 0.5; T0_sp1 = 13;  Qw_sp1 = T0_sp1; V1_sp1 = -12;  V2_sp1 = 18; V3_sp1 = -10;   V4_sp1 = 13;   I_sp1 = 50;
V_sp1(1) = -12.7743548965288; W_sp1(1) = 0.376576642801617;

C_sp2 = 10; E_Ca_sp2 = 20;  E_K_sp2 = - 70; E_Leak_sp2 = - 70; g_Ca_sp2 = 8;  g_K_sp2 = 35; g_Leak_sp2 = 0.1; T0_sp2 = 13;  Qw_sp2 = T0_sp2; V1_sp2 = -12;  V2_sp2 = 13; V3_sp2 = -10;   V4_sp2 = 13;   I_sp2 = 35;
V_sp2(1) = -23.0105113265647; W_sp2(1) = 0.0330094530569983;

[V_sp1,W_sp1,t] = function_ml_1n(E_Ca_sp1,E_K_sp1,E_Leak_sp1,g_Ca_sp1,g_K_sp1,g_Leak_sp1,V1_sp1,V2_sp1,V3_sp1,V4_sp1,T0_sp1,I_sp1,hh,N,t,V_sp1,W_sp1,C_sp1);
[V_sp2,W_sp2,t] = function_ml_1n(E_Ca_sp2,E_K_sp2,E_Leak_sp2,g_Ca_sp2,g_K_sp2,g_Leak_sp2,V1_sp2,V2_sp2,V3_sp2,V4_sp2,T0_sp2,I_sp2,hh,N,t,V_sp2,W_sp2,C_sp2);

figure(1); clf(1);
plot(V_sp1,'k','LineWidth',2)
xlabel('Number of Samples')
ylabel('Membran Potential (mV)')
set(gca,'Fontsize',30)
grid on

x = V_sp2;
u = W_sp2;

% ML_known = load('mat_ml_rk_sp123456.mat');
% x = ML_known.V_sp1;
% u = ML_known.W_sp1;

for i=1:1:150  
    
    m = i; %polinomun derecesi
        
    t = 1:1:length(x);
    
    a = polyfit(t,x,m);
    b = polyfit(t,u,m);
    
    v_sym_t = poly2sym(a,tt); % zaman t ye göre v yi yeniden oluþturur.
    w_sym_t = poly2sym(b,tt);
    
    for k = 1:length(t)
    
        ML_cal_v(k) = polyval(a,t(k));
        ML_cal_w(k) = polyval(b,t(k));
    
    end
    
    m_s_e_v = immse(x,ML_cal_v);
    m_a_e_v = mae(x,ML_cal_v);
    r_m_s_e_v = sqrt(m_s_e_v);
    n_r_m_s_e_v = (r_m_s_e_v/(max(x) - min(x)))*100;
    corr = corr2(x,ML_cal_v);
    
    m_s_e_w = immse(u,ML_cal_w);
    m_a_e_w = mae(u,ML_cal_w);
    r_m_s_e_w = sqrt(m_s_e_w);
    n_r_m_s_e_w = (r_m_s_e_w/(max(u) - min(u)))*100;
    corr = corr2(u,ML_cal_w);

    NRMSE(i) = n_r_m_s_e_v;

end

[Min_vlue,Min_locs] = min(NRMSE);

fig1 = figure('Position',get(0,'Screensize'));
plot(NRMSE(1:100),'LineStyle','-','Marker','o','Color','k','MarkerSize',20,'linewidth',4)
grid on; hold on;
plot(Min_locs,Min_vlue,'Marker','*','Color','r','MarkerSize',30,'linewidth',6)
xline(Min_locs, 'Color', 'r', 'LineWidth', 4)
textLabel = sprintf('Min of %.4f NRMSE at Polynominal Order=%.0f', Min_vlue, Min_locs);
text(0, 55, textLabel, 'fontSize', 50, 'Color', 'r', 'VerticalAlignment','middle')
ylabel('NRMSE')
xlabel('Polynomial Order');
set(gca,'Fontsize',60);
saveas(fig1, 'ml.jpg');

%%
m = find(NRMSE==min(NRMSE)); %polinomun derecesi

% m = 5; %polinomun derecesi
% 
% x = V_sp1;
% u = W_sp1;

t = 1:1:length(x);

a = polyfit(t,x,m);
b = polyfit(t,u,m);

v_sym_t = poly2sym(a,tt); % zaman t ye göre v yi yeniden oluþturur.
w_sym_t = poly2sym(b,tt);

for k = 1:length(t)

    ML_cal_v(k) = polyval(a,t(k));
    ML_cal_w(k) = polyval(b,t(k));

end

m_a_e_v = mae(x,ML_cal_v)
m_s_e_v = immse(x,ML_cal_v)
r_m_s_e_v = sqrt(m_s_e_v)
n_r_m_s_e_v = (r_m_s_e_v/(max(x) - min(x)))*100
corr = corr2(x,ML_cal_v);

m_a_e_w = mae(u,ML_cal_w)
m_s_e_w = immse(u,ML_cal_w)
r_m_s_e_w = sqrt(m_s_e_w)
n_r_m_s_e_w = (r_m_s_e_w/(max(u) - min(u)))*100
corr = corr2(u,ML_cal_w);

fig1 = figure('Position',get(0,'Screensize'));
plot(x,'o','Color','r','linewidth',10);
grid on; hold on;
plot(ML_cal_v,'-.','Color','k','linewidth',10);
ylabel('V, [mV]')
xlabel('Number of Samples');
set(gca,'Fontsize',60);
legend({'Original','Curve Fitting Polynomial Method'},'Location','southwest');
saveas(fig1, 'fhn_sp7_v.jpg');

fig2 = figure('Position',get(0,'Screensize'));
plot(u,'o','Color','r','linewidth',10);
grid on; hold on;
plot(ML_cal_w,'-.','Color','k','linewidth',10);
ylabel('W, [mV]')
xlabel('Number of Samples');
set(gca,'Fontsize',60);
legend({'Original','Curve Fitting Polynomial Method'},'Location','southwest');
saveas(fig2, 'fhn_sp7_w.jpg');

fig3 = figure('Position',get(0,'Screensize'));
plot(x,u,'LineStyle','-','Marker','none','Color','r','MarkerSize',25,'linewidth',4)
grid on; hold on;
plot(ML_cal_v,ML_cal_w,'LineStyle','-','Marker','none','Color','k','MarkerSize',25,'linewidth',4)
ylabel('W')
xlabel('V');
set(gca,'Fontsize',60);
legend({'Original','Curve Fitting Polynomial Method'},'Location','southwest');
saveas(fig3, 'fhn_sp7_vw.jpg');