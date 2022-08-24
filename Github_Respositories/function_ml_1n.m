function [V,W,t] = function_ml_1n(E_Ca,E_K,E_Leak,g_Ca,g_K,g_Leak,V1,V2,V3,V4,T0,I,hh,N,t,V,W,C)

% define function handles
fV=@(t,V,W)  ((I - ((g_Ca * 0.5 * (1 + tanh((V - V1) / V2)) * (V - E_Ca)) + (g_K * W * (V - E_K)) + (g_Leak * (V - E_Leak)))) / C) ;
fW=@(t,V,W)  ((0.5 * (1 + tanh((V - V3)/V4)) - W) / ((T0) / (cosh((V - V3) / (2 * V4)))));

    for i=1:N

        t(i+1) = t(i) + hh;

        k1V = fV(t(i)      ,V(i)           ,W(i)  );
        k1W = fW(t(i)      ,V(i)           ,W(i)  );

        k2V = fV(t(i)+hh/2 , V(i)+hh/2*k1V , W(i)+hh/2*k1W);
        k2W = fW(t(i)+hh/2 , V(i)+hh/2*k1V , W(i)+hh/2*k1W);

        k3V = fV(t(i)+hh/2 , V(i)+hh/2*k2V , W(i)+hh/2*k2W);
        k3W = fW(t(i)+hh/2 , V(i)+hh/2*k2V , W(i)+hh/2*k2W);

        k4V = fV(t(i)+hh   , V(i)+hh  *k3V , W(i)+hh  *k3W);
        k4W = fW(t(i)+hh   , V(i)+hh  *k3V , W(i)+hh  *k3W);

        V(i+1) = V(i) + hh / 6 * (k1V + 2*k2V + 2*k3V + k4V);
        W(i+1) = W(i) + hh / 6 * (k1W + 2*k2W + 2*k3W + k4W);

    end
        
end

