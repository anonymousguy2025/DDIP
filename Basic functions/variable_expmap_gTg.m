function [g,Tg]=variable_expmap_gTg(Gamma)

k      = Gamma(1:3);
theta  = norm(k);

Gammahat  = dinamico_hat(Gamma);
adjGamma  = dinamico_adj(Gamma);

Gammahatp2 = Gammahat*Gammahat;
Gammahatp3 = Gammahatp2*Gammahat;

adjGammap2 = adjGamma*adjGamma;
adjGammap3 = adjGammap2*adjGamma;
adjGammap4 = adjGammap3*adjGamma;

if (theta<=1e-2)
    g  = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]+Gammahat+Gammahatp2/2+Gammahatp3/6;
    
    f1 = 1/2;
    f2 = 1/6;
    f3 = 1/24;
    f4 = 1/120;
    
    Tg  = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]+f1*adjGamma+f2*adjGammap2+f3*adjGammap3+f4*adjGammap4;
else
    
    tp2        = theta*theta;
    tp3        = tp2*theta;
    tp4        = tp3*theta;
    tp5        = tp4*theta;
    
    sintheta   = sin(theta);
    costheta   = cos(theta);
    
    t1 = theta*sintheta;
    t2 = theta*costheta;
    
    g   = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]+Gammahat+...
          (1-costheta)/(tp2)*Gammahatp2+...
          ((theta-sintheta)/(tp3))*Gammahatp3;
    Tg  = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]+((4-4*costheta-t1)/(2*tp2))*adjGamma+...
          ((4*theta-5*sintheta+t2)/(2*tp3))*adjGammap2+...
          ((2-2*costheta-t1)/(2*tp4))*adjGammap3+...
          ((2*theta-3*sintheta+t2)/(2*tp5))*adjGammap4;
end
