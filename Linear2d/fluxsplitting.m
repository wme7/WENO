function [vp,vn] = fluxsplitting(u,f,df,strategy)
% split the flux into right-going and left-going
% OUTPUT:
%   * vp: positive flux, v^{+}, which corresponds to f_{i+1/2}^{-}
%   * vn: negative flux  v^{-}, which corresponds to f_{i+1/2}^{+}

switch strategy
    case 'Upwind' % Godunov/scalar fluxsplit (non-conservative)
        vp = f((u + abs(u))./2); %flux^{+}
        vn = f((u - abs(u))./2); %flux^{-}
    case 'LLF' % Local Lax-Friedrichs
        v = f(u); alpha = abs(df(u));
        vp = 0.5.*(v + alpha.*u); %flux^{+}
        vn = 0.5.*(v - alpha.*u); %flux^{-}
    case 'LF' % (Global) Lax-Friedrichs
        v = f(u); alpha = max(abs(df(u)));
        vp = 0.5.*(v + alpha.*u); %flux^{+}
        vn = 0.5.*(v - alpha.*u); %flux^{-}
    otherwise
        error('only case not available')
end