function [vp,vm]=fluxsplitting_scalar(u,strategy)
% split the flux into right-going and left-going
% OUTPUT:
%   * vp: positive dflux, v^{+}, which corresponds to f_{i+1/2}^{-}
%   * vn: negative dflux  v^{-}, which corresponds to f_{i+1/2}^{+}

switch strategy
    case 'Upwind'
        vp= max(u,0); %dflux^{+}
        vm=-min(u,0); %dflux^{-}
    case 'LF'
        au=max(abs(u(:)));
        vp =  0.5*(u+au); %dflux^{+}
        vm = -0.5*(u-au); %dflux^{-}
    otherwise
        error('only case not available')
end