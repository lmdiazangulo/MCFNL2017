 function [ e ] = computeFluxes(e, alpha)
% function [ e ] = computeFluxes(e)
% Assembles fluxes.

% ------ Preliminar definitions.
K = numel(e);
Np = size(e(1).E, 1);
% ------ Computes fluxes ---------------
% First element.
switch alpha
    case 0
        ZAvR = (e(1).Z + e(2).Z)/2;      % ZAvR := Average impedance on right side.
        YAvR = 1 / ZAvR;                 % YAvR := Average admitance on right side.
        e(1).fluxE(1) = (1/e(1).Z).*( -e(1).Z * e(1).dH(1) - e(1).dE(1) );
        e(1).fluxE(Np)= (1/ ZAvR ).*( +e(2).Z * e(1).dH(2) - e(1).dE(2) );
        e(1).fluxH(1) = (1/e(1).Y).*( -e(1).Y * e(1).dE(1) - e(1).dH(1) );
        e(1).fluxH(Np)= (1/ YAvR ).*( +e(2).Y * e(1).dE(2) - e(1).dH(2) );
        % Bulk of elements.
        for k = 2:(K-1)
            % ----------- NOTE: This can be computed out for efficiency -----------
            ZAvL = (e(k-1).Z + e(k).Z  )/2;
            ZAvR = (e(k).Z   + e(k+1).Z)/2;
            YAvL = 1 / ZAvL;
            YAvR = 1 / ZAvR;
            % ---------------------------------------------------------------------
            e(k).fluxE(1)  = (1/ZAvL).*( -e(k-1).Z*e(k).dH(1)-e(k).dE(1) );
            e(k).fluxE(Np) = (1/ZAvR).*( +e(k+1).Z*e(k).dH(2)-e(k).dE(2) );            
            e(k).fluxH(1)  = (1/YAvL).*( -e(k-1).Y*e(k).dE(1)-e(k).dH(1) );
            e(k).fluxH(Np) = (1/YAvR).*( +e(k+1).Y*e(k).dE(2)-e(k).dH(2) );
        end
        % Last element.
        ZAvL = (e(K-1).Z + e(K).Z)/2;    % ZAvL := Average impedance on left side.
        YAvL = 1 / ZAvL;                 % YAvL := Average admitance on left side.
        e(K).fluxE(1) = (1/ ZAvL ).*( -e(K-1).Z * e(K).dH(1) - e(K).dE(1) );
        e(K).fluxE(Np)= (1/e(K).Z).*( +e(K).Z   * e(K).dH(2) - e(K).dE(2) );
        e(K).fluxH(1) = (1/ YAvL ).*( -e(K-1).Y * e(K).dE(1) - e(K).dH(1) );
        e(K).fluxH(Np)= (1/e(K).Y).*( +e(K).Y   * e(K).dE(2) - e(K).dH(2) );
    case 1
        e(1).fluxE(1) =  - e(1).dH(1) - (1/e(1).Z).*e(1).dE(1);
        e(1).fluxH(1) =  - e(1).dE(1) - (1/e(1).Y).*e(1).dH(1);
        e(1).fluxE(Np)= e(1).dH(2);
        e(1).fluxH(Np)= e(1).dE(2);
        % Bulk of elements.
        for k = 2:(K-1)
            % ---------------------------------------------------------------------
            e(k).fluxE(1)  = -e(k).dH(1);
            e(k).fluxE(Np) = e(k).dH(2);
            e(k).fluxH(1)  = -e(k).dE(1);
            e(k).fluxH(Np) = e(k).dE(2);
        end
        % Last element.
        e(K).fluxE(1) = - e(K).dH(1);
        e(K).fluxH(1) = - e(K).dE(1);
        e(K).fluxE(Np)=   e(K).dH(2); 
        e(K).fluxH(Np)=   e(K).dE(2); 
%         e(K).fluxE(Np)=  e(K).dH(2) - e(K).Y .* e(K).dE(2) ;        
%         e(K).fluxH(Np)=  e(K).dE(2) - e(K).Z .* e(K).dH(2) ;
    otherwise
        ZAvR = (e(1).Z + e(2).Z)/2;      % ZAvR := Average impedance on right side.
        YAvR = 1 / ZAvR;                 % YAvR := Average admitance on right side.
        e(1).fluxE(1) = (1/e(1).Z).*( -e(1).Z * e(1).dH(1) ...
                         - (1-alpha)*e(1).dE(1) );
        e(1).fluxE(Np)= (1/ ZAvR ).*( +e(2).Z * e(1).dH(2) ...
                         - e(1).dE(2) );
        e(1).fluxH(1) = (1/e(1).Y).*( -e(1).Y * e(1).dE(1) ...
                         - e(1).dH(1) );
        e(1).fluxH(Np)= (1/ YAvR ).*( +e(2).Y * e(1).dE(2) ...
                         - e(1).dH(2) );
        % Bulk of elements.
        for k = 2:(K-1)
            % ----------- NOTE: This can be computed out for efficiency -----------
            ZAvL = (e(k-1).Z + e(k).Z  )/2;
            ZAvR = (e(k).Z   + e(k+1).Z)/2;
            YAvL = 1 / ZAvL;
            YAvR = 1 / ZAvR;
            % ---------------------------------------------------------------------
            e(k).fluxE(1)  = (1/ZAvL).*( -e(k-1).Z*e(k).dH(1) ... 
                             - e(k).dE(1) );
            e(k).fluxE(Np) = (1/ZAvR).*( +e(k+1).Z*e(k).dH(2) ...
                             - e(k).dE(2) );
            e(k).fluxH(1)  = (1/YAvL).*( -e(k-1).Y*e(k).dE(1) ...
                             - e(k).dH(1) );
            e(k).fluxH(Np) = (1/YAvR).*( +e(k+1).Y*e(k).dE(2) ...
                             - e(k).dH(2) );
        end
        % Last element.
        ZAvL = (e(K-1).Z + e(K).Z)/2;    % ZAvL := Average impedance on left side.
        YAvL = 1 / ZAvL;                 % YAvL := Average admitance on left side.
        e(K).fluxE(1) = (1/ ZAvL ).*( -e(K-1).Z * e(K).dH(1) ...
                        - e(K).dE(1) );
        e(K).fluxE(Np)= (1/e(K).Z).*( +e(K).Z   * e(K).dH(2) ...
                        - e(K).dE(2) );
        e(K).fluxH(1) = (1/ YAvL ).*( -e(K-1).Y * e(K).dE(1) ...
                        - e(K).dH(1) );
        e(K).fluxH(Np)= (1/e(K).Y).*( +e(K).Y   * e(K).dE(2) ...
                        - e(K).dH(2) );
end