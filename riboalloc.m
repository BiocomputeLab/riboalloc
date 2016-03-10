function out = riboalloc (R_pool, a_e, a_h, N_e, N_h, t_e, t_h, s_h, t_end)
%--------------------------------------------------------------------------
% Model of the ribosome usage between two different types of transcripts:
% (i) an endogenous "average" transcript, and (ii) a heterologous foreign
% transcript. The model consists of the following parameters, provided as
% arguments to this function:
%
%   R_pool:   number of total ribosomes in the cell.
%
%   a_e/a_h:  translation initiation rates for endogenous and heterologous 
%             transcripts.
%
%   N_e/N_h:  number of endogenous and heterologous transcripts.
%   
%   t_e/t_h:  translation time of the endogenous and heterologous 
%             transcripts.
%
%   s_h:      ribosome sites per heterlogous transcript.
%
%   t_end:    length of time to simulate for.
%
%--------------------------------------------------------------------------
% Author:  Thomas E. Gorochowski
% Updated: 19th August 2015
% License: GNU GPLv3 license
%--------------------------------------------------------------------------

    % Delays for the [endogenous, heterologous] number of transcripts
    lags = [t_e, t_h];
    max_t = max(lags);
    
    % Solve the equation over the interval [0, t_end]
    options = ddeset('AbsTol',1e-6);
    out = dde23(@ddeRiboAlloc, lags, @ddeRiboAllocHis, [0, t_end], options);
    
    % The DDE equation to solve
    function dydt = ddeRiboAlloc(t,y,Z)
        % Grab the historical states
        ylag_e = Z(:,1);
        ylag_h = Z(:,2);
        
        % Variable to store derivative
        dydt = zeros(2, 1);
        
        % Ribosome allocation (now)
        R_e = y(1);
        R_h = y(2);
        R_f = R_pool - R_e - R_h;
        % Ribosome allocation (delayed for t_e)
        R_e_de = ylag_e(1);
        R_h_de = ylag_e(2);
        R_f_de = R_pool - R_e_de - R_h_de;
        % Ribosome allocation (delayed for t_h)
        R_e_dh = ylag_h(1);
        R_h_dh = ylag_h(2);
        R_f_dh = R_pool - R_e_dh - R_h_dh;
        
        % Calculate the probability of collision (proportional to density)   
        Pc_h = R_h/(N_h*s_h);
        if Pc_h < 0; Pc_h = 0; end
        if Pc_h > 1; Pc_h = 1; end
        Pc_h_dh = R_h_dh/(N_h*s_h);
        if Pc_h_dh < 0; Pc_h_dh = 0; end
        if Pc_h_dh > 1; Pc_h_dh = 1; end
        
        % Rates on mRNA taking collisions into account
        Rc_h = (s_h/t_h)*(1-Pc_h);
        Rc_h_dh = (s_h/t_h)*(1-Pc_h_dh);
        
        % Effective initation rates (cannot exceed rate on mRNA)
        eff_a_h = a_h ;
        if eff_a_h > Rc_h; eff_a_h = Rc_h; end
        eff_a_h_dh = a_h;
        if eff_a_h_dh > Rc_h_dh; eff_a_h_dh = Rc_h_dh; end
        
        % Handle the history at the start of the simulation
        if t < t_e
            dydt(1) = (R_f * a_e * N_e);
        else
            dydt(1) = (R_f * a_e * N_e) - (R_f_de * a_e * N_e);
        end
        if t < t_h
            dydt(2) = (R_f * eff_a_h * N_h);
        else
            dydt(2) = (R_f * eff_a_h * N_h) - (R_f_dh * eff_a_h_dh * N_h);
        end
    end

    % Provide initial historical conditions
    function s = ddeRiboAllocHis(t)
        s = zeros(2, 1);
    end
end
