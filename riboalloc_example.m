%--------------------------------------------------------------------------
% Author:  Thomas E. Gorochowski
% Updated: 19th August 2015
% License: GNU GPL-3.0 license
%--------------------------------------------------------------------------

% Model parameters
R_pool = 26300.0;    % http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=102015&ver=9&trm=ribosomes%20e.%20coli
a_e    = 0.5/R_pool; % ~0.5 inits/sec (5 to 58 init a min) based on R_pool
a_h    = 0.5/R_pool; % ~0.5 inits/sec (5 to 58 init a min) based on R_pool
N_e    = 1380.0*3;   % http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=100064&ver=9&trm=number%20of%20transcripts%20e.%20coli
                     % ~3 genes per operon http://gbe.oxfordjournals.org/content/5/11/2242.full
N_h    = 100.0;      % Steady-state mRNA count (this is varied)
t_e    = 17.0;       % Average 340aa gene ~ 20 aa/sec
t_h    = 17.0;       % Average 340aa gene ~ 20 aa/sec
s_h    = 340/9;      % Average 340 codons, ribosome ~9 codons footprint
t_end  = 10*60;      % Time to run simulation

% Run the model
%sol = riboalloc(R_pool, a_e, a_h*1000, N_e, 1000.0, t_e, cur_t_h, s_h, t_end);
sol = riboalloc(R_pool, a_e, a_h, N_e, 100.0, t_e, cur_t_h, s_h, t_end);

% Plot the solution
figure;
plot(sol.x,sol.y);
title('RiboAlloc DDE Model');
xlabel('time t');
ylabel('solution y');
legend('R_{e}','R_{h}','Location', 'NorthWest');
