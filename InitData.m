%% Initialize Data
function [input, data] = InitData(settings, ic, cost, ref)

    if nargin<3
        externalICandRef = false;
    else
        externalICandRef = true;
    end

    nx = settings.nx;       % No. of differential states
    nu = settings.nu;       % No. of controls
    nz = settings.nz;       % No. of algebraic states
    ny = settings.ny;        % No. of outputs (references)    
    nyN= settings.nyN;       % No. of outputs at terminal stage 
    np = settings.np;        % No. of parameters (on-line data)
    nc = settings.nc;        % No. of constraints
    ncN = settings.ncN;      % No. of constraints at terminal stage
    N  = settings.N;         % No. of shooting points
    nbx = settings.nbx;      % No. of state bounds
    nbu = settings.nbu;      % No. of control bounds
    nbu_idx = settings.nbu_idx;  % Index of control bounds

    switch settings.model
                      
        case 'InvertedPendulum'
            input.x0 = [0;pi;0;0];    
            input.u0 = zeros(nu,1); 
            input.z0 = zeros(nz,1);
            para0 = 0;  

            Q=repmat([10 10 0.1 0.1 0.01]',1,N);
            QN=[10 10 0.1 0.1]';

            % upper and lower bounds for states (=nbx)
            lb_x = -2;
            ub_x = 2;

            % upper and lower bounds for controls (=nbu)           
            lb_u = -20;
            ub_u = 20;
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];

        case 'ChainofMasses_Lin'
            n=5;
            data.n=n;
            input.x0=zeros(nx,1);
            for i=1:n
                input.x0(i)=7.5*i/n;
            end
            input.u0=zeros(nu,1);
            input.z0 = zeros(nz,1);
            para0=0;
            wv=[];wx=[];
            wu = [0.1 0.1 0.1];
            for i=1:3
                wx = [wx, 25];
                wv = [wv, 0.25*ones(1,n-1)];
            end
            Q = repmat([wx,wv,wu]',1,N);
            QN= [wx,wv]';

            % upper and lower bounds for states (=nbx)
            lb_x = [];
            ub_x = [];

            % upper and lower bounds for controls (=nbu)           
            lb_u = [-1;-1;-1];
            ub_u = [1;1;1];
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];

        case 'ChainofMasses_NLin'
            n=10;
            data.n=n;
            input.x0=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 zeros(1,nx-n)]';
            input.u0=zeros(nu,1);
            input.z0 = zeros(nz,1);
            para0=0;
            wv=[];wx=[];
            wu=[0.01, 0.01, 0.01];
            for i=1:3
                wx=[wx,25];
                wv=[wv,ones(1,n-1)];
            end
            Q = repmat([wx,wv,wu]',1,N);
            QN= [wx,wv]';

            % upper and lower bounds for states (=nbx)
            lb_x = [];
            ub_x = [];

            % upper and lower bounds for controls (=nbu)           
            lb_u = [-1;-1;-1];
            ub_u = [1;1;1];
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];
                                       
        case 'TethUAV'
            
            input.x0=[0; 0; 0; 0; 9.81; 0];%zeros(nx,1);
            input.u0=[0; 0; 0; 0];%zeros(nu,1);
            input.z0 = zeros(nz,1);
            alpha = pi/6;
            para0=[-alpha; alpha];
            % phi phi_dot theta theta_dot 
            q = [200, 1, 200, 0, 0.0001, 0.0001, 1, 1, 1, 0, 5000, 5000];

            qN = q(1:nyN);
            Q = repmat(q',1,N);
            QN = qN';
            
            fR_min = 0;%-inf;
            fR_max = 15;%inf;
            tauR_min = -1.2;%-inf;
            tauR_max = 1.2;%inf;
            fL_min = 0;%-inf;
            fL_max = 10;%inf;
            constr_max = 0;
            constr_min = -inf;
            s1_min = 0;
            s1_max = inf;
            s2_min = 0;
            s2_max = inf;
            
            % upper and lower bounds for states (=nbx) if f1,2 are f_R, tau_R
            lb_x = [fR_min; tauR_min];%0*ones(nbx,1);
            ub_x = [fR_max; tauR_max]; %omegaMax*ones(nbx,1);
            
            % upper and lower bounds for controls (=nbu)           
            lb_u = [s1_min; s2_min];
            ub_u = [s1_max; s2_max];
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [fL_min; constr_min; constr_min];
            ub_g = [fL_max; constr_max; constr_max];            
            lb_gN = [fL_min];
            ub_gN = [fL_max]; 
            
        case 'DiM'	
            input.x0 = zeros(nx,1);    % initial state	
            input.u0 = zeros(nu,1);    % initial control
            input.z0 = zeros(nz,1);
            para0 = 0;  % initial parameters (by default a np by 1 vector, if there is no parameter, set para0=0)	

             %weighting matrices	
            Q=[1200,1200,2000,800,800,5800,... % perceived acc and angular vel	
                    32000*1.1,32000*1.1,1600*1,... %px,py,pz hex	
                    3200*1.1,3200*1.1,2000*1,... %vx, vy, vz hex	
                    4600*1,600*1,... % x,y tri	
                    850*1,850*1,... % vx,vy tri	
                    3700,3000,1500,... % phi, theta, psi hex	
                    750,... % phi tri	
                    0.01,0.0,0.0,... % omega phi,theta,psi hex	
                    500.0,... % omega phi tri	
                    0.0,0.0,0.001,... %ax,ay,az hex %         20*1.1,20*1.1,... % ax,ay tri	
                    0.0,0.01,0.1 ... % alpha phi,theta, psi hex 	
                    ];	
              Q = repmat(Q',1,N);	

              QN=Q(1:nyN,1);	

               % upper and lower bounds for states (=nbx)	
              lb_x = [];	
              ub_x = [];	

               % upper and lower bounds for controls (=nbu)           	
              lb_u = [];	
              ub_u = [];	

               % upper and lower bounds for general constraints (=nc)	
              lb_g=[1.045;1.045;1.045;1.045;1.045;1.045];    % lower bounds for ineq constraints	
              ub_g=[1.3750;1.3750;1.3750;1.3750;1.3750;1.3750];  % upper bounds for ineq constraints	
              lb_gN=[1.045;1.045;1.045;1.045;1.045;1.045];  % lower bounds for ineq constraints at terminal point	
              ub_gN=[1.3750;1.3750;1.3750;1.3750;1.3750;1.3750];  % upper bounds for ineq constraints at terminal point
              
        case 'TurboEngine'
            input.x0 = [1.2; 1.2; 0; 0];
            input.u0 = zeros(nu,1); 
            input.z0 = [1.1; 1.1];
            para0 = [2000; -0.3];  

            Q=repmat([10 1e-7*0.05 1e-6*0.05]',1,N);
            QN=[10]';

            % upper and lower bounds for states (=nbx)
            lb_x = [0;0];
            ub_x = [100;100];

            % upper and lower bounds for controls (=nbu)           
            lb_u = [-800; -800];
            ub_u = [800; 800];
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [0; 0; 0];
            ub_g = [2; 90e3/60; 180e3/60];        
            lb_gN = [0; 0; 0];
            ub_gN = [2; 90e3/60; 180e3/60];

        case 'highFidelityTiltingQuadrotor'
            if( externalICandRef==true )
                % initial conditions
                input.x0    = ic.x0;
                input.u0    = ic.u0;
                input.z0    = ic.z0;
                para0       = ic.para0;
                % cost weight matrices
                Q           = repmat(cost.Q(:), 1, N);
                QN          = cost.QN(:);
                % constraints - state
                lb_x = cost.lb_x;
                ub_x = cost.ub_x;
                % constraints - input
                lb_u = cost.lb_u;
                ub_u = cost.ub_u;
            else
                input.x0 = [
                    zeros(3,1)
                    [1; 0; 0; 0]
                    zeros(3,1)
                    zeros(3,1)
                    zeros(4,1)
                    .33*ones(4,1) % hardcoded hovering normalized spinning rate square
                    ];
                input.u0 = [
                    zeros(4,1)
                    zeros(4,1)
                    ];
                input.z0 = zeros(nz,1);
                para0 = [1; 0; 0; 0];

                Qp = 5*[1 1 1];                 % position weight
                Qeq = [1 1 1];                  % quaternion weight
                Qv = 0.5*[1 1 1];               % linear speed weight
                Qomega = 10*[1 1 1];            % angular speed weight
                Qalpha = 1e-2*[1 1 1 1];        % alpha angle weight
                Qw_bar2 = [1 1 1 1]*1e-1;       % propellers spinning rate square weight
                Qup = [1 1 1 1]*1e-1;           % propellers spinning acceleration square weight
                Qdotalpha = 1e-3*[1 1 1 1];     % alpha dot weight
                Q=repmat([Qp Qeq Qv Qomega Qalpha Qw_bar2 Qup Qdotalpha]',1, N);
                QN=[Qp Qeq Qv Qomega Qalpha Qw_bar2]';

                % upper and lower bounds for states (=nbx)
                lb_x = [
                    -pi/3 .* ones(4,1)  % alpha min: -30 deg
                    0.01.*ones(4,1)     % w_bar2 min: 0.01 of max
                    ];
                ub_x = [
                    pi/3 .* ones(4,1)   % alpha max: 30 deg
                    ones(4,1)           % w_bar2 max: 100% of max
                    ];

                % upper and lower bounds for controls (=nbu)
                lb_u = [
                    -10*ones(4,1)       % w_bar2_dot_min: -10 s^{-1}
                    -6*ones(4,1)        % alpha_dot_min: -6 rad/s
                    ]; % MIN SPEED AT 1%
                ub_u = [
                    10*ones(4,1)        % w_bar2_dot_min: 10 s^{-1}
                    6*ones(4,1)         % alpha_dot_max: 6 rad/s
                    ];
            end
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];
            lb_gN = [];
            ub_gN = [];

    end

    % prepare the data
    
    input.lb = repmat(lb_g,N,1);
    input.ub = repmat(ub_g,N,1);
    input.lb = [input.lb;lb_gN];
    input.ub = [input.ub;ub_gN];
            
    lbu = -inf(nu,1);
    ubu = inf(nu,1);
    for i=1:nbu
        lbu(nbu_idx(i)) = lb_u(i);
        ubu(nbu_idx(i)) = ub_u(i);
    end
                
    input.lbu = repmat(lbu,1,N);
    input.ubu = repmat(ubu,1,N);
    
    input.lbx = repmat(lb_x,1,N);
    input.ubx = repmat(ub_x,1,N);
        
    x = repmat(input.x0,1,N+1);  % initialize all shooting points with the same initial state
    u = repmat(input.u0,1,N);    % initialize all controls with the same initial control
    z = repmat(input.z0,1,N);    % initialize all algebraic state with the same initial condition
    para = repmat(para0,1,N+1);  % initialize all parameters with the same initial para
         
    input.x=x;           % (nx by N+1)
    input.u=u;           % (nu by N)
    input.z=z;           % (nz by N)
    input.od=para;       % (np by N+1)
    input.W=Q;           % (ny by N)
    input.WN=QN;         % (nyN by 1)
     
    input.lambda=zeros(nx,N+1);   % langrangian multiplier w.r.t. equality constraints
    input.mu=zeros(N*nc+ncN,1);   % langrangian multipliers w.r.t. general inequality constraints
    input.mu_u = zeros(N*nu,1);   % langrangian multipliers w.r.t. input bounds
    input.mu_x = zeros(N*nbx,1);  % langrangian multipliers w.r.t. state bounds
    
    %% Reference generation

    switch settings.model

        case 'InvertedPendulum'

            data.REF=zeros(1,nx+nu);

        case 'ChainofMasses_Lin'

            data.REF=[7.5,0,0,zeros(1,3*(n-1)),zeros(1,nu)];

        case 'ChainofMasses_NLin'

            data.REF=[1,0,0,zeros(1,3*(n-1)),zeros(1,nu)];
    
        case 'TethUAV'
            
        	data.REF = zeros(1, ny);
            
        case 'DiM'

             load REF_DiM_2;

             REF_DiM_2 = [REF_DiM_2, zeros(5000,24)];

             data.REF = REF_DiM_2;
             
        case 'TurboEngine'
            
            data.REF=[1.4, 0, 0];

        case 'highFidelityTiltingQuadrotor'
            if( externalICandRef==false )
                p_ref = [1 1 1];
                q_e_ref = [0 0 0];
                v_ref = [0 0 0];
                omega_ref = [0 0 0];
                alpha_ref = [0 0 0 0];
                w_bar2 = .33*ones(1,4); % HARDCODED FOR HOVERING
                u_ref = [0 0 0 0];
                dot_alpha_ref = [0 0 0 0];
                data.REF=[p_ref q_e_ref v_ref omega_ref alpha_ref w_bar2 u_ref dot_alpha_ref];
            else
                data.REF = ref.h;
                data.param = ref.param;
            end

    end
    
end