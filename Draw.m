%% Plot results

% define colors
blue=[0 0 255]/255;
red=[220 20 60]/255;
orange=[255 165 0]/255;
green=[0 205 102]/255;

% start generating pictures
switch settings.model
    
    case 'InvertedPendulum'
        
        figure(1);
        subplot(321)
        plot(time,state_sim(:,1));
        title('p');
        subplot(322)
        plot(time,state_sim(:,2)*180/pi);
        title('\theta');
        subplot(323)
        plot(time,state_sim(:,3));
        title('v');
        subplot(324)
        plot(time,state_sim(:,4)*180/pi);
        title('\omega');
        subplot(3,2,[5 6]);
        title('F');
        stairs(time,controls_MPC(:,1));
        
    case 'ChainofMasses_Lin'
        figure(1);
        subplot(311);
        plot(time,controls_MPC(:,1));
        ylabel('$u_x$','Interpreter','latex');
        subplot(312);
        plot(time,controls_MPC(:,2));
        ylabel('$u_y$','Interpreter','latex');
        subplot(313);
        plot(time,controls_MPC(:,3));
        ylabel('$u_z$','Interpreter','latex');

        n = data.n;
        figure(2);
        plot3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],'Color',red,'LineStyle','--');
        hold on;
        grid on;
        plot3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],'Color',blue,'LineStyle','-');
        scatter3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],10,'MarkerFaceColor','none');
        scatter3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],10,'MarkerFaceColor',red);
        xlabel('X[m]');
        ylabel('Y[m]');
        zlabel('Z[m]');
        
        xlim([-0.5 8]);
        ylim([-1 1]);
        zlim([-6 0.5]);
        
   case 'ChainofMasses_NLin'
        figure(1);
        subplot(311);
        plot(time,controls_MPC(:,1));
        ylabel('$u_x$','Interpreter','latex');
        subplot(312);
        plot(time,controls_MPC(:,2));
        ylabel('$u_y$','Interpreter','latex');
        subplot(313);
        plot(time,controls_MPC(:,3));
        ylabel('$u_z$','Interpreter','latex');

        n = data.n;
        figure(2);
        plot3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],'Color',red,'LineStyle','--');
        hold on;
        grid on;
        plot3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],'Color',blue,'LineStyle','-');
        scatter3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],10,'MarkerFaceColor','none');
        scatter3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],10,'MarkerFaceColor',red);
        xlabel('X[m]');
        ylabel('Y[m]');
        zlabel('Z[m]');
        
        xlim([-1.2 1.2]);
        ylim([-1.2 1.2]);
                
    case 'TethUAV'
         
        phi_ref = input.od(1,1);
        phi_ref = repmat(phi_ref, size(time));
        theta_ref = input.od(2,1);
        theta_ref = repmat(theta_ref, size(time));
        axes_ref = [];
        axes_lim = [];
        
        figure();      
        subplot(221)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,1)),'Color',red);
        plot(time(1:end),rad2deg(phi_ref),'k--');
        title('\phi');
        legend('\phi','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(222)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,2)),'Color',red);
        title('\phi_{dot}');
        legend('\phi_{dot}')
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
         
        subplot(223)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,3)),'Color',red);
        plot(time(1:end),rad2deg(theta_ref),'k--');
        title('\theta');
        legend('\theta','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(224)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,4)),'Color',red);
        title('\theta_{dot}');
        legend('\theta_{dot}');%,'ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        % set axes limits, the same for all the plots:
        maxY = max(axes_lim(:,2));
        minY = min(axes_lim(:,1));
        for i = 1 : length(axes_lim)
            cur_ax = axes_ref(i);
            cur_ax.YLim = [minY maxY];
        end
        
        figure()
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,5),'Color',red);
        title('f1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,6),'Color',red);
        title('f2');
        
        figure();
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,1),'Color',red);
        title('df1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,2),'Color',red);
        title('df2');
        
        figure();
        grid on;
        plot(time(1:end-1),constraints(:,1),'Color',red);
        title('fL');

        % plot time statistics
        
        figure();
        hold on;
        grid on;
        plot(time(2:end-1)', CPT(2:end, 1)); % cpt, tshooting, tcond, tqp
        title('Time statistics');
        xlabel('[s]')
        ylabel('[ms]')
        
        % plot nonlinear cost fcn terms       
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,1)+state_sim(:,3)));
        plot(time', rad2deg(pi/2)*ones(size(time)));
        title('Cost fcn: Avoid singularity');
        legend('\phi + \theta', '\pi/2');
        xlabel('[s]')
        ylabel('[deg]')
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,1)), 'k');
        plot(time, rad2deg(phi_ref), 'k--');
        plot(time', rad2deg(state_sim(:,3)), 'r');
        plot(time, rad2deg(theta_ref), 'r--');
        title('Cost fcn: Attitude behavior close to the ground');
        legend('\phi', '\phi_{ref}', '\theta', '\theta_{ref}');
        xlabel('[s]')
        ylabel('[deg]')
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,2)), 'r');
        plot(time', rad2deg(state_sim(:,4)), 'b');
        plot(time', rad2deg(state_sim(:,1)), 'k');
        plot(time, rad2deg(phi_ref), 'k--');
        title('Cost fcn: \phi and \theta velocities close to the ground');
        leg = legend('$\dot{\phi}$', '$\dot{\theta}$', '$\phi$', '$\phi_{ref}$');
        set(leg,'Interpreter','latex');
        xlabel('[s]')
        ylabel('[deg]')
        
    case 'DiM'	

         samples=size(y_sim,1);	

         figure;	
        title('MCA Tracking of perceived signals');	

         subplot(3,2,1);	
        plot(y_sim(:,1),'r');	
        hold on;	
        plot(data.REF(1:samples,1),'k');	
        title('Longitudinal: $\hat{a}_x$','Interpreter','latex');	

         subplot(3,2,2);	
        plot(y_sim(:,2),'r');	
        hold on;	
        plot(data.REF(1:samples,2),'k');	
        title('Lateral: $\hat{a}_y$','Interpreter','latex');	

         subplot(3,2,3);	
        plot(y_sim(:,3),'r');	
        hold on;	
        plot(data.REF(1:samples,3),'k');	
        title('Vertical: $\hat{a}_z$','Interpreter','latex');	

         subplot(3,2,4);	
        plot(y_sim(:,4),'r');	
        hold on;	
        plot(data.REF(1:samples,4),'k');	
        title('Roll: $\hat{\omega}_{\psi}$','Interpreter','latex');	

         subplot(3,2,5);	
        plot(y_sim(:,5),'r');	
        hold on;	
        plot(data.REF(1:samples,5),'k');	
        title('Pitch: $\hat{\omega}_{\theta}$','Interpreter','latex');	

         subplot(3,2,6);	
        plot(y_sim(:,6),'r');	
        hold on;	
        plot(data.REF(1:samples,6),'k');	
        title('Yaw: $\hat{\omega}_{\phi}$','Interpreter','latex');	


         figure;	
        subplot(3,2,1)	
        plot(y_sim(:,13));	
        hold on;	
        plot(y_sim(:,7),'r');	
        plot(zeros(Tf*100,1),'k');	
        title('Longitudinal displacement')	
        lgd=legend('tripod: $p_{x,T}$','hex: $p_{x,H}$','ref: $p_{x,T}$');	
        set(lgd,'Interpreter','latex');	

        subplot(3,2,2)	
        plot(y_sim(:,14));	
        hold on;	
        plot(y_sim(:,8),'r');	
        plot(zeros(Tf*100,1),'k');	
        title('Lateral displacement')	
        lgd=legend('tripod: $p_{y,T}$','hex: $p_{y,H}$','ref: $p_{y,T}$');	
        set(lgd,'Interpreter','latex');	

        subplot(3,2,3)	
        plot(y_sim(:,9),'r');	
        hold on;	
        plot(zeros(Tf*100,1),'k');	
        title('Vertical displacement: $p_{z,H}$','Interpreter','latex');	

        subplot(3,2,4)	
        plot(y_sim(:,20));	
        hold on;	
        plot(y_sim(:,17),'r');	
        plot(zeros(Tf*100,1),'k');	
        title('Yaw')	
        lgd=legend('tripod: $\phi_T$','hex: $\phi_H$','ref');	
        set(lgd,'Interpreter','latex');	

         subplot(3,2,5)	
        plot(y_sim(:,18),'r');	
        hold on;	
        plot(zeros(Tf*100,1),'k');	
        title('Pitch');	
        lgd=legend('hexpod: $\theta_H$','ref');	
        set(lgd,'Interpreter','latex');	

        subplot(3,2,6)	
        plot(y_sim(:,19),'r');	
        hold on;	
        plot(zeros(Tf*100,1),'k');	
        title('Roll');	
        lgd=legend('hexpod: $\psi_H$','ref');	
        set(lgd,'Interpreter','latex');	

        figure;	
        title('Hex actuator constraints')	
        plot(constraints(:,1:6));	
        hold on;	
        plot(1.045*ones(samples,1),':');	
        plot(1.375*ones(samples,1),':');	
        axis([0 mem.iter 1.0 1.4]);	
        title('Hexpod actuator constraints');
        
    case 'TurboEngine'
        
        figure()
        ax1 = subplot(4,1,1);
        hold on
        stairs(time,state_sim(:,3),'Color',blue);
        stairs(time,state_sim(:,4),'--','Color',red);
        ylabel('actuation / %');
        legend('u1', 'u2');
        grid on;

        ax2 = subplot(4,1,2);
        hold on
        plot(time(1:end-1),constraints(:,1),'Color',blue);
        plot(time(1:end-1),data.REF(1)*ones(1,length(time)-1),'k--');
        ylabel('charging pressure / bar');
        legend('output', 'reference')
        grid on;

        ax3 = subplot(4,1,3);
        hold on
        plot(time(1:end-1),constraints(:,2)*60,'Color',blue);
        plot(time(1:end-1),constraints(:,3)*60,'Color',red);
        plot(time(1:end-1),90e3*ones(1,length(time)-1),'--','Color',blue);
        plot(time(1:end-1),180e3*ones(1,length(time)-1),'--','Color',red);
        ylabel('turbocharger speeds / min^{-1}');
        grid on;
        legend('lp', 'hp', 'limit');

        ax4 = subplot(4,1,4);
        hold on
        plot(time,state_sim(:,1),'Color',blue);
        plot(time,state_sim(:,2),'--','Color',red);
        ylabel('x / -');
        legend('x1', 'x2');
        grid on;
        
        xlabel('Time[s]');

    case 'highFidelityTiltingQuadrotor'
        
        if( exist('fig_trj','var') && isgraphics(fig_trj,'figure'))
            delete(fig_trj);
        end
        fig_trj = figure('Name','Trajectories');        
        t = tiledlayout(4,3);
        
        nexttile(1,[2 1]);  % position
        plot(time,state_sim(:,1:3));
        hold on
        if opt.ref_type == 2
            plot(time,data.REF(1:length(time),1:3),'--');
        elseif opt.ref_type == 0
            plot(time, [data.REF(1).*ones(length(time),1), data.REF(2).*ones(length(time),1), data.REF(3).*ones(length(time),1)], '--')
        end
        grid on
        ylabel('$p$ [m]');
        legend('$p_x$','$p_y$','$p_z$')
        
        nexttile(7,[2 1]); % quaternion
        plot(time, rad2deg(quat2eul(state_sim(:,4:7))));
        %ylim([-1.1 1.1]);
        grid on
        ylabel('$eul$ [deg]');
        legend("yaw","pitch","roll");
        
        nexttile(2,[2 1]);  % linear vel
        plot(time,state_sim(:,8:10));
        grid on
        ylabel('$v$ [m/s]');
        legend('$v_x$','$v_y$','$v_z$')
        
        nexttile(8,[2 1]);  % angular vel
        plot(time,state_sim(:,11:13));
        grid on
        ylabel('$\omega$ [rad/s]');
        legend('$\omega_x$','$\omega_y$','$\omega_z$')
        
        nexttile(3,[1 1]);  % alpha
        plot(time, rad2deg(state_sim(:,14:17)));
        grid on
        ylabel('$\alpha$ [rad]');
        legend('$\alpha_1$','$\alpha_2$','$\alpha_3$','$\alpha_4$')
        
        nexttile(6,[1 1]);  % alpha dot
        plot(time, rad2deg(controls_MPC(:,5:8)));
        grid on
        ylabel('$\dot{\alpha}$ [deg/s]');
        legend('$\dot{\alpha}_1$','$\dot{\alpha}_2$','$\dot{\alpha}_3$','$\dot{\alpha}_4$')
        
        nexttile(9,[1 1]);  % alpha u
        plot(time,controls_MPC(:,1:4));
        grid on
        ylabel('$\dot{\bar{w}}$');
        legend('$\dot{\bar{w}}_1$','$\dot{\bar{w}}_2$','$\dot{\bar{w}}_3$','$\dot{\bar{w}}_4$')

        nexttile(12,[1 1]);  % alpha u
        plot(time,state_sim(:,18:21));
        grid on
        ylabel('$\bar{w}$');
        legend('$\bar{w}_1$','$\bar{w}_2$','$\bar{w}_3$','$\bar{w}_4$')
        
        xlabel(t,'$t$ [s]','interpreter','latex');

end