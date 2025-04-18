classdef studentControllerInterface < matlab.System
    properties (Access = private)

        t_prev = -1;
        theta_d = 0;
        p_prev = 0;
        theta_prev = 0;
        V_servo = 0;

        g = 9.81;
        rg = 0.0254;
        L = 0.4255;
            
        B = [0; 0; 0; 1.5/0.025];
        C = [1, 0, 0, 0; 0, 0, 1, 0];

        Q = [75,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0];
        R = 1;


        QN = [0.025,0,0,0;
              0,0.05,0,0;
              0,0,0.025,0;
              0,0,0,0.05];

        RN = [0.025 0;
              0 0.025];

        QWV = blkdiag([0.025,0,0,0;
                       0,0.05,0,0;
                       0,0,0.025,0;
                       0,0,0,0.05], ...
                       [0.025 0;
                       0 0.025]);

        x_hat = [-0.19; 0.00; 0; 0];
        x_est = [-0.19, 0.00, 0, 0];
        P = [0.00467317845216159	0.00447778693304056	1.79600514031296e-06	-3.34363679122678e-09;
            0.00447778693210613	0.104399968275514	8.18335789697825e-05	-3.66928669353362e-08;
            1.79600044138968e-06	8.18335148621331e-05	0.00447750476234807	2.00500225891637e-05;
            -3.34797749744990e-09	-3.66963222880958e-08	2.00500225883704e-05	0.00124997492495667];
        M = eye(4);

        W = 0.01 * [1, 0, 0, 0;
                     0, 1, 0, 0; 
                     0, 0, 1, 0; 
                     0, 0, 0, 1];
        V = 0.01 * eye(2);
    end
    
    methods(Access = protected)

        function A = linA(obj, x)
            x1 = x(1);
            x2 = x(2);
            x3 = x(3);
            x4 = x(4);

            A = [0, 1,                                            0,                      0;
                 0, 0, 0.0051*cos(x3)*sin(x3)*x4^4 + 0.4183*cos(x3), -0.0102*x4^3*cos(x3)^2;
                 0, 0,                                            0,                      1;
                 0, 0,                                            0,                    -40];

        end
       
        %% Main Controller Interface
        function [V_servo, x_p, P_p] = stepImpl(obj, t, p_ball, theta)
            
            % State Estimation
            xg = generic_SE(obj, t, p_ball, theta);
            xk = kalmanFilter(obj, t, p_ball, theta)';

            % Feedback Controller
            x_p = obj.x_hat;
            P_p = obj.P;
            % V_servo = stepImplP(obj, t, xk);
            V_servo = stepImplLQR(obj, t, xg);
            % V_servo = stepImplLQG(obj, t, xg);

            dV = V_servo - obj.V_servo;
            if abs(dV) > 0.04
                V_servo = 0.04 * sign(dV) + obj.V_servo;
            elseif abs(dV) > 0.005
                V_servo = obj.V_servo;
            end

            if V_servo > 5
               V_servo = 5;
            elseif V_servo < -5
                V_servo = -5;
            end

            obj.V_servo = V_servo;
            obj.t_prev = t;
            obj.p_prev = p_ball;
            obj.theta_prev = theta;

        end

        %% State Estimation: Extended Kalman Filter
        function x = kalmanFilter(obj, t, p_ball, theta)
            % Get some data
            y = [p_ball; theta];
            dt = t - obj.t_prev;

            x = obj.x_hat;

            if dt > 0

                % Calculate kalman states
                A_lin = linA(obj, x);
    
                % Predict
                x_p = obj.x_hat + (A_lin*obj.x_hat + obj.B*obj.V_servo)*dt;
                P_p = obj.P + (A_lin*obj.P+obj.P*A_lin' + eye(4)*1)*dt;
    
                y_p = obj.C*x_p;
                innov = y - y_p;
                S = obj.C*P_p*obj.C'+ eye(2)*1;
    
                K = obj.P*obj.C'/S;
    
                x = x_p + K*innov;
                obj.x_hat = x;
                obj.P = (eye(4) - K*obj.C)*P_p;
            end
        end

        %% State Estimation: Generic Time Stepping
        function x = generic_SE(obj, t, p_ball, theta)
            dt = t - obj.t_prev;

            p_vel = (p_ball-obj.p_prev)/dt;
            theta_vel = (theta-obj.theta_prev)/dt;
            if dt == 0
                p_vel = 0;
                theta_vel = 0;
            end
            
            x = [p_ball, p_vel, theta, theta_vel];
     
            obj.x_est = x;

        end
        
        %% Feedback Controller: P Controller
        function V_servo = stepImplP(obj, t, x)
            p_ball = x(1);
            theta = x(3);

            t_prev = obj.t_prev;
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            k_p = 3;
            theta_d = - k_p * (p_ball - p_ball_ref);

            theta_saturation = 56 * pi / 180;    
            theta_d = min(theta_d, theta_saturation);
            theta_d = max(theta_d, -theta_saturation);

            k_servo = 10;
            V_servo = k_servo * (theta_d - theta);
            
            obj.t_prev = t;
            obj.theta_d = theta_d;
        end

        %% Feedback Controller: LQR
        function V_servo = stepImplLQR(obj, t, x)  
            
            % fetch the previous values
            t_prev = obj.t_prev;
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            dt = t - t-t_prev;
            
            A_lin = linA(obj, x);
            x = x - [p_ball_ref, v_ball_ref, 0, 0];
            
            % coder.extrinsic('lqr')
            K = [8.6603 10.0662 2.4465 0.0586];
            % K = lqr(A_lin, obj.B, obj.Q, obj.R);

            % dP = -(A_lin'*obj.P + obj.P*A_lin - (obj.P*obj.B)/obj.R*(obj.B'*obj.P) + obj.Q);
            % P = dt*dP + obj.P;
            % 
            % K = inv(obj.R)*(obj.B'*P);
            % 
            % obj.P = P;
            
            coder.extrinsic('icare')
%             P = eye(4);
%            [P,~,~] = icare(A_lin,obj.B,obj.Q,obj.R,[],[],[]);
            
            % K = inv(obj.R)*(obj.B'*P);

            % disp(K)

            V_servo = -K * x';
        end

        %% Feedback Controller: LQG 
        function V_servo = stepImplLQG(obj, t, x)
            % fetch the previous values
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            
            A_lin = linA(obj, x);
            x = x - [p_ball_ref, v_ball_ref, 0, 0];
            disp(A_lin)
            
            coder.extrinsic('ss')
            sys = ss(A_lin, obj.B, obj.C, 0);
            
            if abs(x(1)) > 0.01*0.5
                obj.Q = [365,0,0,0;
                         0,75,0,0;
                         0,0,0,0;
                         0,0,0,0];
            else
                obj.Q = [550,0,0,0;
                         0,150,0,0;
                         0,0,0,0;
                         0,0,0,0];
            end

            QXU = blkdiag(obj.Q, obj.R);
            INFO = struct('Kx',[0, 0, 0, 0], ...
                          'Kw',[0, 0, 0, 0] , ...
                          'Ki',[] , ...
                          'L', [0,0; 0,0; 0,0; 0,0] , ...
                          'Mx',[0,0; 0,0; 0,0; 0,0] , ...
                          'Mw',[0,0; 0,0; 0,0; 0,0]);

            coder.extrinsic('lqg')
            [~,INFO] = lqg(sys,QXU,obj.QWV);

            V_servo = -INFO.Kx * x';
            
        end
    end
    
    methods(Access = public)
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end

    end
    
end
