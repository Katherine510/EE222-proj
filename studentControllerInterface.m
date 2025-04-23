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
        P = eye(4);
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

        function K = lqr_custom(obj, A, B, Q, R)
        %LQR_VIA_MATRIX_SIGN_FUNCTION Computes the LQR gain using matrix sign function method
        %
        %   K = LQR_VIA_MATRIX_SIGN_FUNCTION(A, B, Q, R) returns the optimal gain matrix K
        %   for the continuous-time LQR problem:
        %
        %       minimize J = ∫ (x'Qx + u'Ru) dt
        %       subject to dx/dt = Ax + Bu
        %
        %   Inputs:
        %       A - System dynamics matrix (n x n)
        %       B - Input matrix (n x m)
        %       Q - State cost matrix (n x n), symmetric positive semi-definite
        %       R - Input cost matrix (m x m), symmetric positive definite
        %
        %   Output:
        %       K - Optimal state feedback gain matrix (m x n)
        
            % Invert R (assumes R is positive definite)
            R_inv = inv(R);
            G = B * R_inv * B';
        
            % Construct the Hamiltonian matrix Z
            Z = [ A      -G;
                 -Q   -A' ];
        
            % Initialize matrix W for iteration
            % W = Z;
        
            % Newton iteration to compute the matrix sign function
            for i = 1:1000
                Z = Z - 0.5 * (Z - inv(Z));
            end
        
            % Determine the size of the system
            n = size(A, 1);
        
            % Partition W into 4 submatrices
            W11 = Z(1:n, 1:n);
            W12 = Z(1:n, n+1:end);
            W21 = Z(n+1:end, 1:n);
            W22 = Z(n+1:end, n+1:end);
        
            % Solve for the unique positive semidefinite solution P to the Riccati equation
            M = [W12; W22 + eye(n)];
            N = [W11 + eye(n); W21];
            P = M \ (-N);
        
            % Compute the optimal LQR gain
            K = R_inv * B' * P;
        end
       
        %% Main Controller Interface
        function V_servo = stepImpl(obj, t, p_ball, theta)
            
            % State Estimation
            xg = generic_SE(obj, t, p_ball, theta);
            xk = kalmanFilter(obj, t, p_ball, theta)';

            % Feedback Controller
            % V_servo = stepImplP(obj, t, xk);
            V_servo = stepImplLQR(obj, t, xg);
            % V_servo = stepImplLQG(obj, t, xg);

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
            P_p = obj.P + (A_lin*obj.P+obj.P*A_lin' + eye(4)*0.1)*dt;

            y_p = obj.C*x_p;
            innov = y - y_p;
            S = obj.C*P_p*obj.C'+ eye(2)*0.1;

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
            
            coder.extrinsic('lqr')
            K = [0, 0, 0, 0];
            K = lqr_custom(obj, A_lin, obj.B, obj.Q, obj.R);

            % dP = -(A_lin'*obj.P + obj.P*A_lin - (obj.P*obj.B)/obj.R*(obj.B'*obj.P) + obj.Q);
            % P = dt*dP + obj.P;
            % 
            % K = inv(obj.R)*(obj.B'*P);
            % 
            % obj.P = P;
            
            % coder.extrinsic('icare')
            % P = eye(4);
            % [P,~,~] = icare(A_lin,obj.B,obj.Q,obj.R,[],[],[]);
            % 
            % K = inv(obj.R)*(obj.B'*P);

            V_servo = -K * x';
            obj.V_servo = V_servo;
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
            obj.V_servo = V_servo;
        end
    end
    
    methods(Access = public)
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end

    end
    
end
