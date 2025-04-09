classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;
        p_prev = 0;
        theta_prev = 0;

        state = 0;

        V_servo = 0;

        g = 9.81;
        rg = 0.0254;
        L = 0.4255;
            
        B = [0; 0; 0; 1.5/0.025];
        C = [1, 0, 0, 0; 0, 0, 1, 0];

        Q = [260,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0];
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
                      0,0,0,0.05], [0.025 0;
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
        %% Sample Controller: Simple Proportional Controller
        function V_servo = stepImpl(obj, t, p_ball, theta)
        % This is the main function called every iteration. You have to implement
        % the controller in this function, bu you are not allowed to
        % change the signature of this function. 
        % Input arguments:
        %   t: current time
        %   p_ball: position of the ball provided by the ball position sensor (m)
        %
        %   theta: servo motor angle provided by the encoder of the motor (rad)
        % Output:
        %   V_servo: voltage to the servo input.        
            
            xg = generic_SE(obj, t, p_ball, theta);
            V_servo = stepImplLQR(obj, t, xg);
            % V_servo = stepImplLQG(obj, t, xg);

            obj.t_prev = t;
            obj.p_prev = p_ball;
            obj.theta_prev = theta;

            
        end


        %% Course Controller
        function x = kalmanFilter(obj, t, p_ball, theta)
            % Get some data
            y = [p_ball; theta];
            dt = t - obj.t_prev;

            x = obj.x_hat;

            % Calculate kalman states
            x1 = x(1);
            x2 = x(2);
            x3 = x(3);
            x4 = x(4);

            A_lin = [0, 1,                                            0,                      0;
                 0, 0, 0.0051*cos(x3)*sin(x3)*x4^4 + 0.4183*cos(x3), -0.0102*x4^3*cos(x3)^2;
                 0, 0,                                            0,                      1;
                 0, 0,                                            0,                    -40];

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
        
        
        %% Test Controller: LQR Controller
        function V_servo = stepImplLQR(obj, t, x)
        % This is the main function called every iteration. You have to implement
        % the controller in this function, bu you are not allowed to
        % change the signature of this function. 
        % Input arguments:
        %   t: current time
        %   p_ball: position of the ball provided by the ball position sensor (m)
        %
        %   theta: servo motor angle provided by the encoder of the motor (rad)
        % Output:
        %   V_servo: voltage to the servo input.     
            
            % fetch the previous values
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            
            x1 = x(1);
            x2 = x(2);
            x3 = x(3);
            x4 = x(4);

            A_lin = [0, 1,                                            0,                      0;
                 0, 0, 0.0051*cos(x3)*sin(x3)*x4^4 + 0.4183*cos(x3), -0.0102*x4^3*cos(x3)^2;
                 0, 0,                                            0,                      1;
                 0, 0,                                            0,                    -40];

            x = x - [p_ball_ref, v_ball_ref, 0, 0];
            coder.extrinsic('lqr')
            K = [0, 0, 0, 0];
            K = lqr(A_lin, obj.B, obj.Q, obj.R);
            
            
            V_servo = -K * x';
            obj.V_servo = V_servo;
        end

        %% Test Controller: LQG Controller
        function V_servo = stepImplLQG(obj, t, x)
        % This is the main function called every iteration. You have to implement
        % the controller in this function, bu you are not allowed to
        % change the signature of this function. 
        % Input arguments:
        %   t: current time
        %   p_ball: position of the ball provided by the ball position sensor (m)
        %
        %   theta: servo motor angle provided by the encoder of the motor (rad)
        % Output:
        %   V_servo: voltage to the servo input.        
            
            % fetch the previous values
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            
            x1 = x(1);
            x2 = x(2);
            x3 = x(3);
            x4 = x(4);

            A_lin = [0, 1,                                            0,                      0;
                 0, 0, 0.0051*cos(x3)*sin(x3)*x4^4 + 0.4183*cos(x3), -0.0102*x4^3*cos(x3)^2;
                 0, 0,                                            0,                      1;
                 0, 0,                                            0,                    -40];
            x = x - [p_ball_ref, v_ball_ref, 0, 0];
            
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
            coder.extrinsic('lqg')
            [KLQG,INFO] = lqg(sys,QXU,obj.QWV);


            V_servo = -INFO.Kx * x';
            obj.V_servo = V_servo;
        end
    end
    
    methods(Access = public)
        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)        
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end

    end
    
end
