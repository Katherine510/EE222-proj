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

        A = 0;
        B = 0;
        C = 0;
        Q = 0;
        R = 0;

        QN = 0;
        RN = 0;
        QWV = 0;

        x_hat = [-0.19; 0.00; 0; 0];
        M = 0;
        W = 0;
        V = 0;

        V_servo = 0;
    end
    methods(Access = protected)

        %% Sample Controller: Simple Proportional Controller
        function V_servo = stepImplP(obj, t, p_ball, theta)
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
            
            t_prev = obj.t_prev;
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            % Decide desired servo angle based on simple proportional feedback.
            k_p = 3;
            theta_d = - k_p * (p_ball - p_ball_ref);

            % Make sure that the desired servo angle does not exceed the physical
            % limit. This part of code is not necessary but highly recommended
            % because it addresses the actual physical limit of the servo motor.
            theta_saturation = 56 * pi / 180;    
            theta_d = min(theta_d, theta_saturation);
            theta_d = max(theta_d, -theta_saturation);

            % Simple position control to control servo angle to the desired
            % position.
            k_servo = 10;
            V_servo = k_servo * (theta_d - theta);
            
            % Update class properties if necessary.
            obj.t_prev = t;
            obj.theta_d = theta_d;

            
        end


        %% Course Controller
        function x = kalmanFilter(obj, t, p_ball, theta)
            % Get some data
            y = [p_ball; theta];
            dt = t - obj.t_prev;
            
            % Calculate kalman states
            syms x1 x2 x3 x4;
            A_lin = double(subs(obj.A, [x1, x2, x3, x4], obj.x_hat'));
            
            L = -obj.M*obj.C'*obj.V;
            dx_hat = A_lin*obj.x_hat + obj.B*obj.V_servo + L * (y - obj.C*obj.x_hat);
            % obj.x_hat = obj.x_hat + dx_hat*dt;
            obj.x_hat = [p_ball, dx_hat(1), theta, dx_hat(3)]';
            

            dM = A_lin*obj.M + obj.M*A_lin' + obj.W - obj.M*obj.C'*obj.V*obj.C*obj.M;
            obj.M = obj.M + dM*dt;
            x = obj.x_hat';

        end

        function x = generic_SE(obj, t, p_ball, theta)
            dt = t - obj.t_prev;

            p_vel = (p_ball-obj.p_prev)/dt;
            theta_vel = (theta-obj.theta_prev)/dt;

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
            
            syms x1 x2 x3 x4;
            A_lin = double(subs(obj.A, [x1, x2, x3, x4], x));
            x = x - [p_ball_ref, v_ball_ref, 0, 0];
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
            
            syms x1 x2 x3 x4;
            A_lin = double(subs(obj.A, [x1, x2, x3, x4], x));
            x = x - [p_ball_ref, v_ball_ref, 0, 0];

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
            [KLQG,INFO] = lqg(sys,QXU,obj.QWV);


            V_servo = -INFO.Kx * x';
            obj.V_servo = V_servo;
        end
    end
    
    methods(Access = public)
        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)
            
            xg = generic_SE(obj, t, p_ball, theta);
            xk = kalmanFilter(obj, t, p_ball, theta);
            V_servo = stepImplLQR(obj, t, xg);
            % V_servo = stepImplLQG(obj, t, xg);
            theta_d = obj.theta_d;


            obj.t_prev = t;
            obj.p_prev = p_ball;
            obj.theta_prev = theta;
        end


        function setupMODULE(obj)
            disp("You can use this function for initializaition.");
            g = 9.81;
            rg = 0.0254;
            L = 0.4255;
            K = 1.5;
            t = 0.025;
            
            syms x1 x2 x3 x4;
            eq = [ x2; 
                   (5/7)*(rg/L)*g*sin(x3) - (5/7)*(rg/L)^2 * x4^4 * (cos(x3))^2; 
                   x4; 
                   -(x4/t)];
            obj.A  = jacobian(eq, [x1, x2, x3, x4]);
            obj.B = [0; 0; 0; K/t];
            obj.C = [1, 0, 0, 0; 0, 0, 1, 0];

            obj.Q = [260,0,0,0;
                     0,0,0,0;
                     0,0,0,0;
                     0,0,0,0];
            obj.R = 1;


            obj.QN = [0.025,0,0,0;
                      0,0.05,0,0;
                      0,0,0.025,0;
                      0,0,0,0.05];
            obj.RN = [0.025 0;
                      0 0.025];
            obj.QWV = blkdiag(obj.QN, obj.RN);




            obj.x_hat = [-0.19; 0.00; 0; 0];
            obj.M = eye(4);

            obj.W = 0.01 * [1, 0, 0, 0;
                     0, 1, 0, 0; 
                     0, 0, 1, 0; 
                     0, 0, 0, 1];
            obj.V = 0.01 * eye(2);
        end

    end
    
end
