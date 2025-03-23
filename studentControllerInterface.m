classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;
        extra_dummy1 = 0;
        extra_dummy2 = 0;
        state = 0;
    end
    methods(Access = protected)
        % function setupImpl(obj)
        %    disp("You can use this function for initializaition.");
        % end
        


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
        function V_servo = stepCourseImpl(obj, t, p_ball, theta)
            % This controller is used to initially center the ball
            %   - The idea is to use a lighter controller to first center 
            %     the ball on the desired value before fully
            %     folowing the trajectory
            k_p = 0.5;
            k_servo = 20;
            theta_saturation = 7 * pi / 180;

            t_prev = obj.t_prev;
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);

            theta_d = - k_p * (p_ball - p_ball_ref);
            theta_d = min(theta_d, theta_saturation);
            theta_d = max(theta_d, -theta_saturation);

            V_servo = k_servo * (theta_d - theta);
            obj.theta_d = theta_d;
        end

        
        %% Fine Controller
        function V_servo = stepFineImpl(obj, t, p_ball, theta)
            % This controller is used for the actual trajectory tracking
            %   - Once we are approximately around where the desired
            %   location is, we will then run a more advanced controller to
            %   actually follow the desired trajectory
            k_p = 22;
            k_servo = 100;
            theta_saturation = 56 * pi / 180;

            t_prev = obj.t_prev;
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);

            theta_d = - k_p * (p_ball - p_ball_ref);
            theta_d = min(theta_d, theta_saturation);
            theta_d = max(theta_d, -theta_saturation);

            V_servo = k_servo * (theta_d - theta);
            obj.theta_d = theta_d;
        end
        
        
        %% Test Controller: Simple Proportional Controller
        function V_servo = stepImplFunnyPD(obj, t, p_ball, theta)
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
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            if obj.state == 0
                V_servo = stepCourseImpl(obj, t, p_ball, theta);

                % Once we are close enough, we switch controllers
                obj.state = abs((p_ball - p_ball_ref)) < 0.01;
            elseif obj.state == 1
                V_servo = stepFineImpl(obj, t, p_ball, theta);
            end
            obj.t_prev = t;
        end
    end
    
    methods(Access = public)
        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)        
            V_servo = stepImplFunnyPD(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
    
end
