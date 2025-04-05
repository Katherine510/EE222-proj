# UC Berkeley EE222/ME237 Nonlinear Systems Ball and Beam Project

EE222/ME237 Nonlinear Systems, Spring 2025 Starter code and instructions for the course project.

## Project Results:

### LQR Controller

| Metric  | Sinusoidal | Square | 
| ------------- | ------------- | ------------- |
| Average Tracking Error  | 0.0004  | 0.0015 |
| Average Energy Consumption  | 0.0216  | 0.0741 |
| Tracking Cost  | 0.72  | 2.64 |
| Energy Cost  | 0.11  | 0.37 |
| Total Score  | 0.83  | 3.01 |


[LQR Sinusoidal video](https://youtu.be/mS_5MZMwuOk) |
[LQR Square video](https://youtu.be/Lpghpy8DOZk)


### LQG Controller

| Metric  | Sinusoidal | Square | 
| ------------- | ------------- | ------------- |
| Average Tracking Error  | 0.0004  | 0.0014 |
| Average Energy Consumption  | 0.0253  | 0.0860 |
| Tracking Cost  | 0.70  | 2.60 |
| Energy Cost  | 0.13  | 0.43 |
| Total Score  | 0.83  | 3.03 |


[LQG Sinusoidal video](https://youtu.be/EFeHQGZdgwg) |
[LQG Square video](https://youtu.be/mWMFsFKKsPc)


## Project Overview

This project involves designing and testing nonlinear controllers for a ball and beam system. The objective is to develop controllers that stabilize the ball at a desired position on the beam. You will first implement your controllers in MATLAB simulations and later test them on physical hardware.

## Understanding the Problem

To gain a full understanding of the problem and project expectations, please refer to the following documents in this repository:

[`EE_222_Course_Project.pdf`](EE_222_Course_Project.pdf) – Overview of the project and system model. (Disregard the due dates and GitHub link in this older document)

[`EE222 Lab Feedback and FAQ.pdf`](EE222_Lab_Feedback_and_FAQ.pdf) – Common issues and recommendations.

[`EE222_Lab_Part_1_Simulation.pdf`](EE222_Lab_Part_1_Simulation.pdf) – Instructions for running simulation.

[`EE222_Lab_Part_2_Hardware_Testing.pdf`](EE222_Lab_Part_2_Hardware_Testing.pdf) – Instructions for hardware testing. (To be updated)

## Code Instructions

### Prerequisites

Install MATLAB and Simulink using the Berkeley academic license.

### Getting Started

Clone or fork this repository.

Run `setup.m` or manually add the repository and its subfolders to the MATLAB path.

Modify only studentControllerInterface.m to implement your controller.

To test your controller:

Run `run_matlab_ball_and_beam.m` for a MATLAB-based simulation.

Run `run_simulink_ball_and_beam.m` for a Simulink-based simulation.
