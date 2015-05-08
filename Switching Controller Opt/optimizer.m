clear *; close all; clc

global best_yet initials
best_yet = 0;
consts = get_consts();

% States: x=[y,z,th,psi,dy,dz,dth,dpsi,m]'
fuel = consts.m_nofuel+0.7*consts.max.m_fuel;
initials = [10 150 0 0 0 0 0 0 fuel;
    10 150 pi/2 0 0 0 0 0 fuel;
    10 25 0 0 0 0 0 0 fuel;
    10 25 pi/4 0 0 0 0 0 fuel;
    100 25 pi/4 0 0 0 0 0 fuel;
    100 450 pi 0 0 0 0 0 fuel;
    10 450 pi 0 0 0 0 0 fuel;
    100 1500 pi 0 0 0 0 0 fuel;
    100 1500 pi/2 0 0 0 0 0 fuel;
    100 700 pi/2 0 0 0 0 0 fuel];


opts.MaxFunEvals  = 50000;
opts.TolX = 1e-2;
opts.StopFitness = 1e-4;
opts.CMA.active = 1;
opts.Restarts = 0;
opts.LogPlot = 'off';
opts.StopOnStagnation = 'on';
opts = cmaes('defaults', opts);

sigma = [.001 .001 10 .001 .001 .001 1 .001 0.001 0.1 0.1...
         .001  10  1  .001 .001 .001 .001 .001 .001 .001 .1]';

x0 = [0.001 0.001 100 0.001 0.001 0.001 10 0.001 0.001 .5 2.0... 
      0.1 100 10 0.000001 0.001 0.001 .001 0.001 0.001 .1 10];

[XMIN,FMIN,COUNTEVAL,STOPFLAG,OUT,BESTEVER] = cmaes('costfun',x0,sigma,opts);
