% Exact solution for Buckley leverett with two Shock-Rarefaction consists
% of 7 segments.  A brief description of each line segment:
%
% Segment 1:
%
%    straight line from (-1,0) to (-0.5, 0 ).
%
% Segment 2:
%
%    Rarefaction connecting (-0.5,0) to (ls, qsl).
%
%    The naming convention here is ls = 'left shock', and qsl = 'qstar_left'.
%
%    With qsl = 0.1339745962155613, we get ls = -0.5 + dt * fp( qsl ).
%
%    The rarefaction can be plotted by simply plotting the points:
% 
%          [-0.5 + dt*fp( q ), q],   where q = linspace(0,qsl).
%
% Segment 3:
%
%    straight line connecting (ls, qsl) to (ls, 1).
%
% Segment 4:
%
%    straight line connecting (ls,1) to (0,1)
%
% Segment 5:
%
%    rarefaction connecting (0,1) to (rs, qsr).
%
%    rs = location of 'right shock', and qsr is the solution to the R-H
%    conditions.
%
%    Again, these points can be plotted by graphing
%
%           [dt*fp(q), q], with q a sampling of the interval (qsr, 1).
%
% Segment 6:
%
%    straight line connecting (rs, qsr) to (rs, 0 ).
%
% Segment 7:
%
%    straight line connecting (rs,0) to (1,0)
%

qs_left_problem  = 0.1339745962155613;
qs_right_problem = 0.5;

q = linspace(qs_right_problem, 1.0 );
dt = 0.4;

figure(1);
clf;

ql = linspace(0., qs_left_problem );
left_shock = -0.5 + dt * buckley_fp( qs_left_problem );
fl = -0.5 + dt*buckley_fp( ql );

% segment 1
plot( [-1 -0.5], [0 0], '-r' )
hold on;

% segment 2
plot( fl, ql, '-r' );
hold on;

% segment 3
plot( fl, ql, '-r' );
plot( [left_shock left_shock], [qs_left_problem 1.0], '-r' );
hold on;

% segment 4
plot( [left_shock 0], [1 1], '-r' );


% segment 5
plot( dt*buckley_fp(q), q, '-r'  );
hold on;

% segment 6
plot( [ dt*buckley_fp(qs_right_problem) dt*buckley_fp(qs_right_problem) ], ...
      [ qs_right_problem 0 ], '-r' );
hold on;

% segment 7
plot( [dt*buckley_fp(qs_right_problem) 1.0], [0 0], '-r' );

hold on;
plot( [-1 1], [0 0], '--k' );

axis( [-1 1 -0.1 1.1] );
hold off;

