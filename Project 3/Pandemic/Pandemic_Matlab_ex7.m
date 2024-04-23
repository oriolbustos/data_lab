
function Result = Pandemic_Matlab_ex7(gap, confined_pc)
% This function is an example of how to use Covid19_Simulator
% Use it wisely!
% Clear screen
clc;

%----------------------------------------%
%               PARAMETERS               %
%----------------------------------------%

% Maximum number of iterations
max_iter = 1000;
% Size of the grid  (it's a square grid)
size_x = 100; % size_x x size_x is the size of the grid in poitns
% Density of population
density = 0.15;
% Number of people in the simulation
people = round((size_x)^2*density);
% Border parameters:
inner_border = 1; % 0 no inner border, 1 with inner border
%gap = 7; % oberture in the imperfect border
% Health parameters:
sick_pc = 50; % percentage of sick people (all grid, or within inner border)
max_sick_time = 30; % time of recovery after being infected (in average)
% Movement parameters:
confinement = 0; % 0 no confinement, 1 confinment. 
%confined_pc = 50;% degree of confinement in percetage (people not moving)
% Plot parameters
do_plot = 1; % 0 plot the result, 1 don't plot it.
plot_epidemic_curves_check = false;

Input = struct('max_iter', max_iter, 'people', people,...,
                 'size_x', size_x,...,
                 'inner_border', inner_border,...,
                 'gap', gap,...,
                 'sick_pc', sick_pc,...,
                 'max_sick_time', max_sick_time,...,
                 'confinement', confinement,...,
                 'confined_pc',confined_pc,...,
                 'do_plot', do_plot);

%----------------------------------------%
%               FUNCTION USE             %
%----------------------------------------%

Result = Covid19_Simulator(Input, plot_epidemic_curves_check);

end
%--------------------------------------------------------------------------

function Result = Covid19_Simulator(Input, plot_epidemic_curves_check)
% This function simulates an epidemy caused by a virus in different
% scenarios: 
%   -  No measure against the epidemy.
%   -  Contention of the plague focus.
%   -  Isolation (or confinement) of people individually.
%   -  Contention plus isolation.
% The outcome of the simulation is a epidemilogy plot with the evolution
% of the population across the epidemy (non-infected, infected,
% recovered)using epidemiology curves.
%
% Input Arguments:
% All the inputs arguments are stored in a structure called Input.
% Input.max_iter: Max. iteration in the loop that simulates time evolution.
% Input.people: Number of people in the simulation.
% Input.size_x: Max. grid value.
% Input.inner_border: logical. 1 means that there is an inner borders,
%                        0 means that there isn't any.
% Input.gap: interger and scalar. Size of the 'hole' in the inner border.
% Input.sick_pc: real. Percentage of sick people.
% Input.max_sick_time: An integer. Maximum time infected.
% Input.confinement: logical. 0 for no confinement, 1 for confinement.
% Input.confined_pc: real. Percentage of sick people.
% Input.do_plot: logical. Plot(1) or not plot(the result)

% Output Arguments:
% All the inputs arguments are stored in a structure called Result.
% Result.X_border: Vector. x coordinates of the border.
% Result.Y_border: Vector. y coordinates of the border.
% Result.X_people_time: Matrix. x Coordinates of people across time.
% Result.Y_people_time: Matrix. y Coordinates of people across time.
% Result.health_color_time: 3-D array. Health status color across time.
% Result.infection_matrix: Matrix. 
                           % Who(row) infected whom (columns)
                           % 1 infected, 0 non-infected
% Result.non_infected: Vector. Number of non-infected cases along time.
% Result.infected : Vector. Number of infected cases along time.
% Result.recovered: Vector. Number of recovered cases along time.

% message to the user
disp(' ')
disp('%----------------------%')
disp('%    Covid_Simulator   %')
disp('%----------------------%')
disp(' ')

%----------------------------------------%
%            EXTRACT PARAMETERS          %
%----------------------------------------%

% Extract from input structure
max_iter = Input.max_iter;
people =  Input.people;
size_x = Input.size_x; 
inner_border = Input.inner_border;
gap = Input.gap; 
sick_pc = Input.sick_pc;
max_sick_time = Input.max_sick_time; 
confinement = Input.confinement; 
confined_pc = Input.confined_pc;
do_plot = Input.do_plot; 

%----------------------------------------%
%              SOME CHECKS               %
%----------------------------------------%

% Check if input values are valid
make_checks(max_iter, size_x, gap, inner_border, people, sick_pc, confinement, confined_pc, max_sick_time, do_plot);

%----------------------------------------%
%            INITIALIZATION              %
%----------------------------------------%

% Message to the user
disp('Now Initializing')
% Define a 2D grid (the patient's space)
[X, Y, X_rshp, Y_rshp] = make_grid(size_x);
% Create the borders for this space
[X_border, Y_border, X_inner, Y_inner] = make_border(X, Y, inner_border, gap);
% Create people within the space but not in the border at time 0. 
[X_people, Y_people] = make_people(people, X_border, Y_border, X_rshp, Y_rshp);
% Create people's health status, associate colors, set time since infection.
[health_status, health_color, sick_time, recover_time] = make_health(X_people, Y_people, inner_border, Y_inner, people, sick_pc, max_sick_time);
% Create people movement. A fraction of people can be quiet, following a quarantine. 
[vX_people, vY_people, y_inner_max] = make_movement(X_people, inner_border, Y_inner, people, confinement, confined_pc);

%----------------------------------------%
%               MAIN LOOP                %
%----------------------------------------%

% Message to the user
disp('Now Simulating')
% Create the vectors with people's position, accross time
X_people_time = zeros(max_iter, size(X_people, 2));
Y_people_time = zeros(max_iter, size(X_people, 2)); % X has the same dimensios as Y
% Create arrays with people's health status color, across time
health_color_time = zeros(size(health_color, 1),size(health_color, 2), max_iter);
% Create the vectors to count non_infected, infected and recovered people
non_infected = zeros(1, max_iter);
infected = zeros(1, max_iter);
recovered = zeros(1, max_iter);
% Create the infection matrix (who infects whom)
infection_matrix = zeros(people);

for h = 1: max_iter
    % New position of people for each iteration
    X_people_time(h, :) = X_people;
    Y_people_time(h, :) = Y_people;
    % New health status color for each iteration
    health_color_time(:, :, h) = health_color;
    % Count cases in this iteration
    non_infected(h) = sum(health_status == 1);
    infected(h) = sum(health_status == 2);
    recovered(h) = sum(health_status == 3);

    % PLOT at each iteration %
    %figure()
    %scatter(X_people_time(h,:), Y_people_time(h,:),[], health_color_time(:, :,1), 'Filled');
   %hold on
   % scatter(X_border, Y_border,[], [0.91 0.41 0.17],'*'); %orange frontier
   % xlabel('Coordinates in X axis')
   % ylabel('Coordinates in Y axis')
   % title(sprintf('State of the Simulation at Time Index %d',h))
   % M(h)= getframe;

    % Stop the loop in case no people is infected
    if infected(h) == 0
        last_h = h;
        break;
    end
    % Update the position of the people on the grid
    [X_people, Y_people, vX_people, vY_people] = update_position(X_people, Y_people, vX_people, vY_people, X_border, Y_border, X_inner, Y_inner);
    % Update health status and healt color from healt to sick
    [health_status, health_color, infection_matrix, sick_time] = update_sick(X_people, Y_people, y_inner_max, health_status, health_color, infection_matrix, sick_time);
    [health_status, health_color, sick_time] = update_recover(health_status, health_color, sick_time, recover_time, plot_epidemic_curves_check);
end

% Message to the user
disp('Now Generating Results')

%----------------------------------------%
%               HEALTH COLORS            %
%----------------------------------------%
% Select only the indexes where the number of infected > 0
% plus the first iteration where infected == 0
health_color_time = health_color_time(:, :, 1: h);

%----------------------------------------%
%            RESHAPE POSITIONS           %
%----------------------------------------%

% Select only the indexes where the number of infected > 0
% plus the first iteration where infected == 0
X_people_time = X_people_time(1: last_h, :);
Y_people_time = Y_people_time(1: last_h, :);
% Reshape the variables to be plotted properly
X_people_time = reshape(X_people_time, last_h, size(X_people, 2));
Y_people_time = reshape(Y_people_time, last_h, size(Y_people, 2));

%----------------------------------------%
%             EPIDEMIC CURVES            %
%----------------------------------------%

% Select only the indexes where the number of infected > 0
% plus the first iteration where infected == 0
non_infected = non_infected(1: last_h);
infected = infected(1: last_h);
recovered = recovered(1: last_h);
% Plot the epidemic curve according to the epidemic scenario
if plot_epidemic_curves_check
    plot_epidemic_curves(do_plot, non_infected, infected, recovered, inner_border, confinement)
end
%----------------------------------------%
%                RESULTS                 %
%----------------------------------------%
Result = struct('X_border', X_border, 'Y_border', Y_border ,...,
                 'X_people_time', X_people_time,...,
                 'Y_people_time', Y_people_time,...,
                 'health_color_time', health_color_time,...,
                 'infection_matrix', infection_matrix,...,
                 'non_infected',non_infected,...,
                 'infected', infected,...,
                 'recovered', recovered);
%Result = struct('X_border', X_border, 'Y_border', Y_border ,...,
%                 'X_people_time', X_people_time,...,
%                 'Y_people_time', Y_people_time,...,
%                 'health_color_time', health_color_time,...,
%                 'infection_matrix', infection_matrix,...,
%                 'non_infected',non_infected,...,
%                 'infected', infected,...,
%                 'recovered', recovered, ...
%                 'M', M);
             
% Message to the user
disp('End of the Simulation')
end

%--------------------------------------------------------------------------

%----------------------------------------%
%               FUNCTIONS                %
%----------------------------------------%

function answer = isneg(value)
% Function for detecting if the value is negative or not
    if value < 0
        answer = 1;
    else
        answer = 0;
    end
end 
function answer = isint(value)
% Function for detecting if the value is integer or not
    if floor(value)==value && value~=Inf
        answer = 1;
    else
        answer = 0;
    end
end
function make_checks(max_iter, size_x, gap, inner_border, people, sick_pc, confinement, confined_pc, max_sick_time, do_plot)
 % max_iter integer greater than 0
 if ~isint(max_iter) || isneg(max_iter) || max_iter == 0 || ischar(max_iter)
    error('The input variable '',max_iter'' must be an integer greater than 0.')
end
 % size_x integer greater than 0
if ~isint(size_x) || isneg(size_x) || size_x == 0 || ischar(size_x)
    error('The input variable ''size_x'' must be an integer greater than 0.')
end
% Allowed data type for gap.
if ~isint(gap) || isneg(gap) || gap == 0 || ischar(gap)
    error('The input variable ''gap'' must be an integer greater than 0.')
end
% Inner_border either 1 or 0
if inner_border ~= 0 && inner_border ~= 1
    error('The input variable ''inner_border'' must be either 0 or 1.')
end
if ~isint(people) || isneg(people) || people == 0 || ischar(people)
    error('The input variable ''people'' must be an integer greater than 0.')
end
% Sick_pc interge bewteen 1 and 100:
if ~isreal(sick_pc) || isneg(sick_pc) || sick_pc == 0 || ischar(sick_pc) || sick_pc > 100
    error('The input variable ''sick_pc'' must be an integer greater than 0 and lower or equal than 100.')
end
% Confinement either 1 or 0
if confinement ~= 0 && confinement ~= 1
    error('The input variable ''confinement'' must be either 0 or 1.')
end
% Confinement_pc integer bewteen 1 and 100:
if ~isreal(confined_pc) || isneg(confined_pc) || confined_pc == 0 || ischar(confined_pc) ||  confined_pc > 100
    error('The input variable ''confined_pc'' must be an integer greater than 0 and lower or equal than 100.')
end
%  max_sick_time integer greater than 0
if ~isint( max_sick_time) || isneg(max_sick_time) ||  max_sick_time == 0 || ischar(max_sick_time)
    error('The input variable ''max_sick_time'' must be an integer greater than 0.')
end
% do_plot has to be either 0 or 1
if do_plot ~= 0 && do_plot ~= 1
    error('The input variable ''do_plot'' must be equal to 0 or 1.')
end
 end

% To create the grid
function [X, Y, X_rshp, Y_rshp] = make_grid(size_x)
% This function creates a folded(X, Y) and a 
% unfolded grid (X_rshp, Y_rshp) to simulate the space
% where the patients will move, get infected, and recover.

% Note: the grid starts at the corner (0,0).

% Input Arguments:
% size_x: Max. grid value.
% Output Arguments:
% X: A matrix. x coordinates of the grid.
% Y: A matrix. x coordinates of the grid.
% X_rshp: A vector. x coordinates of the grid (unfolded).
% Y_rshp: A vector. y coordinates of the grid (unfolded).

% Grid creation grid
% Unfolded grid:
[X, Y] = meshgrid(1: size_x, 1: size_x);
% Folded grid:
X_rshp = reshape(X,1,size(X, 1) * size(X, 2));
Y_rshp = reshape(Y,1,size(Y, 1) * size(Y, 2));
end

% To create the border
function [X_border, Y_border, X_inner, Y_inner] = make_border(X, Y, inner_border, gap)
% This function creates the borders for the patients' space. 
% Outer and (optionally)inner borders for the simulation are defined.
% The inner border is not totally hermetic so the virus can escape from 
% the confinement region. This is controlled by the parameter 'gap'

% Input Arguments:
% size_x: Max. grid value.
% X: A matrix. x coordinates of the grid.
% Y: A matrix. x coordinates of the grid.
% inner_border: logical. 1 means that there is an inner borders,
%                        0 means that there isn't any.
% gap: interger and scalar. Size of the 'hole' in the inner border.
% Output Arguments:
% X_border: A vector. x coordinates of total border(unfolded).
% Y_border: A vector. y coordinates of total border(unfolded).
% X_inner: A vector. x coordinates of the inner border(unfolded).
% Y_inner: A vector. y coordinates of the inner border(unfolded).

% Initialization for the inner border and some checks:
y_inner_lim = 20;
x = X(1, :);
L = length(x);
L_half = floor(L/2) + 1;
gap_half = floor(gap/2);
% Allowed range of values for gap.
if (L_half - gap_half) < x(1) || (L_half + gap_half) > x(end)
    error('The selected value of ''gap'' is out of the range of the variable')
end
% Outer border X,Y coordinates
X1 = X(1, :);
Y1 = Y(1, :);
X2 = X(end, :);
Y2 = Y(end, :);
X3 = X(:, 1);
Y3 = Y(:, 1);
X4 = X(:, end);
Y4 = Y(:, end);
% Inner border coordinates
if inner_border == 1
    X_inner = [x(1, 1: L_half - gap_half), x(1, L_half + gap_half: end)];
    Y_inner = repmat(y_inner_lim, 1, length(X_inner));
else
    X_inner = [];
    Y_inner = [];
end
% Total border
X_border = [X1, X2, X3', X4', X_inner];
Y_border = [Y1, Y2, Y3', Y4', Y_inner];
end

% To create the people
function [X_people, Y_people] = make_people(people, X_border, Y_border, X_rshp, Y_rshp)
% This function creates the people that is going to be in movement afterwards inside
% the grid (with the exception of the borders).

% Input Arguments:
% people: integer. Number of people to generate.
% X_border: A vector. x coordinates of total border(unfolded).
% Y_border: A vector. y coordinates of total border(unfolded).
% X_rshp: A vector. x coordinates of the grid (unfolded).
% Y_rshp: A vector. y coordinates of the grid (unfolded).
% Output Arguments:
% X_people: A vector. x coordinates of the preople on the grid (unfolded).
% Y_people: A vector. y coordinates of the preople on the grid (unfolded).

% Initialization of the function and some checks:

% Find the indexes of coordinates that are in the border (all the grid)
% Annotate them with an 1. (1 means is in the border, otherwise 0)
length_border = length(X_border);
length_people = length(X_rshp);
no_border_lgl = ones(1, length_people);
for i = 1: length_people
    for j = 1: length_border
        condition_x = X_border(j);
        condition_y = Y_border(j);
        if X_rshp(i) == condition_x && Y_rshp(i) == condition_y
            no_border_lgl(i) = 0;
            break;
        end
    end
end
% Select the coordinates that are no in the border.
no_border = find(no_border_lgl == 1);
X_rshp_nb = X_rshp(no_border);
Y_rshp_nb = Y_rshp(no_border);

% Unsort them
perm_ind_X_rshp_nb = randperm(size(X_rshp_nb, 2));
perm_ind_Y_rshp_nb = randperm(size(Y_rshp_nb, 2));
% Select only a number of indexes from them equal to 'people'.
sel_ind_X_rshp_nb = perm_ind_X_rshp_nb(1: people);
sel_ind_Y_rshp_nb = perm_ind_Y_rshp_nb(1: people);
% Select the coordinates corresponding to 'people'.
X_people = X_rshp_nb(sel_ind_X_rshp_nb);
Y_people = Y_rshp_nb(sel_ind_Y_rshp_nb);
end

% To create health status
function [health_status, health_color, sick_time, recover_time] = make_health(X_people, Y_people, inner_border, Y_inner, people, sick_pc, max_sick_time)
% This function creates the initial health state for the people in the
% study: Than includes:
%   - Assigning a label corresponding to the patient's status:
%       1(non-infected), 2(infected), 3(recovered)
%   - Assigning a color to this statues:
%       blue(non-infected), magenta(infected), green(recovered)-
%   - Assigning the time since infection (for infected people). 
%       This is done just because we don't want inicial infected 
%       people becoming recovered exactly at the same time.
%
% Note: It is possible that infected people at time 1 recover. That means
% that they were infected in the past, before the study starts and they
% recover by the time that the study starts.

% Input Arguments:
% X_people: A vector. x coordinates of the preople on the grid (unfolded).
% Y_people: A vector. y coordinates of the preople on the grid (unfolded).
% inner_border: logical. 1 means that there is an inner borders,
%                        0 means that there isn't any.
% Y_inner: A vector with y coordinate for the inner border.
% people: integer. Number of people to generate.
% sick_pc: real. Percentage of sick people.
% max_sick_time: Maximum time being sick (in average) 
% Output Arguments:
% health_status: A vector with patients' status.
% health_color: A matrix. Color associated(rgb) associated to healt_status
% sick_time: A vector. Time since infection.
% recover_time: A vector. Time to recover from infection

% Inicialization of the function:
x = 1: length (X_people);
health_status = ones(1, length(x));
health_color = repmat([0 0 1], size(x,2), 1);

% If stament to differenciate cases: without(0) or with(1) inner borders
if inner_border == 0
    sick_people = round(people * sick_pc / 100);
    perm_ind_health = randperm(size(x, 2));
    sel_ind_sick = perm_ind_health(1: sick_people);
    health_status(sel_ind_sick) = 2;
    health_color(sel_ind_sick, :) = repmat([1 0 1], size(health_color(sel_ind_sick, :),1), 1);
elseif inner_border == 1
    aux_inner = find(Y_people < Y_inner(1));
    people_inner = length(aux_inner);
    sick_people = round(people_inner * sick_pc / 100);
    sel_ind_sick = aux_inner(1: sick_people);
    health_status(sel_ind_sick) = 2;
    health_color(sel_ind_sick, :) = repmat([1 0 1], size(health_color(sel_ind_sick, :),1), 1);
end

% 0 time means not infected yet
% betwen 1 and recover_time(i) means you are infected
% beyond max_sick_time recovered 
mean_sick_time_ini = 5;
std_sick_time_ini = 5;
sick_time = zeros(size(health_status));
prop_sick_time = abs(round(std_sick_time_ini .* randn(length(sel_ind_sick),1) + mean_sick_time_ini));
prop_sick_time(prop_sick_time == 0) = 1; % Infected people can't have sick_time < 0, in such a case -> 1
sick_time(sel_ind_sick) = prop_sick_time; 

% Randomize max_sick_time -> recover_sick_time
std_max_sick_time = 5;
prop_recover_time = abs(round(std_max_sick_time .* randn(length(sick_time),1) + max_sick_time));
prop_recover_time(prop_recover_time == 0) = 1; % Infected people can't have sick_time < 0, in such a case -> 1
recover_time = prop_recover_time; 
end

% To select the people who is moving
function [vX_people, vY_people, y_inner_max] = make_movement(X_people, inner_border, Y_inner, people, confinement, confined_pc)
% This function selects which people is moving in the grid. That simulates
% if people to which extend people remains quiet (i.e. at home) following
% a quarantine. Then initalizes the movement for the selected people.
%
% Input Arguments:
% X_people: A vector. x coordinates of the people on the grid (unfolded).
% people: integer. Number of people to generate.
% confinement: logical. 0 for no confinement, 1 for confinement.
% confined_pc: real. Percentage of sick people.
% Output Arguments:
% vX_people: A vector. Direction of x movement for the people on the grid (unfolded).
% vY_people: A vector. Direction of y movement for the people on the grid (unfolded).
% y_inner_max: Integer or NaN. Integer equal to Y_inner(1) if exists. If
% not set as NaN (this happens when there is no inner border).

% Inicializialization of the function:
% Select which people is going to move.
x = 1: length (X_people);
if confinement == 1
    moving_people = round(people * (100 - confined_pc) / 100);
    perm_ind_moving = randperm(size(x, 2)); 
    sel_ind_moving = perm_ind_moving(1: moving_people);
elseif confinement == 0
    sel_ind_moving = 1: length(x);
end
% Set the initial speed of people (8 directions).
speed_cond = randi(8,1,length(X_people));
vX_people = zeros(1, length(X_people));
vY_people = vX_people;
for i = sel_ind_moving
    switch(speed_cond(i))
        case(1)
            vX_people(i) = 1;
            vY_people(i) = 1;
        case(2)
            vX_people(i) = 1;
            vY_people(i) = 0;
        case(3)
            vX_people(i) = 1;
            vY_people(i) = -1;
        case(4)
            vX_people(i) = 0;
            vY_people(i) = -1;
        case(5)
            vX_people(i) = -1;
            vY_people(i) = -1;
        case(6)
            vX_people(i) = -1;
            vY_people(i) = 0;
        case(7)
            vX_people(i) = -1;
            vY_people(i) = 1;
        case(8)
            vX_people(i) = 0;
            vY_people(i) = 1;
    end
end
% Taking into account the inner frontier
if inner_border == 1
    y_inner_max = Y_inner(1);
elseif inner_border == 0
    y_inner_max = NaN;
end
end

% Update people position
function [X_people, Y_people, vX_people, vY_people] = update_position(X_people, Y_people, vX_people, vY_people, X_border, Y_border, X_inner, Y_inner)
% This function updates the position and direction of movement in the grid.

% Input Arguments:
% X_people: A vector. x coordinates of the people on the grid (unfolded).
% Y_people: A vector. y coordinates of the people on the grid (unfolded).
% vX_people: A vector. Direction of x movement for the people on the grid (unfolded).
% vY_people: A vector. Direction of y movement for the people on the grid (unfolded).
% X_border: A vector. x coordinates of total border(unfolded).
% Y_border: A vector. y coordinates of total border(unfolded).
% X_inner: A vector with x coordinate for the inner border.
% Y_inner: A vector with y coordinate for the inner border.
% Output Arguments:
% X_people: A vector. New x coordinates of the people on the grid (unfolded).
% Y_people: A vector. New y coordinates of the people on the grid (unfolded).
% vX_people: A vector. New Direction of x movement for the people on the grid (unfolded).
% vY_people: A vector. New Direction of y movement for the people on the grid (unfolded).

for i = 1: length(X_people)
    % Inicialization of the function and some checks:
    % Check if we have inner border or not
    if isempty(X_inner) == 1
        cond_inner_1 = 0;
        cond_inner_2 = 0;
    else
        cond_inner_1 = length(intersect(X_people(i), X_inner));
        cond_inner_2 = Y_people(i) == Y_inner(1);
    end
    % Outer border
    if X_people(i) == min(X_border) || X_people(i) == max(X_border)
        % the person i is a horizontal border (external)
        if Y_people(i) == min(Y_border) || Y_people(i) == max(Y_border)
            % the person i is also in a vertical border (external)
            % (so the person is in a corner)
            vX_people(i) = - vX_people(i);
            vY_people(i) = - vY_people(i);
        elseif cond_inner_1 == 1 && cond_inner_2 == 1
            % The person i is in a point of the grid that belongs
            % to the inner and also to the external border
            vX_people(i) = - vX_people(i);
            vY_people(i) = - vY_people(i);
        else
            % The person i is in a point of the external border X ---
            vX_people(i) = - vX_people(i);
            vY_people(i) = + vY_people(i);
        end
    elseif Y_people(i) == min(Y_border) || Y_people(i) == max(Y_border)
        % The person i is in a point of the external border Y |
        vX_people(i) = + vX_people(i);
        vY_people(i) = - vY_people(i);
    elseif cond_inner_1 == 1 && cond_inner_2 == 1
         % The person i is in a point only from the inner border
        vX_people(i) = + vX_people(i);
        vY_people(i) = - vY_people(i);
    end
end
% Update the position
X_people = X_people + vX_people;
Y_people = Y_people + vY_people;
end

% Update health (non-infected to infected)
function  [health_status, health_color, infection_matrix, sick_time] = update_sick(X_people, Y_people, y_inner_max, health_status, health_color, infection_matrix, sick_time)
% This function updates the health status, the health_color and starts the
% sick_time counter  when someone is infected.

% Input Arguments:
% X_people: A vector. x coordinates of the people on the grid (unfolded).
% Y_people: A vector. y coordinates of the people on the grid (unfolded).
% y_inner_max: Integer or NaN. Integer equal to Y_inner(1) if exists. If
% it's an NaN (this happens when there is no inner border).
% health_status: A vector with patients' status.
% health_color: A matrix. Color associated(rgb) associated to healt_status.
% infection_matrix: Matrix. 
                    % Who(row) infected whom (columns)
                    % 1 infected, 0 non-infected.
% sick_time: A vector. Time since infection.
% Output Arguments:
% health_status: New vector with patients' status.
% health_color: A matrix. New color associated(rgb) associated to healt_status.
% infection_matrix: Matrix updated. 
                    % Who(row) infected whom (columns)
                    % 1 infected, 0 non-infected
% sick_time: A vector. Set to 1 to ner sick people.

% In orther to save time, the loop that controls the what people is in
% contact can changes the number of elements to iterate in each operation.
% We do this combining an infinte while loop with an abort condition with a
% for loop.
ind_people = 1: length(X_people);
    l = 0;
    while 1
        l = l + 1;
        old_length = length(ind_people);
        for i = ind_people
            coincidences_zeros = zeros(1, length(ind_people));
            k = 0;
            for j = ind_people
                if i ~= j % we don't count twice the same person
                    if Y_people(i) ~= y_inner_max && Y_people(j) ~= y_inner_max % people in the border can't infect.
                        if X_people(i) == X_people(j) && Y_people(i) == Y_people(j)
                            if health_status(i) == 3 || health_status(j) == 3
                            elseif health_status(i) == 2 && health_status(j) == 1 % i infects j
                                infection_matrix(i ,j) = 1; % infection ij
                                health_status(j) = 2; % change status to sick
                                health_color(j, :) = [1 0 1]; % change to color magenta
                                sick_time(j) = 1; % start the sick time counter
                                k = k + 1;
                                coincidences_zeros(k) = j;
                            elseif health_status(i) == 1 && health_status(j) == 2 %j infects i
                                infection_matrix(j ,i) = 1; % infection ij
                                health_status(i) = 2; %change status to sick
                                health_color(i, :) = [1 0 1]; % change to color magenta
                                sick_time(i) = 1; % start the sick time counter
                                k = k + 1;
                                coincidences_zeros(k) = j; 
                            end
                        end
                    end
                end
                
            end
            % inner for loop break condition:
            coincidences = find(coincidences_zeros ~= 0);
            aux = setdiff(ind_people, [i coincidences]); 
            if length(aux) < length(ind_people) - 1
                ind_people = aux;
                break;
            end
        end
        % inner outer while loop break condition:
        new_length = length(ind_people);
        if new_length == old_length
            break;
        end
    end
end

% Update health (infected to recovered)
function [health_status, health_color, sick_time] = update_recover(health_status, health_color, sick_time, max_sick_time, plot_epidemic_curves_check)
% This function updates the health status, the health_color and starts the
% sick_time counter  when someone is recovered.

% Input Arguments:

% health_status: A vector with patients' status.
% health_color: A matrix. Color associated(rgb) associated to healt_status.
% sick_time: A vector. Time since infection.
% max_sick_time: An integer. Maximum time infected.
% Output Arguments:
% health_status: New vector with patients' status.
% health_color: A matrix. New color associated(rgb) associated to healt_status.
% sick_time: A vector. Increment of 1 unity for sick people.

% Randomize max_sick_time -> recover_sick_time
std_max_sick_time = 5;
prop_recover_sick_time = abs(round(std_max_sick_time .* randn(length(sick_time),1) + max_sick_time));
prop_recover_sick_time(prop_recover_sick_time == 0) = 1; % Infected people can't have sick_time < 0, in such a case -> 1
recover_sick_time = prop_recover_sick_time; 

for i = 1: length(sick_time)
    if sick_time(i) > 0 && sick_time(i) < recover_sick_time(i)
        sick_time(i) = sick_time(i) + 1; % increment the time a person is sick
    elseif sick_time(i) >= recover_sick_time(i)% it may happen that infected people is recovered at first iteration
        health_status(i) = 3;  %change status to recovered
        health_color(i, :) = [0 1 0]; % change to color green
    end
end
end

% Plot the epidemic curves
function plot_epidemic_curves(do_plot, non_infected, infected, recovered, inner_border, confinement)

% This function plots the result of the simulation (the epidemic curves)
% with an appropiate title

% Initial check:
% Should we plot?
if do_plot == 0
    return
end

total_pop = non_infected(1) + infected(1) + recovered(1);
% Select the case_study to set the title.
if inner_border == 0 && confinement == 0
    case_study = 'Epidemic Curves: No measures';
elseif inner_border == 0 && confinement == 1
    case_study = 'Epidemic Curves: Isolation';
elseif inner_border == 1 && confinement == 0
    case_study = 'Epidemic Curves: Contention';
elseif inner_border == 1 && confinement == 1
    case_study = 'Epidemic Curves: Isolation and Contention';
end
% Create the plot
figure()
plot(100 * non_infected / total_pop,'b','LineWidth',2)
hold on
plot(100* infected / total_pop,'m','LineWidth',2)
hold on
plot(100 * recovered / total_pop,'g','LineWidth',2)
grid on
xlabel('Time (a.u.)')
ylabel('Population')
title(case_study)
legend('Non-Infected','Infected','Recovered','Location','NorthWest')
end