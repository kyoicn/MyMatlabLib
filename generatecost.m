%% cost distribution

function cost = generatecost(m, n, cost_distribution)
cost = zeros(m, n);

%% cost distribution method
% 1: station-random
% 2: station-uniform
% 3: station-fixed
% 4: random cost
% 5: fixed

if (strcmp(cost_distribution, 'stationr'))
    cost_code = 1;
elseif (strcmp(cost_distribution, 'stationu'))
    cost_code = 2;
elseif (strcmp(cost_distribution, 'stationf'))
    cost_code = 3;
elseif (strcmp(cost_distribution, 'random'))
    cost_code = 4;
elseif (strcmp(cost_distribution, 'fixed'))
    cost_code = 5;
elseif (strcmp(cost_distribution, 'uniform'))
    cost_code = 6;
end
    

%% station distribution
% 1: fixed
% 2: random
% 3: uniform + random

% station_distribution = 3;

%% stations
% assuming we have several base station in the area, according to the
% distance between one position and its nearset base station, the cost for
% sensing differs. Here we simply use an alpha-fading model.

ns = 36;
stations = zeros(ns, 2);

if (cost_code == 1) % station-random
    for i = 1 : ns
        stations(i, 1) = ceil(rand() * m);
        stations(i, 2) = ceil(rand() * n);
    end
    
elseif (cost_code == 2) % station-uniform
    row = floor(sqrt(ns));
    rem = ns - row^2;
    gapr = floor(m / (row - 1)) - 1;
    gapc = floor(n / (row - 1)) - 1;
    for i = 1 : row^2
        stations(i, 1) = 1 + gapr * floor((i - 1) / row);
        stations(i, 2) = 1 + gapc * mod(i - 1, row);
    end
    for i = row^2 + 1 : ns
        stations(i, 1) = ceil(rand() * m);
        stations(i, 2) = ceil(rand() * n);
    end    
elseif (cost_code == 3) % station-fixed
    stations = [
     1     1
     1    32
     1    63
    32     1
    32    32
    32    63
    63     1
    63    32
    63    63
    62    27
    16    43
    45    54
    60    55
    34    25
     5    21];    
end

%% Cost

if (cost_code <= 3) % station-based
    
    alpha = 2;
    basicc = 1;
    penalty = 1;
    cost = ones(m, n) * basicc;
    
    for i = 1: m
        for j = 1 : n
            mi = -1;
            for ii = 1 : ns
                if (mi == -1 || mi > sqrt((i - stations(ii, 1))^2 + (j - stations(ii, 2))^2))
                    mi = sqrt((i - stations(ii, 1))^2 + (j - stations(ii, 2))^2);
                end
            end
            cost(i, j) = cost(i, j) + penalty * mi^alpha;
        end
    end    
 
elseif (cost_code == 4) % random cost
    cost = rand(m, n);
    
elseif (cost_code == 5) % fixed
    % as fixed
elseif (cost_code == 6) % uniform
    cost = ones(m, n);
    
end