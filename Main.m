%% Short Hand %%

% "_0" label is time 0
% "_m" label indicates a matrix
% "_v" label indicates a vector
% "dist" = distribution

%% User Inputs %%

% population_0 = input('Input an initial population vector and press enter: ');
% number_generations = input('Input a number of generations and press enter: ');
% leslie_matrix = input('Input a Leslie Matrix and press enter: ');
% lineage_count = input('Input a number of lineages to track and press enter: ');

%% Example Inputs %%

%EXAMPLE 1
% population_0 = [2, 2, 2]; %example initial population vector
% number_generations = 10; %example number of generations
% leslie_matrix = [0 1 1.1; 0.6 0 0; 0 0.5 0]; %example leslie matrix
% lineage_count = 2; %example number of lineages to track

%EXAMPLE 2
population_0 = [30, 25, 15, 30]; %example initial population vector
number_generations = 200; %example number of generations
leslie_matrix = [0 1 1.1 1.2; 0.6 0 0 0; 0 0.5 0 0; 0 0 0.4 0]; %sample leslie matrix
lineage_count = 2; %example number of lineages to track 
%life_table = [0 0.6 0; 1 0.5 1; 2 0.25 1.1; 3 0 1.2];

lambda = eig(leslie_matrix);

if lambda(1)> 0
    scaling = lambda(1);
else
    scaling = lambda(2);
end

leslie_matrix = leslie_matrix./scaling;

%EXAMPLE 3
% population_0 = [56 72 16 84 23];
% number_generations = 100; 
% leslie_matrix = [0 0.5 0.75 1 1.25; 0.9 0 0 0 0; 0 0.75 0 0 0; 0 0 0.6 0 0; 0 0 0 0.4 0];
% lineage_count = 2;

%EXAMPLE 4
% population_0 = [506 823 234 348 294 1844];
% number_generations = 50;
% leslie_matrix = [0 0.1 0.1 0.1 0.1 0; 0.95 0 0 0 0 0; 0 0.92 0 0 0 0; 0 0 0.87 0 0 0; 0 0 0 0.56 0 0; 0 0 0 0 0.3 0];
% lineage_count = 2; 

%% Establish Age Distribution Matrix (age_dist_m) %%

time = 1:number_generations; %creates a time vector using the set number of generations

age_dist_m = zeros(length(population_0), length(time)); %creates a matrix with age distributions at each time step
age_dist_m(:,1) = population_0; %sets the initial age distributions to the initial population values in population_0

total_population_0 = sum(population_0); %calculates the total population at time = 0 by adding the values of population_0
total_population_v = zeros(1, length(time)); %initializes a vector of total population values
total_population_v(1) = total_population_0; %sets the first entry of the total population vector equal to the initial total population

for i = 2:length(time)
    age_dist_m(:,i) = round(leslie_matrix*age_dist_m(:,i-1)); %applies the leslie matrix to the previous age distribution for each time step
    
    total_population_v(i) = sum(age_dist_m(:,i)); %adds the total population at time step i to the total population vector
end

%% Establish a Matrix of Individuals (individuals_m) %%

%max_population = max(total_population_v); %sets the column dimension of the individuals matrix by determining the maximum population value from the total population vector

% individuals_m = -1*ones(number_generations,max_population); %creates a matrix with number_generations rows and max_population columns where every entry is -1
% 
% for i = 1:length(time)
%     individuals_m(i,1:age_dist_m(1,i)) = 0; %boundary case, set the first n individuals to be age zero (using the first index of the age_dist_m matrix)
%     
%     for j = 2:length(population_0)
%         index = 1+sum(age_dist_m(1:j-1,i)); %determine the number of individuals already assigned to an age, "+1" (MATLAB indices are inclusive)
%         
%         individuals_m(i,index:(index-1+age_dist_m(j,i))) = j-1; %assign individuals of the next age group to the newest open spaces in the individuals_m matrix
%     end
% end

%% Extract the Terminal Population and Choose Number of Individuals to Track  %%

%terminal_population = extract_terminal_population2(age_dist_m); %returns the portion of the last row which contains individuals

age_i = -1; %Choose an within the possible range to choose two individual in the same age class. Enter -1 for two random individuals. 
initial_values = terminal_indices(lineage_count, age_dist_m, age_i); %function which returns two random indices from the final row of individuals. 

genealogy_m = -1*ones(number_generations, lineage_count, 2); %initialize the 3-D genealogy matrix
genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); %set the front row to the indices and the back row to the ages specified in the initial_values matrix

%% Track the lineages to an MRCA %%

[mrca, complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);
if isequal(mrca, number_generations)
    disp('no mrca')
else
    mrca
end
complete_genealogy
coal_events


%% Track the lineages to an MRCA %%



%disp(terminal_population);
% ez we'll start with just two hard-coded pairs to find the mrca of.

% lineage_a_current_age = 0;
% lineage_b_current_age = 0;
% 
% 
% mrca = calculate_mrca(lineage_a_current_age, lineage_b_current_age, ...
%     terminal_population, life_table, number_generations, age_dist_m, individuals_m);
% 
% disp('calculate_mrca returns value : ')
% disp(mrca)

    
    

