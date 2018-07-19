
%% Set Variables %%

rng('shuffle'); % random seed for random number generator

% declare main glaobal variables for program
total_pop_N = 500; % size of population for all age classes

number_generations = 2400; % number of generations

burn_in_gens = 102; % number of generations for burn in of population growth

lineage_count = 2; % number of lineages to sample to determine time to MRCA 

iterations = 500; % number of iterations of sampling from population

fprintf('----------------------------------------------------\n');
fprintf('Simulation of time to MRCA in an age-structured coalescent\n\n');

t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % get date and time
date_string = datestr(t); % convert date and time to string
fprintf('%s\n\n',date_string);

fprintf('Initial total population size: %g\n\n', total_pop_N);

%% open file with life table, get Leslie matrix for population growth

% provide path name to life table file
file_path_name = 'sage_grouse_life_table.xlsx';
fprintf('Life table file: %s\n\n',file_path_name);

% scaling factor to adjust population growth rate; 
% set = 1/lambda for same lambda as in the life table
% when scale = 1.0 there is no population growth (lambda = 1.0)
% when scale > 1.0 there is negative population growth rate (lambda < 1.0)
% when scale < 1.0 there is positive population growth rate (lambda > 1.0)
scale = 1.0;

[life_table_m,leslie_matrix,age_classes,orig_lambda,mod_lambda,CV_fecundity,G,alpha,AL,mean_fecundity] = life_to_leslie_5(file_path_name,scale);
fprintf('life table has a population growth rate of %f\n', orig_lambda);
fprintf('(set scale to 1/lambda or = %f for original rate)\n\n', 1/orig_lambda);

fprintf('after rescaling, the population growth rate is %f\n\n', mod_lambda);

fprintf('number of age classes: %g\n\n', age_classes);

fprintf('coefficient of variation in age-specific fecundity: %f\n\n', CV_fecundity);
%% Establish Age Distribution Matrix (age_dist_m)

time = 1:number_generations; % create a time vector using the number of generations
 
uniform_age_cohort_N = round(total_pop_N/(age_classes)); % initial uniform size for all age cohorts
 
population_0 = zeros(1,age_classes);
 
population_0(1:(age_classes)) = uniform_age_cohort_N;
 
total_population_0 = sum(population_0);
 
 
if total_population_0 - total_pop_N > 0
    
    population_0(1) = population_0(1) + (total_population_0 - total_pop_N); 
    
elseif total_population_0 - total_pop_N < 0
    
     population_0(1) =  population_0(1) + (total_pop_N - total_population_0);
     
end

% create the demographic matrix of population size for each age class over time
age_dist_m = create_age_dist_m(number_generations, population_0, leslie_matrix, burn_in_gens); 

%% Sample lineages %%

% sample pairs of random age lineages
age_i = -1; % set lineages sampled in present to random ages

no_mrca_random = 0; % counter for number of times no MRCA was found

mrca_random = zeros(1,iterations); % allocate space for results

for iter=1:iterations

    initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % samples lineage from all lineages in the present. 

    genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!

    genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

    [mrca,complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);

    if mrca == number_generations
        no_mrca_random = no_mrca_random + 1; % increment counter
    else
        mrca_random(1,iter) = mrca;
        
    end
end % for iter

    % compute mean time to MRCA without iterations that did not experience
    % coalescence
    num_non_zero_sims_random = iterations - no_mrca_random;
    sum_mrca_random = sum(mrca_random);
    mean_random = sum_mrca_random/num_non_zero_sims_random;
    median_random = median(mrca_random(1,1:num_non_zero_sims_random));


% sample pairs of age zero lineages
age_i = 0; % set lineages sampled in present to age zero

no_mrca_zero = 0; % counter for number of times no MRCA was found

mrca_zero = zeros(1,iterations); % allocate space for results

%Creates a mutaton matrix to track possible mutation events

%The first row contains the number of age zero individuals in the first lineage at time iter and 
%the second row contains the number of age zer individuals in the second lineage at time iter. 
%The third row contains the time to mrca because the number of age zero individuals will be smaller if an mrca is reached early in the lineage trace. "

mutation_m = zeros(lineage_count+1, iterations);

for iter=1:iterations

    initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % samples lineage from all lineages in the present. 

    genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!

    genealogy_m(end,:,1) = initial_values(1,:);
    genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

    [mrca,complete_genealogy,coal_events, age_zero_counter] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);

    if mrca == number_generations
        no_mrca_zero = no_mrca_zero + 1; % increment counter
    else
        mrca_zero(1,iter) = mrca;
    end
    
    mutation_m(1:lineage_count,iter) = age_zero_counter;
    mutation_m(end, iter) = mrca;
    
end % for iter

    % compute mean time to MRCA without iterations that did not experience
    % coalescence
    num_non_zero_sims_zero = iterations - no_mrca_zero;
    sum_mrca_zero = sum(mrca_zero);
    mean_zero = sum_mrca_zero/num_non_zero_sims_zero;
    median_zero = median(mrca_zero(1,1:num_non_zero_sims_zero));


    figure;

    hold on;
    subplot(2,1,1);
    hist(mrca_random(1,1:num_non_zero_sims_random)); 
    xlabel('time to MRCA - random lineage pairs')
    ylabel('Count')

    subplot(2,1,2);
    hist(mrca_zero(1,1:num_non_zero_sims_zero)); 
    xlabel('time to MRCA - age zero lineage pairs')
    ylabel('Count')

    suptitle('Distributions of coalescence times');
    hold off;
    
%savefig('histogram_constant_Atlantic_cod_extreme_shortAL.fig');
    
    figure;
    hold on;
    subplot(1,2,1);
    boxplot(mrca_zero(1,1:num_non_zero_sims_zero), 'Labels',{'age zero lineage pairs'});
    
    subplot(1,2,2);
    boxplot(mrca_random(1,1:num_non_zero_sims_random), 'Labels',{'random age lineage pairs'});
    hold off;

%savefig('boxplot_constant_Atlantic_cod_extreme_shortAL.fig');

    fprintf('Summary of simulation:\n\n');
        
    fprintf('Number of iterations was %i\n', iterations);
    fprintf('Maximum number of generations was %i\n', number_generations);
    fprintf('Total population size was %i\n', total_pop_N);
    
    fprintf('\n');
    
    fprintf('zero age lineages coalescence times:\n');
    fprintf('average: %g \n', mean_zero);
    fprintf('median: %g \n', median_zero);
    fprintf('number of simulations that reached coalescence: %g of %g \n', num_non_zero_sims_zero, iterations);
    
    fprintf('\n');
   
    fprintf('random age lineages coalescence times:\n');
    fprintf('average: %g \n', mean_random);
    fprintf('median: %g \n', median_random);
    fprintf('number of simulations that reached coalescence: %g of %g \n', num_non_zero_sims_random, iterations);

fprintf('----------------------------------------------------\n');

%Save entire workspace as a ".mat" file

% output_filename = "/Users/BrettCotler/Desktop/Output_Data/atlantic_cod/constant_Atlantic_cod_extreme_shortAL.mat";
% save(output_filename);


function [parent] = sample_lineage(age)


    % sample one lineage given an age class and a generation
    age = parent_age_class; % age of lineage to find
    
    generation = numgen - 1; % generation of lineage to find
    
    found_one = 0;
    while found_one == 0

        index = randi(total_pop_sizes(generation,1)); % get random integer 1:total population size at time numgen
        
        if pop_in_present(1,index) == age % lineage must be given age 
            parent = index;
            found_one = 1;
        else
            found_one = 0; % keep looking for age zero individual
        end
    end % while

    disp(parent);
    
end

