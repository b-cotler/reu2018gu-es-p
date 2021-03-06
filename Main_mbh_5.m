
%% Set Variables %%

rng('shuffle'); % random seed for random number generator

% declare main global variables for program
total_pop_N = 500; % size of population for all age classes

number_generations = 2000; % number of generations

burn_in_gens = 130; % number of generations for burn in of population growth

lineage_count = 6; % number of lineages to sample to determine time to MRCA 

iterations = 500; % number of iterations of sampling from population

fprintf('----------------------------------------------------\n');
fprintf('Simulation of time to MRCA in an age-structured coalescent\n\n');

t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % get date and time
date_string = datestr(t); % convert date and time to string
fprintf('%s\n\n',date_string);

fprintf('Initial total population size: %g\n\n', total_pop_N);

%% open file with life table, get Leslie matrix for population growth %%

% provide path name to life table file
file_path_name = 'Atlantic_cod_life_table.xlsx';
fprintf('Life table file: %s\n\n',file_path_name);

% scaling factor to adjust population growth rate; 
% set = 1/lambda for same lambda as in the life table
% when scale = 1.0 there is no population growth (lambda = 1.0)
% when scale > 1.0 there is negative population growth rate (lambda < 1.0)
% when scale < 1.0 there is positive population growth rate (lambda > 1.0)
scale = 1;

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
[age_dist_m, burn_in_m,total_population_v,steady_state_total, steady_state_vector] = create_age_dist_m(number_generations, population_0, leslie_matrix, burn_in_gens,mod_lambda); 

%% Sample lineages %%

% sample pairs of random age lineages
age_i = -1; % set lineages sampled in present to random ages

no_mrca_random = 0; % counter for number of times no MRCA was found

mrca_random = zeros(1,iterations); % allocate space for results

mutation_m_random = zeros(lineage_count+1,iterations); %create a matrix to record mutation opportunities

random_mutation_avg = zeros(1,lineage_count); %create a vector to measure weighted average of mutation opportunities

first_coal_event_random = zeros(iterations,3); %create a matrix to keep track of the first coalescent event that occurs (to be used when K>2)
first_coal_event_random_data = zeros(iterations,1);
no_first_coal_random = 0;

for iter=1:iterations

    initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % samples lineage from all lineages in the present. 

    genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!

    genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

    [mrca,complete_genealogy,coal_events,age_zero_counter, total_multiple_mergers] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);

    if mrca == number_generations
        no_mrca_random = no_mrca_random + 1; % increment counter
    else
        mrca_random(1,iter) = mrca;
        
    end
    
    first_coal_event_random(iter,:) = coal_events(1,:); %take the first coalescent event that occurs from the coalescent events matrix
    if isequal(first_coal_event_random(iter,3), -1)
        no_first_coal_random = no_first_coal_random+1;
    else
        first_coal_event_random_data(iter) = first_coal_event_zero(iter,3);
    end
    
    % accumulate mutation opportunities and weight by mrca
    mutation_m_random(1:lineage_count,iter) = age_zero_counter; 
    mutation_m_random(end, iter) = mrca;
    for i = 1:lineage_count
    random_mutation_avg(i) = random_mutation_avg(i) + age_zero_counter(i)/mrca;
    end
    

end % for iter

    % compute mean time to MRCA without iterations that did not experience
    % coalescence
    num_non_zero_sims_random = iterations - no_mrca_random;
    sum_mrca_random = sum(mrca_random);
    mean_random = sum_mrca_random/num_non_zero_sims_random;
    median_random = median(mrca_random(1,1:num_non_zero_sims_random));
    
    % establish final weighted mutations both cumulative and per lineage
    cumulative_random_mutation_avg = 0;
    for i = 1:lineage_count
        random_mutation_avg(i) = random_mutation_avg(i)/iterations;
        cumulative_random_mutation_avg = cumulative_random_mutation_avg + random_mutation_avg(i);
    end
    cumulative_random_mutation_avg = cumulative_random_mutation_avg/lineage_count;
    
% sample pairs of age zero lineages
age_i = 0; % set lineages sampled in present to age zero

no_mrca_zero = 0; % counter for number of times no MRCA was found

mrca_zero = zeros(1,iterations); % allocate space for results

%Creates a mutaton matrix to track possible mutation events

%The first row contains the number of age zero individuals in the first lineage at time iter and 
%the second row contains the number of age zer individuals in the second lineage at time iter. 
%The third row contains the time to mrca because the number of age zero individuals will be smaller if an mrca is reached early in the lineage trace. "

mutation_m_zero = zeros(lineage_count+1, iterations);
zero_mutation_avg = zeros(1,lineage_count); %create a vector to measure weighted average of mutation opportunities
first_coal_event_zero = zeros(iterations,3); %create a matrix to keep track of the first coalescent event that occurs (to be used when K>2)
first_coal_event_zero_data = zeros(iterations,1);
no_first_coal_zero = 0;

for iter=1:iterations

    initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % samples lineage from all lineages in the present. 

    genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!

    genealogy_m(end,:,1) = initial_values(1,:);
    genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

    [mrca,complete_genealogy,coal_events, age_zero_counter, total_multiple_mergers] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);

    if mrca == number_generations
        no_mrca_zero = no_mrca_zero + 1; % increment counter
    else
        mrca_zero(1,iter) = mrca;
    end
    
    first_coal_event_zero(iter,:) = coal_events(1,:); %take the first coalescent event that occurs from the coal_events matrix
    if isequal(first_coal_event_zero(iter,3), -1)
        no_first_coal_zero = no_first_coal_zero+1;
    else
        first_coal_event_zero_data(iter) = first_coal_event_zero(iter,3);
    end
        
        
    
    mutation_m_zero(1:lineage_count,iter) = age_zero_counter;
    mutation_m_zero(end, iter) = mrca;
    for i = 1:lineage_count
    zero_mutation_avg(i) = zero_mutation_avg(i) + age_zero_counter(i)/mrca;
    end
    
end % for iter

    % compute mean time to MRCA without iterations that did not experience
    % coalescence
    num_non_zero_sims_zero = iterations - no_mrca_zero;
    sum_mrca_zero = sum(mrca_zero);
    mean_zero = sum_mrca_zero/num_non_zero_sims_zero;
    median_zero = median(mrca_zero(1,1:num_non_zero_sims_zero));

    % establish final weighted mutations both cumulative and per lineage
    cumulative_zero_mutation_avg = 0;
    for i = 1:lineage_count
        zero_mutation_avg(i) = zero_mutation_avg(i)/iterations;
        cumulative_zero_mutation_avg = cumulative_zero_mutation_avg + zero_mutation_avg(i);
    end
    cumulative_zero_mutation_avg = cumulative_zero_mutation_avg/lineage_count;
    
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
    
%    savefig('histogram_growing_Atlantic_cod_highCvf.fig');
    
    figure;
    
    hold on;
    subplot(1,2,1);
    boxplot(mrca_zero(1,1:num_non_zero_sims_zero), 'Labels',{'age zero lineage pairs'});
    
    subplot(1,2,2);
    boxplot(mrca_random(1,1:num_non_zero_sims_random), 'Labels',{'random age lineage pairs'});
    hold off;

%    savefig('boxplot_growing_Atlantic_cod_highCvf.fig');

%BREAK

%     figure;
% 
%     hold on;
%     subplot(2,1,1);
%     hist(first_coal_event_random_data(1:(iterations-no_first_coal_random),1)); 
%     xlabel('time to MRCA - random lineage pairs')
%     ylabel('Count')
% 
%     subplot(2,1,2);
%     hist(first_coal_event_zero_data(1:(iterations-no_first_coal_zero),1)); 
%     xlabel('time to MRCA - age zero lineage pairs')
%     ylabel('Count')
% 
%     suptitle('Distributions of coalescence times');
%     hold off;
%     
% %    savefig('histogram_growing_Atlantic_cod_highCvf.fig');
%     
%     figure;
%     hold on;
%     subplot(1,2,1);
%     boxplot(first_coal_event_zero_data(1:(iterations-no_first_coal_zero),1), 'Labels',{'age zero lineage pairs'});
%     
%     subplot(1,2,2);
%     boxplot(first_coal_event_random_data(1:(iterations-no_first_coal_random),1), 'Labels',{'random age lineage pairs'});
%     hold off;
% 
% %    savefig('boxplot_growing_Atlantic_cod_highCvf.fig');

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

output_filename = "/Users/BrettCotler/Desktop/Output_Data/atlantic_cod_k6constant_cauchy.mat";
save(output_filename);
% 
% 
