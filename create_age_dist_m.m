function [age_dist_m] = create_age_dist_m(number_generations, population_0, leslie_matrix,burn_in_gens)
%Inputs: 
%1. number_generations - The preset number of generations that we are testing the data on
%2. population_0 - a vector with the initial population of individuals for
%each age class
%3. leslie_matrix - a leslie matrix corresponding to the input life table,
%generated using the life_to_leslie function
%4. burn-in-gens - the number of generations required for the population to
%stabilize

% This function generates a burn-in matrix to allow the population to
% stabilize, and an age distribution matrix which contains the distribution
% of individuals among the age classes over the preset number of
% generations. The age distribution matrix begins with the final column of
% the burn-in matrix. This function also displays warnings if the
% population size of an age class grows too large, and a warning of when
% the population stabilizes to ensure the number of burn-in generations is
% correctly set. 


%Outputs:
%1. age_dist_m - The only output is the demographic matrix/age distribution matrix



time = 1:number_generations; %creates a time vector using the set number of generations

total_population_0 = sum(population_0); %calculates the total population at time = 0 by adding the values of population_0
total_population_v = zeros(1, length(time)); %initializes a vector of total population values
total_population_v(1) = total_population_0; %sets the first entry of the total population vector equal to the initial total population

%% Burn-In %%
warning = 0;
burn_in_m = zeros(length(population_0),burn_in_gens);
burn_in_m(:,1) = population_0;
for i = 2:burn_in_gens
    burn_in_m(:,i) = round(leslie_matrix*burn_in_m(:,i-1));
    total_population_v(i) = sum(burn_in_m(:,i));
    
    if isequal(burn_in_m(:,i),burn_in_m(:,i-1)) && isequal(warning,0)
        fprintf("The Population reaches a steady state at burn-in t = %d \n", i);
        burn_in = false;
        warning = 1;
    elseif isequal(i,burn_in_gens) && isequal(warning,0)
        burn_in = false;
        fprintf("The Population has not reached steady state after %d burn-in generations \n", burn_in_gens);
        warning = 1;
    end
size_check_v = burn_in_m(:,i) > (2^50)*ones(size(burn_in_m,1),1); 
    if isequal(size_check_v,zeros(size(burn_in_m,1),1) == 0) && (isequal(warning,0))
        fprintf("The Population is too large at burn-in t = %d \n",i); 
        burn_in = false;
        warning = 1;
    end
end

%% Matrix Creation %%
if burn_in == false

age_dist_m = zeros(length(population_0), length(time)); %creates a matrix with age distributions at each time step
age_dist_m(:,1) = burn_in_m(:,end); %sets the initial age distributions to the initial population values in population_0
warning2 = 0;
for i = 2:length(time)
    age_dist_m(:,i) = round(leslie_matrix*age_dist_m(:,i-1)); %applies the leslie matrix to the previous age distribution for each time step
    
    total_population_v(i) = sum(age_dist_m(:,i)); %adds the total population at time step i to the total population vector
size_check_v = age_dist_m(:,i) > (2^50)*ones(size(age_dist_m,1),1); 
    if isequal(size_check_v,zeros(size(age_dist_m,1),1) == 0) && (isequal(warning2,0))
        fprintf("The Population is too large at t = %d \n",i); 
        warning2 = 1;
    end
    if isequal(age_dist_m(:,i),age_dist_m(:,i-1)) && isequal(warning2,0)
        age_dist_m = age_dist_m(:,1:i);
        fprintf("The Population reaches a steady state at t = %d \n", i);
        warning2 = 1;
    end
end
end