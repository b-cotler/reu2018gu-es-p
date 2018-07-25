function [indices] = terminal_indices_2(lineage_count, age_dist_m, age)
%Inputs:
%1. lineage_count - the preset number of lineages that the user wants to
%sample
%2. age_dist_m - the age distribution/demographic matrix containing the
%population of each age class over the number of generations
%3. age -vector, index should be -1 if #lineage_count random individuals are to be
%selected, or a given age >= 0 if #lineage_count individuals from that
%specific age class are to be sampled. for example age = [1 0] will sample
%one age 1 and one age 0 individual. age = [-1 -1 -1] will sample three
%random age individuals. age = [0 0] will sample two age 0 individuals.

%This function randomly selects indices in the final generation of the
%demographic matrix to begin the lineage trace. If the age input is -1 it
%randomly selects #lineage_count individuals randomly (without
%replacement). If the age input is 0 or higher it randomly selects
%individuals from that given age class

%Outputs:
%1. indices - indices is a 2xlineage_count matrix which contains both the
%index of each selected individual and its age. The first row contains the
%indices and the second row contains the ages.

indices = zeros(2, lineage_count); %initializes the indices matrix which will fill the top surface of the genealogy matrix
[age_classes,generations] = size(age_dist_m); %measures size of age_dist_m
options_lower = [1]; %sets a lower bound vector for each age class
options_upper = [age_dist_m(1,end)]; %sets an upper bound vector for each age class
for a = 2:size(age_dist_m)
    options_lower = [options_lower sum(age_dist_m(1:a-1,end))+1]; %construct a vector of lower bounds
    options_upper = [options_upper sum(age_dist_m(1:a,end))]; %construct a vector of upper bounds
end
scaling = sum(age_dist_m(1:age_classes,generations)); %scale factor for age class sizes
probabilities = age_dist_m(1:age_classes,generations)./scaling; %scales the age class sizes as probabilities
stored_age = age; %stores the age for each iteration of the loop

%% Specified Same Age %%
for k = 1:lineage_count
    if isequal(age(k),-1)
        age_index = rand; %chooses random number between zero and one
        for a = 1:age_classes
            if age_index <= sum(probabilities(1:a)) %check if the random number is less than the values of the previous probabilities combined
                new_age = a-1; %if the random number falls in the corresponding range of an age class, assign the individual to that age class
                break %the loop is broken when an age is assigned
            end
        end
    else
    new_age = age(k);
    end
    indices(1,k) = options_lower(new_age+1) + round((options_upper(new_age+1)-options_lower(new_age+1))*rand);
    indices(2,k) = new_age;
    if k>1
        for r = 1:k-1
            while isequal(indices(1,k),indices(1,r))
                age_index = rand; %chooses random number between zero and one
                for a = 1:age_classes
                    if age_index <= sum(probabilities(1:a)) %check if the random number is less than the values of the previous probabilities combined
                        new_age = a-1; %if the random number falls in the corresponding range of an age class, assign the individual to that age class
                        break %the loop is broken when an age is assigned
                    end
                end
                indices(1,k) = options_lower(new_age+1) + round((options_upper(new_age+1)-options_lower(new_age+1))*rand);
                indices(2,k) = new_age;
                
            end
        end
    end
    age = stored_age;
end

