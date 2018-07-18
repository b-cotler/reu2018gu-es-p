function [mrca, complete_genealogy, coal_events, age_zero_counter] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m)

%Inputs:
%1. genealogy_m - a number_generations x lineage_count x 2 matrix. At the
%time of input, the matrix should be filled with a null value of -1 at
%every index except it should be filled with the indices and ages of
%individuals in it's final surface/row.
%2. leslie_matrix -  a leslie matrix corresponding to the input life table,
%generated using the life_to_leslie function
%3. age_dist_m - the age distribution/demographic matrix containing the
%population of each age class over the number of generations

%this function traces a preset number of lineages (lineage_count) back in
%time to find a most recent common ancestor. The function is written
%iteratively and it's loops are broken if a common ancestor is reached.

%Outputs:
%1. mrca - the time to most recent common ancestor. Mrca is equal to the maximum
%number of generations if no most recent common ancestor is reached.

%2. complete_genealogy - a completed genealogy_m. Page one contains the indices of
%individuals in each lineage for every generation until a most recent
%common ancestor is reached. Page two contains the ages of individuals in
%each lineage for every generation until a most recent common ancestor is
%reached. When an MRCA is reached, the lineages will converge and all
%remaining cells will be filled with "-1".


%3. coal_events - coal events is a vector (when lineage_count = 2) or a matrix (when lineage_count > 2)
% where the first column is records the lineage lost in a coalescent event,
% the second column is the remaining lineage after a coalescent event, and
% the third column is the generation count when the coalescent event
% occurs.

%Establish the dimensions of the genealogy matrix to set bounds for
%iteration
generations = size(genealogy_m,1); 
lineages = 1:size(genealogy_m,2); 

%Create a counter to record the number of coalescent events
coal_count = 0;

%Establish a matrix for recording coalescent events
coal_events = -1*ones(length(lineages)-1,3);

%Create a counter to measure the time to most recent common ancestor (MRCA)
time_count = 1; 

%create a counter for age zero individuals to better understand
%polymorphism/mutation
age_zero_counter = zeros(1,length(lineages)); 

for g = generations:-1:2 %iterate over the generations
    
    %Establish vectors to record the bounds of each age class for each
    %generation
    options_lower = [1]; %initializes lower bound vector for each age class
    options_upper = [age_dist_m(1,g-1)]; %initializes upper bound vector for each age class
    for a = 2:size(age_dist_m)
        options_lower = [options_lower sum(age_dist_m(1:a-1,g-1))+1]; %concatenates vector of lower bounds
        options_upper = [options_upper sum(age_dist_m(1:a,g-1))]; %concatenates vector of upper bounds
    end
    %iterate over each lineage
    for q = 1:length(lineages) %iterate over the number of lineages
        k = lineages(q);
        
        %consider the case where the individual at generation g and lineage
        % k is age zero and must pick a parent. 
        if genealogy_m(g,k,2) == 0 
            age_zero_counter(k) = age_zero_counter(k)+1;
            fecundities = age_dist_m(:,g-1).*transpose(leslie_matrix(1,:)); %create a vector of expected number of offspring
            scaling = sum(fecundities); % determines the value to scale the fecundities
            fecundities = fecundities./scaling; %scales the fecundities vector so the sum is one
            parent_age_index = rand; %chooses a uniformly dist random number between 0 and 1
            for a = 1:length(fecundities)
                if parent_age_index <= sum(fecundities(1:a)) %check if the random number is less than the values of the previous scaled fecundities combined
                    parent_age = a-1; %if the random number is in the desired range, assign the parent the age corresponding to that range
                    break %the loop is broken when a parent age is assigned
                end
            end
            
         %choose an individual in the given population through random
         %number generation and assign it the corresponding age
            genealogy_m(g-1,k,1) = options_lower(parent_age+1) + round((options_upper(parent_age+1)-options_lower(parent_age+1))*rand);
            genealogy_m(g-1,k,2) = parent_age; 
            
        %case where the individual is not a newborn and age-1 ancestor must be chosen
        else   
            %establish the age of the individual by thinking about the age
            %in the generation before. 
            age_old = genealogy_m(g,k,2); %set an age variable equal to the current age in lineage k            
            genealogy_m(g-1,k,1) = options_lower(age_old) + round((options_upper(age_old)-options_lower(age_old))*rand); %choose a random individual in the age_old-1 generation
            genealogy_m(g-1,k,2) = age_old-1; %set the age of the individual in the previous generation
            
            %check to make sure the same individual hasn't been assigned
            %twice. If so, reassign the individual in the lineage with the
            %higher k value
            if k>1
                for r = 1:k-1
                    if isequal(genealogy_m(g-1,r,2),0)
                        break
                    else
                        loop_check = 0;
                        while isequal(genealogy_m(g-1,r,1),genealogy_m(g-1,k,1)) && isequal(genealogy_m(g,r,2),genealogy_m(g,k,2))
                            genealogy_m(g-1,k,1) = options_lower(age_old) + round((options_upper(age_old)-options_lower(age_old))*rand);
                           
              %Create a check to break the loop if there is an infinite
              %loop
                            loop_check = loop_check+1;
                            if isequal(loop_check,100)
                                disp("Warning: Infinite Loop");
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    
    %check for coalescent events, record and remove a lineage if a
    %coalescent event has occurred
    for r = 1:length(lineages-1)
        for s = r+1:length(lineages)
            if (isequal(genealogy_m(g-1,r,1), genealogy_m(g-1,s,1))) %if there is a coalescent event
                lineages(s)=[]; %remove the later lineage
                coal_count = coal_count+1; %add to the number of coal events
                coal_events(coal_count,:) = [s r time_count]; %add to coal events matrix
            end
        end
    end
    
    %set conditions to break the loop - lineages == 1 or time_count ==
    %generations
    time_count = time_count + 1; %add one to time
    if isequal(length(lineages),1)
        break %breaks the loop if the number of lineages reaches 1(at an mrca)
    end
    if isequal(time_count,generations)
        break %breaks the loop if the maximum number of generations have been reached.
    end
end
complete_genealogy = genealogy_m;
mrca = time_count;