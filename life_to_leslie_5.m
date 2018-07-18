function [life_table_m,leslie_matrix,ages,orig_lambda,mod_lambda,CV_fecundity,G,alpha, AL,mean_fecundity] = life_to_leslie_5(file_path_name,scale)

    % A function to input a life table and produce a Leslie matrix that can be used to simulate population demography over time. 

    % This function requires a life table in a .csv or an .xlsx file in the form of an number of ages x 3 matrix. 
    % Column 1 contains the ages as integers, Column 2 contains age-specific survival as a real number, Column 3 contains age-specific fecundity as a real number. 
    % The oldest age class in the life table should always have zero probability of
    % survival, although the age-specific fecundity need not be zero. 
    % Life tables should NOT include a final age class where both survival and fecundity are zero.
    %
    % An example life table file would contain:

    %     Example life table in .xlsx file (row with title string...)		
    %     Age	Sx	Bx (row with column headings)
    %     1	0.9	0.000
    %     2	0.6	0.429
    %     3	0.4	0.310
    %     4	0.0	0.120

    % Inputs:
    %
    % file_path_name - path name to life table file
    % scale - scalar rela number that is used to rescale the population
    % growth rate. To achieve a population growth rate of one, all elements of the Leslie matrix are divided by the leading eigenvalue.
    % If the leading eigenvalue is multiplied by a scale value, the growth
    % rate can be altered to be more or less than one and the rate of
    % growth or shrinkage altered.
    % when scale = 1.0 - population growth rate of exactly 1.0
    % when scale > 1.0 - negative population growth rate
    % when scale < 1.0 - positive population growth rate
    
    % Outputs:
    % life_table_m - the numeric values of the life table, an ages x three matrix of real numbers
    % leslie_matrix - corresponding Leslie Matrix that is square ages x ages,
    %	rescaled to constant population size over time if rescale = true.
    % ages - integer number of age classes in life table
    % orig_lambda - real number pop growth rate of original life table
    % mod_lambda - real number pop growth rate after rescaling
    % CV_fecundity - scalar real number coefficient of variance in age-specific
    %   fecundity
    % G - scalar real number generation length, G, or mean age of a parent
    % alpha - scalar real number age when bx first > zero
    % AL - integer adult lifespan

    
    % use importdata to read life table file since it does not require
    % knowledge of file type, works for Excel and for .csv
    file_contents = importdata(file_path_name);
    
    % data cell contains numeric life table
    life_table_m = file_contents.data;

    % determine the number of rows of the life_table to create a correctly sized leslie matrix
    [ages] = size(life_table_m,1); 

    % create a square leslie matrix with dimension = number of ages in the life table
    leslie_matrix = zeros(ages,ages); 

    for i = 1:ages-1 
       
        leslie_matrix(1,i) = life_table_m(i,3); %set the first row of the leslie matrix to the age-specific fecundity coefficients
       
        leslie_matrix(i+1,i) = life_table_m(i,2); %set the diagonal (beginning at (2,1)) equal to the age-specific survival coefficients
    end

    leslie_matrix(1,end) = life_table_m(end,3);
    
    % get eigenvalues of the original Leslie matrix, leading postive eigenvalue is equal to pop growth rate
    lambda = eig(leslie_matrix); 

    if lambda(1) > 0
        orig_lambda = lambda(1);
    else
        orig_lambda = lambda(2);
    end

    
    % rescale Leslie matrix to constant population size over time
    leslie_matrix = leslie_matrix./(orig_lambda*scale);

    
    % get eigenvalues of Leslie matrix after rescaling 
    lambda = eig(leslie_matrix); 
    
    if lambda(1) > 0
        mod_lambda = lambda(1);
    else
        mod_lambda = lambda(2);
    end

    
    % compute the generation length, G, or mean age of a parent
    products = life_table_m(:,2).*life_table_m(:,3); % get product of (Sx)x(Bx) for each age class
    G = sum(products); % sum over all age classes

    % find age when bx first > zero
    counter = 1; % counter
    while life_table_m(counter,3) == 0
        counter = counter + 1;
    end
    alpha = counter; % age of reproductive maturity
    
    % compute adult lifespan
    AL = ages - alpha + 1;
    
    % compute coefficient of variation of age-specific fecundity
    var_fecundity = var(leslie_matrix(1,:));
    mean_fecundity = mean(leslie_matrix(1,:));
    CV_fecundity = sqrt(var_fecundity)/mean_fecundity;
    
    
    

