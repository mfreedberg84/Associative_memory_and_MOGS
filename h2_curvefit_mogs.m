
% Curvefit program
% ----------------------------------------------
data = Correct_sequence_data;
[totalrows, columns] = size(data);
outfile = 'Correct_sequence_data_output.txt';
% ----------------------------------------------
    
ycol = columns;
xcol = ycol - 1;
groupingcols = 1:(ycol-2);
groupdata = rand(1,ycol-2); 
%initialize grouping data to random number to kick off starting a dataset

currentrow = 0;
output = [];
   
warning off all

%disp('Sub     R'); Uncomment this and line 50 to see model fits

% Analyze each participants data individually
while currentrow < totalrows
    currentrow = currentrow+1;
    groupdata = data(currentrow,groupingcols);
    startrow = currentrow;
    stop = 0;
    while stop == 0
        if sum(groupdata ~= data(currentrow,groupingcols)) >= 1
            stop = 1;
        else
            currentrow = currentrow + 1;
            if currentrow > totalrows 
                stop = 1;
            end
        end
    end
    currentrow = currentrow - 1;
    
    x = data(startrow:currentrow,xcol);
    y = data(startrow:currentrow,ycol);
    %now we have a dataset, initiate the fit
	
	%get start data from exceptions file if we have it
    [bestparams, rr, ls, exitflag] = h3_exponentialfit(x,y);

    %fprintf('%4.0f %4.2f\n',groupdata(1,1),rr); uncomment to see model
    %fits
    output = [output; groupdata, bestparams, rr,ls];
    
end    %while currentrow <=totalrows

                        
