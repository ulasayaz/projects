fprintf('(((((((((((((cond: %d (%g), test: %d, countSuccess: %d, countFailure: %d)))))))))))))\n',cond,conditions(cond),test,countSuccess,countFailure)

pass = false;
while pass == false
    try
        results(test,cond) = example_TV2_2D_theory_function(settings);
        pass = true;
        countSuccess = countSuccess+1;
    catch ME
        ME.identifier
        countFailure = countFailure + 1;
    end
    close all
end