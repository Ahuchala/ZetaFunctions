-- sage passes in arguments to command line to run in M2
-- this processes them

parseInput = () -> (
    function_name = value(scriptCommandLine#1);
    n = value(scriptCommandLine#2);
    polynomials = toString scriptCommandLine#3;
    if length scriptCommandLine > 4 then return function_name(n,polynomials, value(scriptCommandLine#4));

    -- print scriptCommandLine#1;
    return function_name(n,polynomials);
)
