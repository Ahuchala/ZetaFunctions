-- sage passes in arguments to command line to run in M2
-- this processes them

parseInput = () -> (
    function_name = value(scriptCommandLine#1);
    n = value(scriptCommandLine#2);
    polynomials = toString scriptCommandLine#3;
    return function_name(n,polynomials);
)
