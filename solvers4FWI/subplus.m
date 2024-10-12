function output = subplus(input)
%SUBPLUS Calculate positive part of function
%   Replacement for function in link below without Curve Fitting Toolbox
%   https://www.mathworks.com/help/curvefit/subplus.html
    output = input.*(input>0);
end

