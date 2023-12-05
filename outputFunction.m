function stop = outputFunction(x,optimValues,state)

% outputFunction for optimizers (fminunc, fmincon etc), saving intermediate results 
% in global variable "outputFcn_global_data" for later access. 
%
% It is not supposed for live updates during optimization but for 
% later inspection, which is much more efficient. 
%
% Usage 
%   options = optimoptions( ... , 'OutputFcn',@outputFcn_global ); 
%   [XOpt,fval,exitflag,output] = fminunc(@fun, X0, options); 
%   outputFcn_global_data(k).x 
%
% Link to original code: https://www.mathworks.com/matlabcentral/fileexchange/80831-outputfcn_global

stop = false;
global outputGlobal
switch state
  case 'init'
    outputGlobal = struct();
    outputGlobal.x = x;
    outputGlobal.eval = optimValues.fval;
  case 'iter'
    ind = length(outputGlobal)+1;
    outputGlobal(ind).x = x;
    outputGlobal(ind).eval = optimValues.fval;
  case 'done'

  otherwise
    error('yikes. you fucked up, idiot.')
end
end
