function varargout = lorentzfit(x,y,varargin)
% LORENTZFIT fits a single- or multi-parameter Lorentzian function to data
%
%   LORENTZFIT(X,Y) returns YPRIME(X), a Lorentzian fit to the data
%   found using LSQCURVEFIT. The function Y(X) is fit by the model:
%       YPRIME(X) = P1./((X - P2).^2 + P3) + C.
%
%   [YPRIME PARAMS RESNORM RESIDUAL] = LORENTZFIT(X,Y) returns YPRIME(X) 
%   values in addition to fit-parameters PARAMS = [P1 P2 P3 C]. The RESNORM
%   and RESIDUAL outputs from LSQCURVEFIT are also returned.
%
%   [...] = LORENTZFIT(X,Y,P0) can be used to provide starting
%   values (P0 = [P01 P02 P03 C0]) for the parameters in PARAMS.
%
%   [...] = LORENTZFIT(X,Y,P0,BOUNDS) may be used to define lower 
%   and upper bounds for the possbile values for each parameter in PARAMS.
%       BOUNDS = [LB1 LB2 LB3 LB4;
%                 UB1 UB2 UB3 UB4].
%   If the user does not wish to manually define values for P0, it may be
%   enetered as an empty matrix P0 = []. In this case, default values will
%   be used. The default bounds for all parameters are (-Inf,Inf).
%
%   [...] = LORENTZFIT(X,Y,P0,BOUNDS,NPARAMS) may be used to specify the
%   number of parameters used in the Lorentzian fitting function. The 
%   number of parameters defined in P0 and BOUNDS must match the function 
%   specified by NPARAMS. If the user does not wish to manually define 
%   values for P0 or BOUNDS, both may be enetered as empty matricies: 
%   P0 = []; BOUNDS = [].
%
%   -NPARAMS options
%       
%           '1'     - Single parameter Lorentzian (no constant term)
%                     L1(X) = 1./(P1(X.^2 + 1))
%
%           '1c'    - Single parameter Lorentzian (with constant term)
%                     L1C(X) = 1./(P1(X.^2 + 1)) + C
% 
%           '2'     - Two parameter Lorentzian (no constant term)
%                     L2(X) = P1./(X.^2 + P2)
%
%           '2c'    - Two parameter Lorentzian (with constant term)
%                     L2C(X) = P1./(X.^2 + P2) + C
%
%           '3'     - Three parameter Lorentzian (no constant term)
%                     L3(X) = P1./((X - P2).^2 + P3)
%
% [DEFAULT] '3c'    - Three parameter Lorentzian (with constant term)
%                     L3C(X) = P1./((X - P2).^2 + P3) + C
%
%   [...] = LORENTZFIT(X,Y,P0,BOUNDS,NPARAMS,OPTIONS) defines the OPTIONS
%   array for the MATLAB function LSQCURVEFIT. OPTIONS can be set using the
%   following command:
%
%       OPTIONS = optimset('PARAM1',VALUE1,'PARAM2',VALUE2,...);
%
%   See the help documentation for OPTIMSET for more details.
%   
% 
%   X and Y must be the same size, numeric, and non-complex. P0 and BOUNDS
%   must also be numeric and non-complex. NPARAMS is a character array.
%
%   Examples: 
%       x = -16:0.1:35;
%       y = 19.4./((x - 7).^2 + 15.8) + randn(size(x))./10;
%       [yprime1 params1 resnorm1 residual1] = lorentzfit(x,y,[20 10 15 0]);
%       figure; plot(x,y,'b.','LineWidth',2)
%       hold on; plot(x,yprime1,'r-','LineWidth',2)
%
%       [yprime2 params2 resnorm2 residual2] = lorentzfit(x,y,[],[],'3');
%       figure; plot(x,y,'b.','LineWidth',2)
%       hold on; plot(x,yprime2,'r-','LineWidth',2)
%
%   See also: lsqcurvefit.

% Jered R Wells
% 11/15/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.6 (2015/07/31)
%
% REF: http://www.home.uos.de/kbetzler/notes/fitp.pdf
%
%   UPDATES
%   v1.5 - 2015/07/22 - jrw
%       Added INPUTCHECK and OPTIONS. Rearranged SWITCH loop.
%   v1.6 - 2015/07/31 - jrw
%       Added stopping thresholds for fitting based on magnitude of input data.

%% CHECK INPUTS
narginchk(2,6);
nargoutchk(0,4);
fname = 'lorentzfit';

% Checked required inputs
inputcheck(x,{'numeric'},{'real','nonnan','nonempty','finite'},fname,'X',1);
inputcheck(y,{'numeric'},{'real','nonnan','nonempty','finite','size',size(x)},fname,'Y',2);

% Set defaults for optional inputs
p3 = ((max(x(:))-min(x(:)))./10).^2;
p2 = (max(x(:))+min(x(:)))./2;
p1 = max(y(:)).*p3; 
c = min(y(:));

optargs = {[],[],'3c',optimset('TolFun',max(mean(y(:))*1e-6,1e-15),'TolX',max(mean(x(:))*1e-6,1e-15))};

numvarargs = length(varargin);
for ii = 1:numvarargs; if isempty(varargin{ii}); varargin{ii} = optargs{ii}; end; end

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[p0,bounds,nparams,options] = optargs{:};

% Check optional inputs
if ~isempty(p0)
    inputcheck(p0,{'numeric'},{'real','nonnan','vector'},fname,'P0',3);
end
if ~isempty(bounds)
    inputcheck(bounds,{'numeric'},{'real','nonnan','nrows',2},fname,'BOUNDS',4);
    lb = bounds(1,:);
    ub = bounds(2,:);
else
    lb = [];
    ub = [];
end
inputcheck(nparams,{'char'},{},fname,'NPARAMS',5);
inputcheck(options,{'struct'},{},fname,'OPTIONS',6);

%% PROCESS
switch lower(nparams)
    case '1'
        % Define P0, LB, UB
        if isempty(p0);
            p1 = max(y(:)); p0 = p1;
        elseif numel(p0)~=1
            error 'P0 must be empty or have one element for NPARAMS = ''1'''; 
        end

        if isempty(bounds)
            lb = -Inf; ub = Inf;
        elseif ~all(size(bounds)==[2 1]) 
            error 'BOUNDS must be empty or it must be a 2x1 matrix for NPARAMS = ''1''';
        else
            lb = bounds(1,:); ub = bounds(2,:);
        end

        if any(lb>=ub)
            error 'Lower bounds must be less than upper bounds'; 
        end

        [params,resnorm,residual] = lsqcurvefit(@lfun1,p0,x,y,lb,ub,options);
        yprime = lfun1(params,x);

    case '1c'
        % Define P0, LB, UB
        if isempty(p0);
            p1 = max(y(:)); p0 = [p1 c];
        elseif numel(p0)~=2
            error 'P0 must be empty or have two elements for NPARAMS = ''1c'''; 
        end

        if isempty(bounds)
            lb = [-Inf,-Inf]; ub = [Inf,Inf];
        elseif ~all(size(bounds)==[2 2]) 
            error 'BOUNDS must be empty or it must be a 2x2 matrix for NPARAMS = ''1c''';
        else
            lb = bounds(1,:); ub = bounds(2,:);
        end

        if any(lb>=ub)
            error 'Lower bounds must be less than upper bounds'; 
        end

        [params,resnorm,residual] = lsqcurvefit(@lfun1c,p0,x,y,lb,ub,options);
        yprime = lfun1c(params,x);
    case '2'
        % Define P0, LB, UB
        if isempty(p0);
            p2 = ((max(x(:))-min(x(:)))./10).^2; 
            p1 = max(y(:)).*p2;
            p0 = [p1 p2];
        elseif numel(p0)~=2
            error 'P0 must be empty or have two elements for NPARAMS = ''2'''; 
        end

        if isempty(bounds)
            lb = [-Inf,-Inf]; ub = [Inf,Inf];
        elseif ~all(size(bounds)==[2 2]) 
            error 'BOUNDS must be empty or it must be a 2x2 matrix for NPARAMS = ''2''';
        else
            lb = bounds(1,:); ub = bounds(2,:);
        end

        if any(lb>=ub)
            error 'Lower bounds must be less than upper bounds'; 
        end

        [params,resnorm,residual] = lsqcurvefit(@lfun2,p0,x,y,lb,ub,options);
        yprime = lfun2(params,x);
    case '2c'
        % Define P0, LB, UB
        if isempty(p0);
            p2 = ((max(x(:))-min(x(:)))./10).^2; 
            p1 = max(y(:)).*p2;
            p0 = [p1 p2 c];
        elseif numel(p0)~=3
            error 'P0 must be empty or have three elements for NPARAMS = ''2c'''; 
        end

        if isempty(bounds)
            lb = [-Inf,-Inf,-Inf]; ub = [Inf,Inf,Inf];
        elseif ~all(size(bounds)==[2 3]) 
            error 'BOUNDS must be empty or it must be a 2x3 matrix for NPARAMS = ''2c''';
        else
            lb = bounds(1,:); ub = bounds(2,:);
        end

        if any(lb>=ub)
            error 'Lower bounds must be less than upper bounds'; 
        end

        [params,resnorm,residual] = lsqcurvefit(@lfun2c,p0,x,y,lb,ub,options);
        yprime = lfun2c(params,x);
    case '3'
        % Define P0, LB, UB
        if isempty(p0);
            p0 = [p1 p2 p3];
        elseif numel(p0)~=3
            error 'P0 must be empty or have three elements for NPARAMS = ''3'''; 
        end

        if isempty(bounds)
            lb = [-Inf,-Inf,-Inf]; ub = [Inf,Inf,Inf];
        elseif ~all(size(bounds)==[2 3]) 
            error 'BOUNDS must be empty or it must be a 2x3 matrix for NPARAMS = ''3''';
        else
            lb = bounds(1,:); ub = bounds(2,:);
        end

        if any(lb>=ub)
            error 'Lower bounds must be less than upper bounds'; 
        end

        [params,resnorm,residual] = lsqcurvefit(@lfun3,p0,x,y,lb,ub,options);
        yprime = lfun3(params,x);
    case '3c'
        % Define P0, LB, UB
        if isempty(p0);
            p0 = [p1 p2 p3 c];
        elseif numel(p0)~=4
            error 'P0 must be empty or have four elements for NPARAMS = ''3c'''; 
        end

        if isempty(bounds)
        elseif ~all(size(bounds)==[2 4]) 
            error 'BOUNDS must be empty or it must be a 2x4 matrix for NPARAMS = ''3c''';
        else
            lb = bounds(1,:); ub = bounds(2,:);
        end

        if any(lb>=ub)
            error 'Lower bounds must be less than upper bounds'; 
        end

        [params,resnorm,residual] = lsqcurvefit(@lfun3c,p0,x,y,lb,ub,options);
        yprime = lfun3c(params,x);
otherwise
    warning('MATLAB:lorentzfit:default','Fitting default (''3c'') Lorentzian model.')
    [params,resnorm,residual] = lsqcurvefit(@lfun3c,p0,x,y,lb,ub,options);
    yprime = lfun3c(params,x);
end

varargout = {yprime,params,resnorm,residual};

end % MAIN

function F = lfun1(p,x)
F = 1./(p.*(x.^2 + 1));
end % LFUN1

function F = lfun1c(p,x)
F = 1./(p(1).*(x.^2 + 1)) + p(2);
end % LFUN1C

function F = lfun2(p,x)
F = p(1)./(x.^2 + p(2));
end % LFUN2

function F = lfun2c(p,x)
F = p(1)./(x.^2 + p(2)) + p(3);
end % LFUN2C

function F = lfun3(p,x)
F = p(1)./((x-p(2)).^2+p(3));
end % LFUN3

function F = lfun3c(p,x)
F = p(1)./((x-p(2)).^2+p(3)) + p(4);
end % LFUN3C

function varargout = inputcheck(A,varargin)
% INPUTCHECK checks the validity of input array with VALIDATEATTRIBUTES 
%
% INPUTCHECK(A,CLASSES,ATTRIBUTES) validates that array A belongs
% to at least one of the specified CLASSES and has all of the specified
% ATTRIBUTES. If A does not meet the criteria, MATLAB issues a formatted
% error message. INPUTCHECK(A,CLASSES,ATTRIBUTES) is equivalent to
% VALIDATEATTRIBUTES(A,CLASSES,ATTRIBUTES) with the exception that custom 
% CLASSES and ATTRIBUTES can be easily implemented by editing the function.
%
% Some additional input error-checking is also provided. Also, empty arrays
% may be defined if any of the input arguments is irrelevant.
%
% INPUTCHECK(A,CLASSES,ATTRIBUTES,ARGINDEX) includes the
% position of the input in your function argument list as part of any
% generated error messages.
% 
% INPUTCHECK(A,CLASSES,ATTRIBUTES,FUNCNAME) includes the
% specified function name in generated error identifiers.
% 
% INPUTCHECK(A,CLASSES,ATTRIBUTES,FUNCNAME,VARNAME) includes the
% specified variable name in generated error messages.
% 
% INPUTCHECK(A,CLASSES,ATTRIBUTES,FUNCNAME,VARNAME,ARGINDEX)
% includes the specified information in the generated error messages or
% identifiers.
%
% [V,ME] = INPUTCHECK(A,CLASSES,ATTRIBUTES,FUNCNAME,VARNAME,ARGINDEX,VERBOSE)
% toggles the error message generation state. VERBOSE = TRUE will produce 
% error messages which halt function execution and return control to the 
% command line if array A is not valid. Otherwise, the validity state V 
% will be output as logical TRUE and the error catch variable ME is empty.
% VERBOSE = FALSE will output validity V as either TRUE or FALSE in 
% addition to reporting errors in ME without halting function execution 
% when A is not valid. By default, VERBOSE = TRUE.
% 
%   Input Arguments:
%
%   A          Any class of array.
%
%   CLASSES    Cell array of strings that specify valid classes for array A.
%              For example, if CLASSES = {'logical','cell'}, A must be a
%              logical array or a cell array. The string 'numeric' is an
%              abbreviation for the classes uint8, uint16, uint32, uint64,
%              int8, int16, int32, int64, single, double. CLASSES may 
%              include MATLAB built-in or custom classes:
%
%             'numeric'         Any value for which the isnumeric function 
%                               returns true, including int8, int16, int32, 
%                               int64, uint8, uint16, uint32, uint64, 
%                               single, or double.
%             'single'          Single-precision number.
%             'double'          Double-precision number.
%             'int8'            Signed 8-bit integer.
%             'int16'           Signed 16-bit integer.
%             'int32'           Signed 32-bit integer.
%             'int64'           Signed 64-bit integer.
%             'uint8'           Unsigned 8-bit integer.
%             'uint16'          Unsigned 16-bit integer.
%             'uint32'          Unsigned 32-bit integer.
%             'uint64'          Unsigned 64-bit integer.
%             'logical'         Logical true or false.
%             'char'            Character or string.
%             'struct'          Structure array.
%             'cell'            Cell array.
%             'function_handle'	Scalar function handle.
%
%             Examples: 
%                   % Define CLASSES using the following syntax
%                   classes = {};
%                   classes = {'double'};
%                   classes = {'int8','int16'};
%
%   ATTRIBUTES Cell array that contains descriptions of valid attributes
%              for array A. For example, if ATTRIBUTES = {'real','finite'}
%              A must contain only real and finite values.
%
%              Supported attributes include:
%   
%             Attributes that describe the size and shape of array A:
%             '2d'              Two-dimensional array, including scalars, 
%                               vectors, matrices, and empty arrays.
%             'column'          Column vector, N-by-1.
%             'row'             Row vector, 1-by-N.
%             'scalar'          Scalar value, 1-by-1.
%             'vector'          Row or column vector, or a scalar value.
%             'size', [d1,...,dN]	Array with dimensions d1-by-...-by-dN. 
%                                   If you do not want to check a 
%                                   particular dimension, specify NaN for 
%                                   that dimension, such as [3,4,NaN,2].
%               'ndims', N      Array with N dimensions.
%             'numel', N        Array with N elements.
%               'numel>', N     Array with more than N elements.
%               'numel<', N     Array with fewer than N elements.
%             'ncols', N        Array with N columns.
%               'ncols>', N     Array with more than N columns.
%               'ncols<', N     Array with fewer than N columns.
%             'nrows', N        Array with N rows.
%               'nrows>', N     Array with more than N rows.
%               'nrows<', N     Array with fewer than N rows.
%             'nonempty'        No dimensions equal zero.
%             'nonsparse'       Array that is not sparse.
%
%             Attributes that specify valid ranges for values in A:
%             '>', N            All values are greater than N.
%             '>=', N           All values are greater than or equal to N.
%             '<', N            All values are less than N.
%             '<=', N           All values are less than or equal to N.
%
%             Attributes that check types of values in A, where A is a 
%             numeric or logical array:
%             'binary'          Array of ones and zeros.
%             'even'            All elements are even integers 
%                               (includes zero).
%             'odd'             All elements are odd integers.
%             'integer'         All elements are integer-valued.
%             'real'            All elements are real.
%             'finite'          All elements are finite.
%             'nonnan'          No elements equal to NaN (Not a Number).
%             'nonnegative'     All elements are nonnegative.
%             'nonzero'         All elements are nonzero.
%             'positive'        All elements are positive.
%
%              Some attributes also require numeric values. For those 
%              attributes, the numeric value or vector must immediately 
%              follow the attribute name string. For example,
%              {'>=', 5, '<=', 10, size, [3 4 2]} checks that all
%              values of A are between 5 and 10, and that A is 3-by-4-by-2.
%
%             Examples: 
%                   % Define ATTRIBUTES using the following syntax
%                   attributes = {};
%                   attributes = {'2d'};
%                   attributes = {'2d','>',5,'numel',1000};
%
%   ARGINDEX   Positive integer that specifies the position of the input
%              argument.
%
%   FUNCNAME   String that specifies the function name. If you specify an
%              empty string, '', FUNCNAME is ignored.
%
%   VARNAME    String that specifies input argument name. If you specify an
%              empty string, '', VARNAME is ignored.
%
%   VERBOSE    Logical that toggles the error-generation state of
%              INPUTCHECK.
%
%   Output Arguments:
%
%   V          (Logical) validity state of A.
%
%   ME         Class MException object for error capture.
%
%   Example: Create a three dimensional array and then check for the
%            attribute '2d'.
%
%       A = [ 1 2 3; 4 5 6 ];
%       B = [ 7 8 9; 10 11 12];
%       C = cat(3,A,B);
%       inputcheck(C,{'numeric'},{'2d'},'my_func','my_var',2)
%
%   This code throws an error and displays a formatted message:
%
%       Expected input number 2, my_var, to be two-dimensional.
%
% [REF: http://www.mathworks.com/help/techdoc/ref/validateattributes.html]
%
% See also: validateattributes.

%
% Jered R Wells
% 08/15/2012
%
% v1.0 (08/17/2012)

%% Initialize variables
V = true;
ME = MException('','');
%#ok<*SPERR>

classes = {};
attributes = {};
funcname = '';
varname = '';
argindex = [];
verbose = true;

%% Input check
if nargin==1
    % This case is not very useful as no error checking actually occurs!
    warning('MATLAB:inputcheck:noCrit','No criteria for array validation')
    classes = {};
    attributes = {};
    validateattributes(A,classes,attributes);
elseif nargin>1
    for ii = 1:length(varargin)
        switch ii
            case 1
                classes = varargin{ii};
                if isempty(classes); classes = {}; end
                if ischar(classes); classes = {classes}; end
                validateattributes(classes,{'cell'},{},'INPUTCHECK','CLASSES',2);
            case 2
                attributes = varargin{ii};
                if isempty(attributes); attributes = {}; end
                validateattributes(attributes,{'cell'},{},'INPUTCHECK','ATTRIBUTES',3);
            case 3
                funcname = varargin{ii};
                if isempty(funcname); funcname = ''; end
                if isnumeric(funcname)&&length(varargin)==3; 
                    funcname = '';
                    argindex = varargin{ii};
                    if isempty(argindex)
                        argindex = []; 
                    else
                        validateattributes(argindex,{'numeric'},{'positive','integer','numel',1},'INPUTCHECK','ARGINDEX',3);
                    end
                end
                validateattributes(funcname,{'char'},{},'INPUTCHECK','FUNCNAME',4);
            case 4
                varname = varargin{ii};
                if isempty(varname); varname = ''; end
                validateattributes(varname,{'char'},{},'INPUTCHECK','VARNAME',5);
            case 5
                argindex = varargin{ii};
                if isempty(argindex)
                    argindex = []; 
                else
                    validateattributes(argindex,{'numeric'},{'positive','integer','numel',1},'INPUTCHECK','ARGINDEX',6);
                end
            case 6
                verbose = varargin{ii};
                validateattributes(verbose,{'logical'},{'numel',1},'INPUTCHECK','VERBOSE',7);
            otherwise
                error 'Too many input arguments'
        end % SWITCH
    end % ii
else
    error 'Array A must be defined'
end % IF

% Format strings
if isempty(funcname); funcnameIC = '';
else funcnameIC = ['Error using ',funcname];
end
if isempty(argindex); argindexIC = '';
else argindexIC = ['number ',num2str(argindex),',',' '];
end
if isempty(varname); varnameIC = '';
else varnameIC = [varname,',',' '];
end
    

%% Process
% Check for custom CLASSES
% N/A

% Check for custom ATTRIBUTES
ii = 1;
customAttributes = true(length(attributes));
while ii<=length(attributes)
    if ischar(attributes{ii})
        tmpName = lower(attributes{ii});
        switch tmpName
            % Create a new case for each custom ATTRIBUTE
            % Be sure to check for the requirement of N
            case 'ndims'
                % Check for N
                ii = ii + 1;                % Advance ii if N is required
                tmpN = attributes{ii};
                validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','ndims',3);
                % Validate array A
                if ndims(A)~=tmpN
                    % Formatted string
                    MSGID = 'MATLAB:inputcheck:incorrectNdims';
                    ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of dimensions equal to %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                    if verbose
                        error(ERRMSG)
                    else
                        ME = MException(MSGID, ERRMSG);
                        V = false;
                        break
                    end
                end
            case 'numel>'
                % Check for N
                ii = ii + 1;                % Advance ii if N is required
                tmpN = attributes{ii};
                validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','numel>',3);
                % Validate array A
                if ~(numel(A)>tmpN)
                    % Formatted string
                    MSGID = 'MATLAB:inputcheck:incorrectNumel';
                    ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of elements greater than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                    if verbose
                        error(ERRMSG)
                    else
                        ME = MException(MSGID, ERRMSG);
                        V = false;
                        break
                    end
                end
            case 'numel<'
                % Check for N
                ii = ii + 1;                % Advance ii if N is required
                tmpN = attributes{ii};
                validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','numel<',3);
                % Validate array A
                if ~(numel(A)<tmpN)
                    % Formatted string
                    MSGID = 'MATLAB:inputcheck:incorrectNumel';
                    ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of elements less than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                    if verbose
                        error(ERRMSG)
                    else
                        ME = MException(MSGID, ERRMSG);
                        V = false;
                        break
                    end
                end
            case 'ncols>'
                % Check for N
                ii = ii + 1;                % Advance ii if N is required
                tmpN = attributes{ii};
                validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','ncols>',3);
                % Validate array A
                if ~(size(A,2)>tmpN)
                     % Formatted string
                    MSGID = 'MATLAB:inputcheck:incorrectNcols';
                    ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of columns greater than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                    if verbose
                        error(ERRMSG)
                    else
                        ME = MException(MSGID, ERRMSG);
                        V = false;
                        break
                    end
                end
            case 'ncols<'
                % Check for N
                ii = ii + 1;                % Advance ii if N is required
                tmpN = attributes{ii};
                validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','ncols<',3);
                % Validate array A
                if ~(size(A,2)<tmpN)
                    % Formatted string
                    MSGID = 'MATLAB:inputcheck:incorrectNcols';
                    ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of columns less than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                    if verbose
                        error(ERRMSG)
                    else
                        ME = MException(MSGID, ERRMSG);
                        V = false;
                        break
                    end
                end
            case 'nrows>'
                % Check for N
                ii = ii + 1;                % Advance ii if N is required
                tmpN = attributes{ii};
                validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','nrows>',3);
                % Validate array A
                if ~(size(A,1)>tmpN)
                    % Formatted string
                    MSGID = 'MATLAB:inputcheck:incorrectNrows';
                    ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of rows greater than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                    if verbose
                        error(ERRMSG)
                    else
                        ME = MException(MSGID, ERRMSG);
                        V = false;
                        break
                    end
                end
            case 'nrows<'
                % Check for N
                ii = ii + 1;                % Advance ii if N is required
                tmpN = attributes{ii};
                validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','nrows<',3);
                % Validate array A
                if ~(size(A,1)<tmpN)
                    % Formatted string
                    MSGID = 'MATLAB:inputcheck:incorrectNrows';
                    ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of rows less than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                    if verbose
                        error(ERRMSG)
                    else
                        ME = MException(MSGID, ERRMSG);
                        V = false;
                        break
                    end
                end
            otherwise
                % If TMPNAME is not a custom attribute, assign FALSE to the
                % appropriate element in CUSTOMATTRIBUTES
                customAttributes(ii) = false;
        end % SWITCH
    else
        % If ATTRIBUTES{ii} is not a custom attribute, assign FALSE to the
        % appropriate element in CUSTOMATTRIBUTES
        customAttributes(ii) = false;
    end % IF
    ii = ii + 1;
end % ii
attributes = attributes(~customAttributes);

if verbose
    validateattributes(A,classes,attributes,funcname,varname,argindex)
else
    try
        validateattributes(A,classes,attributes,funcname,varname,argindex)
    catch ME
        V = false;
    end
end % IF

varargout = {V,ME};

end % MAIN
