

function varargout = derivative5(im, varargin)

    varargin = varargin(:);
    varargout = cell(size(varargin));
    
   
    secondDeriv = false;    
    for n = 1:length(varargin)
        if length(varargin{n}) > 1
            secondDeriv = true;
            break
        end
    end
    
    if ~secondDeriv
      
        p = [0.037659  0.249153  0.426375  0.249153  0.037659];
        d1 =[0.109604  0.276691  0.000000 -0.276691 -0.109604];
    else         
       
        p  = [0.030320  0.249724  0.439911  0.249724  0.030320];
        d1 = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
        d2 = [0.232905  0.002668 -0.471147  0.002668  0.232905];
    end

    gx = false;
    
    for n = 1:length(varargin)
      if strcmpi('x', varargin{n})
          varargout{n} = conv2(p, d1, im, 'same');    
          gx = true;   
          gxn = n;
      elseif strcmpi('y', varargin{n})
          varargout{n} = conv2(d1, p, im, 'same');
      elseif strcmpi('xx', varargin{n})
          varargout{n} = conv2(p, d2, im, 'same');    
      elseif strcmpi('yy', varargin{n})
          varargout{n} = conv2(d2, p, im, 'same');
      elseif strcmpi('xy', varargin{n}) | strcmpi('yx', varargin{n})
          if gx
              varargout{n} = conv2(d1, p, varargout{gxn}, 'same');
          else
              gx = conv2(p, d1, im, 'same');    
              varargout{n} = conv2(d1, p, gx, 'same');
          end
      else
          error(sprintf('''%s'' is an unrecognized derivative option',varargin{n}));
      end
    end
    
