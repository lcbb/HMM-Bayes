function list = struct2vector(structure,field,varargin)

list = reshape([structure.(field)],[size(structure(1).(field)) length(structure)]);

if nargin == 3
    if size(list,1)==1
        list = squeeze(list(1,varargin{1},:));
    else
        list = squeeze(list(varargin{1},1,:));
    end
    
elseif nargin == 4
    list = squeeze(list(varargin{1},varargin{2},:));

end


end