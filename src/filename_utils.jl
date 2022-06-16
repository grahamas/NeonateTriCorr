function fn2str(func::Function)
    strip(string(func), ['!'])
end
function obj2str(obj)
    strip(string(obj), ['(', ')'])
end