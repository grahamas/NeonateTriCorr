function fn2str(func::Function)
    strip(string(func), ['!'])
end