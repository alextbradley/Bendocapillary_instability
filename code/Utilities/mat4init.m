function yinit = mat4init(x, omega)
%initialising initial guess for BVP4c
yinit = [sin(omega*x)
    cos(omega*x)
    -sin(omega*x)
    -cos(omega*x)
    sin(omega*x)
    cos(omega*x)];
end