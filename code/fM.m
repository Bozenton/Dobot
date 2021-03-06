function M = fM(l2,l3,m1,m2,m3,r,th2,th3)
M = reshape( ...
    [(l2.^2.*m2)./4.0+(l2.^2.*m3)./2.0+l3.^2.*m3+(m1.*r.^2)./2.0+l3.^2.*m3.*cos(th2.*2.0+th3.*2.0)+(l2.^2.*m2.*cos(th2.*2.0))./4.0+(l2.^2.*m3.*cos(th2.*2.0))./2.0+l2.*l3.*m3.*cos(th3)+l2.*l3.*m3.*cos(th2.*2.0+th3),0.0,0.0, ...
    0.0,    (l2.^2.*m2)./2.0+l2.^2.*m3+l3.^2.*m3.*2.0+l2.*l3.*m3.*cos(th3).*2.0,l3.*m3.*(l3.*2.0+l2.*cos(th3)), ...
    0.0,    l3.*m3.*(l3.*2.0+l2.*cos(th3)),     l3.^2.*m3.*2.0],...
    [3,3]);

end

