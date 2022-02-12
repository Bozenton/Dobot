function V = fV(dth1,dth2,dth3,g,l2,l3,m2,m3,th2,th3)

V = [dth1.*(dth2.*l2.^2.*m2.*sin(th2.*2.0)+dth2.*l2.^2.*m3.*sin(th2.*2.0).*2.0+dth2.*l3.^2.*m3.*sin(th2.*2.0+th3.*2.0).*4.0+dth3.*l3.^2.*m3.*sin(th2.*2.0+th3.*2.0).*4.0+dth2.*l2.*l3.*m3.*sin(th2.*2.0+th3).*4.0+dth3.*l2.*l3.*m3.*sin(th2.*2.0+th3).*2.0+dth3.*l2.*l3.*m3.*sin(th3).*2.0).*(-1.0./2.0);...
    (l2.*m2.*cos(th2))./2.0+l2.*m3.*cos(th2)+l3.*m3.*cos(th2+th3)-g.*l3.*m3.*cos(th2+th3)-(g.*l2.*m2.*cos(th2))./2.0-g.*l2.*m3.*cos(th2)+(dth1.^2.*l2.^2.*m2.*sin(th2.*2.0))./4.0+(dth1.^2.*l2.^2.*m3.*sin(th2.*2.0))./2.0+dth1.^2.*l3.^2.*m3.*sin(th2.*2.0+th3.*2.0)-dth3.^2.*l2.*l3.*m3.*sin(th3)+dth1.^2.*l2.*l3.*m3.*sin(th2.*2.0+th3)-dth2.*dth3.*l2.*l3.*m3.*sin(th3).*2.0;...
    (l3.*m3.*(cos(th2+th3).*2.0-g.*cos(th2+th3).*2.0+dth1.^2.*l3.*sin(th2.*2.0+th3.*2.0).*2.0+dth1.^2.*l2.*sin(th3)+dth2.^2.*l2.*sin(th3).*2.0+dth1.^2.*l2.*sin(th2.*2.0+th3)))./2.0];

end

