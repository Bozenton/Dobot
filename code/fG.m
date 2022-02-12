function G = fG(l2,l3,m2,m3,th2,th3)

G = [0.0;
    l2.*m2.*cos(th2).*(-1.0./2.0)-l2.*m3.*cos(th2)-l3.*m3.*cos(th2+th3);
    -l3.*m3.*cos(th2+th3)];

end

