function T = Tmat(alp, a, d, th)
    sa = sin(alp);
    ca = cos(alp);
    st = sin(th);
    ct = cos(th);
    T = [ ct    -st    0    a
         st*ca  ct*ca  -sa  -sa*d
         st*sa  ct*sa  ca   ca*d
          0      0     0      1  ];
end

