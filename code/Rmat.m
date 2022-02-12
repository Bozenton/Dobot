function R = Rmat(th, alp)
    sa = sin(alp); ca = cos(alp);
    st = sin(th); ct = cos(th);
    R = [ ct    -st    0
        st*ca  ct*ca  -sa
        st*sa  ct*sa  ca];
end

