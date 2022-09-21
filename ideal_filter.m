function y = ideal_filter(t, T1, f1)
  y = zeros(size(t));
  idx1 = -1*T1 <= t & t <= T1;
  y(idx1) = f1;
  y(~(idx1)) = 0;
  end