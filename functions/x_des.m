function x = x_des(that,p1,p2)
% Design cl / clcd / alpha [ - / deg / - ] as a function of t/c [%]
x = NaN(length(that));
that(that < 24.1) = 24.1;
that(that > 100)  = 100;
for i = 1:length(that)
    if that(i) <= 48
        x(i) = p1(1)*that(i)^3 + p1(2)*that(i)^2 + p1(3)*that(i) + p1(4);
    else
        x(i) = p2(1)*that(i) + p2(2);
    end
end
end