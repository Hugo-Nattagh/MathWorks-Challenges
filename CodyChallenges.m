% Mathworks Cody Challenges
% Practicing MatLab Techniques
a = 7;
b = 4;
c = 1:10;
d = -15:25;
e = [1 2 3 4 5 4 6 7 8 9];
f = [2 5 7 9];
g =[-4 6 3+4i 1+i 0];
h = [1 2 5 4 9 8 6 5 2 1 1 5 8 8 4 4 4 4 7 5 2];
i = [0 28 4 5 0 9 1 0 8 1 14 57 0];
j = [1 4 8 17 3];
k = [0 1 0 0 0 1 1 1 0 0 1 1 1 1 0 0 0 0 1];
i = [NaN 28 4 NaN 5 0 9 NaN 1 0 8 NaN 1 14 57 0];

A = [1 2 3; 4 5 6];
B = [1 2 3; 4 5 6; 7 8 9];
C = [1 2 3 4 5 6; 7 8 9 10 11 12; 13 14 15 16 17 18];
D = [1 2 3 4 5 6; 7 8 9 NaN 11 12; 13 NaN 15 NaN 17 18];
E = [1 1 0 0; 0 0 0 5; 2 0 0 0; 0 1 1 3];
F = [0 0; 1 0; 2 2; 0 1];

z = mostDistant(F)

% 1 - Times 2

function y = times2(x)
    y = 2*x;
end

% 2 - Make the Vector [1 2 3 4 5 6 7 8 9 10]

function x = OneToTen
    x = 1:10;
end

% 3 - Find the Sum of all the Numbers of the Input Vector

function y = vecsum(x)
    y = sum(x);
end

% 4 - Add Two Numbers

function c = add_two_numbers(a,b)
    c = a+b;
end

% 5 - Determine if an ipunt if odd

function tf = is_it_odd(n)
    if mod(n,2) == 0
        tf = 'Not Odd';
    else
        tf = 'Odd';
    end
end

% 6 - Select every other Element of a Vector

function y = everyOther(x)
    y = x(1, 1:2:end);
end

% 7 - Pizza!

function y = pizza(z, a)
    y = pi*(z^2)*a;
end

% 8 - Is my Wife right?

function out = wiferight(in)
    out = 'yes';
end

% 9 - Column Removal

function B = column_removal(A, n)
    if n == 1
        B = A(:, n+1:end);
    elseif n == length(A)
        B = A(:, 1:n-1);
    else
        B = [A(:, 1:n-1), A(:, n+1:end)];
    end
end

% 10 - Triangle Numbers

function t = triangle(n)
    t = n;
    for i = 1:n-1
        n = n-1;
        t = t + n;
    end
end

% 11 - Find all Elements less than 0 or greater than 10 and replace them
% with NaN

function y = cleanUp(x)
    y = x
    for i = 1:length(x)
        if x(:, i) < 0
            y(:, i) = NaN;
        elseif x(:, i) > 10
            y(:, i) = NaN;
        end
    end            
end

% 12 - Make a Checkerboard Matrix

function a = checkerboard(n)
    a = ones(n);
    r = 0;
    t = 0;
    for i = 1:n
        t = ~t;
        a(i,1) = t;
        r = ~t;
        for j = 2:n
            a(i, j) = r;
            if r == 0
                r = 1;
            else
                r = 0;
            end
        end
    end      
end

% 13 - Swap the First and the Last Column

function B = swap_ends(A)
    F = A(:,1);
    L = A(:,end);
    M = A(:,2:end-1);
    B = [ L, M, F ];
end

% 14 - Determine wether a Vector is Monotonically Increasing

function tf = mono_increase(x)
    for i = 1:length(x)-1
        if x(i) < x(i+1)
            tf = true;
        else
            tf = false;
            break
        end        
    end
end

% 15 - Fibonacci Sequence

function f = fib(n)
    Phi = (1+sqrt(5))/2;
    phi = -1/Phi
    f = (1/sqrt(5))*((Phi^n)-(phi^n))
end

% 16 - Create Times-Tables

function m = timestables(n)
    m = ones(n);
    for i = 1:n
        p = 0;
        for j = 1:n
            p = p + i;
            m(i,j) = p;
        end
    end
end

% 17 - Finding Perfect Squares

function b = isItSquared(a)
    b = false;
    for i = 1:length(a)
        n = a(i)^2;
        for j = 1:length(a)
            if n == a(j)
                b = true;
            end
        end
    end
end

% 18 - Remove any Row in which a NaN appears

function B = remove_nan_rows(A)
    [row col] = size(A);
    t = []
    for i = 1:row
        for j = 1:col
            if isnan(A(i,j))
                t = [t,i];
                break
            end
        end
    end
    A(t,:) = [];
    B = A;
end

% 19 - Find the numeric mean of the prime numbers in a matrix

function out = meanOfPrimes(in)
    [row col] = size(in);
    t = 0
    x = 0
    for i = 1:row
        for j = 1:col
            if isprime(in(i,j))
                x = x + 1;
                t = t + in(i,j);
            end
        end
    end
    out = t/x
end

% 20 - Who has the most change?

function b = most_change(a)
    [row col] = size(a);
    t = [];
    for i = 1:row
        x = a(i,1) * 0.25;
        y = a(i,2) * 0.1;
        z = a(i,3) * 0.05;
        xx = a(i,4) * 0.01;
        xy = x + y + z + xx;
        t = [t, xy];
    end
    b = max(t)
end

% 21 - Most nonzero elements in a row

function r = fullest_row(a)
    [row col] = size(a);
    t =[];
    for i = 1:row
        x = 0;
        for j = 1:col
            if a(i,j) == 0
                x = x + 1;
            end
        end
    t = [t, x];
    end
    r = find(t==min(t))
end

% 22 - Return the 3n+1 Sequence for n

function c = collatz(n)
    t = [];
    while n ~= 1
        if mod(n,2) == 0
            n = n/2;
            t = [t,n];
        else
            n = (3*n)+1;
            t = [t,n];
        end
    end
    c = t
end

% 23 - Summing Digits

function b = sumDigits(n)
    b = 0
    x = 2^n
    y = num2str(x)
    for i = 1:length(y)
        b = b + str2num(y(i))
    end
end

% 24 - Remove the Vowels

function s2 = refcn(s1)
    s2 = strrep(s1, 'a', '');
    s2 = strrep(s2, 'e', '');
    s2 = strrep(s2, 'i', '');
    s2 = strrep(s2, 'o', '');
    s2 = strrep(s2, 'u', '');
end

% 25 - Sort a list of Complex Numbers based on how far they are from the
% Origin

function zSorted = complexSort(z)
    t = [];
    for i = 1:length(z)
        x = real(z(i));
        y = imag(z(i));
        p = sqrt(x^2 + y^2);
        t = [t, p];
    end
    zSorted = sort(t, 'descend')
end

% 26 - Which Values Occur Exactly Three Times

function y = threeTimes(x)
    t = [];
    for i = 1:length(x)
        if ~ismember(x(i), t)
            z = length(find(x==x(i)));
            if z == 3
                t = [t,x(i)];
            end
        end
    end
    y = t
end

% 27 - Bullseye Matrix

function a = bullseye(n)
    t = (n+1)/2;
    m(1:n,1:n)= t;
    for i = 1:t
        m(1+i:n-i,1+i:n-i)= t-i;
    end
    a = m
end

% 28 - The Glodbach Conjecture

function [p1,p2] = goldbach(n)
    for i = 1:n
        if isprime(i)
            p1 = n-i
            if isprime(p1)
                p2 = i
                break
            end
        end
    end
end

% 29 - Return the Largest number that i adjacent to a Zero

function y = nearZero(x)
    for i = 1:length(x)
        if x(i) == 0
            if i == 1
                y = x(i+1)
            elseif i == length(x)
                if x(i-1) > y
                    y = x(i-1)
                end
            else
                if x(i-1) > y
                    y = x(i-1)
                end
                if x(i+1) > y
                    y = x(i+1)
                end
            end
        end
    end
end

% 30 - Target Sorting

function b = targetSort(a,t)
    n = 0;
    while n ~= 1
        n = 1;
        for i = 1:length(a)-1
            x = abs(t - a(i))
            y = abs(t - a(i+1));
            if y > x
                z = a(i);
                a(i) = a(i+1);
                a(i+1) = z;
                n = 0;
            end
        end
    end
    b = a;
end

% 31 - Nearest Numbers

function [index1 index2] = nearestNumbers(A)
    y = abs(max(A)-min(A));
    for i = 1:length(A)
        for j = i+1:length(A)
            x = abs(A(i)-A(j));
            if x < y
                y = x;
            end
        end
    end
    for i = 1:length(A)
        for j = i+1:length(A)
            x = abs(A(i)-A(j));
            if x == y
                index1 = find(A==A(i));
                index2 = find(A==A(j));
            end
        end
    end
end

% 32 - Pascal's Triangle

function y = pascalTri(n)
    t =[];
    for i = 1:n
        t = [t,i];
    end
    for j = n-1:-1:1
        t = [t,j];
    end
    if n == 0
        y = 1;
    else
        y = t;
    end    
end

% 33 - Remove all the Consonants

function s2 = refcn2(s1)
    s1(regexp(s1, '[^aeiou]')) = [];
    s2 = s1;
end

% 34 - Find the longest Sequence of 1's in a Binary Sequence

function y = lengthOnes(x)
    y = 0;
    if x(1) == 1
        t = 1;
    else
        t = 0;
    end    
    for i = 1:length(x)-1
        if x(i+1) == 1
            t = t + 1;
            if t > y
                y = t;
            end
        else
            t = 0;
        end
    end
end

% 35 - Cell Joiner

function out_str = cellstr_joiner(in_cell, delim)
    out_str = join(in_cell, delim);        
end

% 36 - Find the Alphabetic word product

function p = word_product(s)
    p = 1
    alph = 'abcdefghijklmnopqrstuvwxyz';
    alphcap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    for i = 1:length(s)
        if contains(alph, s(i))
            for j = 1:length(alph)
                if s(i) == alph(j)
                    p = p * j
                end
            end
        elseif contains(alphcap, s(i))
            for k = 1:length(alph)
                if s(i) == alphcap(k)
                    p = p * k
                end
            end
        end
    end
end

% 37 - Pangrams!

function tf = isPangram(s)
    t = '';
    s = upper(s);
    for i = 1:length(s)
        if ~contains(t, s(i))
            t = [t,s(i)];
        end
    end
    tf = length(t) == 26;
end

% 38 - Binary Numbers
% To do!

function A = binary_numbers(n)
    A = dec2bin(n,3)
end

% 39 - Counting Money

function b = moneySum(a)
    b = 0;
    for i = 1:length(a)
        a(i) = erase(a(i), '$');
        a(i) = erase(a(i), ',');
        str = cell2mat(a(i));
        num = str2num(str);
        b = b + num;
    end
end

% 40 - Which Doors are Open?
% Not good, find how to loop over the index

function y = which_doors_open(n)
    a = false;
    b = false;
    c = false;
    t = [a, b , c];
    for i = 1:n
        if t(i) == true
            t(i) = false;
        else
            t(i) = true;
        end
        for j = 1:n
            if t(i+(j*n)) == true
                t(i+(j*n)) = false;
            else
                t(i+(j*n)) = true;
            end
        end
    end
end

% 41 - Replace NaNs with the number that appears to its left in the row

function y = replace_nans(x)
    y = [];
    if isnan(x(1))
        x(1) = 0;
    end
    for i = 1:length(x)
        if i == 1
            if isnan(x(i))
                y = [0];
            end
        end
        if isnan(x(i))
            y = [y,x(i-1)];
        else
            y = [y,x(i)];
        end
    end
end

% 42 - Making Change

function b = makingChange(a)
    b = [];
    spitta = [100 50 20 10 5 2 1 0.5 0.25 0.1 0.05 0.01];
    x = a
    for i = 1:length(spitta)
        if x > spitta(i)
            n = fix(x/spitta(i));
            b = [b,n];
            x = mod(x,spitta(i));
        else
            b = [b,0];
        end
    end
end

% 43 - Balanced Number

function tf = isBalanced(n)
    nstr = num2str(n);
    niseven = mod(length(nstr),2) == 0;
    first = 0;
    last = 0;
    x = fix(length(num2str(n))/2);
    for i = 1:x
        first = first + str2num(nstr(i));
        if niseven == true
            last = last + str2num(nstr(x+i));
        else
            last = last + str2num(nstr(x+1+i));
        end
    end
    tf = first == last;
end

% 44 - Quote Doubler

function s2 = quote_doubler(s1)
    s2 = replace(s1, "'", "''");
end

% 45 - De-Dupe

function b = dedupe(a)
    b = []
    for i = 1:length(a)
        if ~ismember(a(i),b)
            b = [b, a(i)]
        end
    end
end

% 46 - Elapsed Time

function elapsed = elapsed_time(d1,d2)
    m = datetime(d1, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
    n = datetime(d2, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
    elapsed = n-m;
end

% 47 - Reverse Run-Length Encoder

function y = RevCountSeq(x)
    y =[];
    if mod(length(x),2)==0
        for i = 1:2:length(x)
            for j = 1:x(i)
                y = [y,x(i+1)];
            end
        end
    end
end

% 48 - Return a List by Number of Occurrences
% Not good enough, need sorting

function y = popularity(x)
    y = unique(x);
end

% 49 - Trimming Spaces

function b = removeSpaces(a)
    b = strtrim(a);
end

% 50 - Find the two most distant points

function ix = mostDistant(p)
    x = 0;
    [row col] = size(p);
    for i = 1:row
        for j = i+1:row
            a = sqrt((p(i,1)-p(j,1))^2 + (p(i,2)-p(j,2))^2)
            if a > x
                x = a;
                ix = [i,j];
            end
        end
    end
end















