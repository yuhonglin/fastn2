%% configuration
d=0.101;
N = 10;
nTest = 20000;

%% generate data

S_list = cell(nTest);
mu_list = cell(nTest);

for i = 1:nTest
    H = rand(100,N);
    S_list{i} = H'*H/100;
    mu_list{i} = mean(H);
end

%% test speed
tic
for i = 1:nTest
    old = maxsrn2(mu_list{i}',S_list{i},d);
end
toc

tic
for i = 1:nTest
    new = fastn2mex(mu_list{i}',S_list{i},d);
end
toc

%% test correctness
difference = nan*ones(1,nTest);
for i = 1:nTest
    old = maxsrn2(mu_list{i}',S_list{i},d);
    new = fastn2mex(mu_list{i}',S_list{i},d);
    difference(i) = norm(old-new);
end

%[old, new]

%[mu*old/sqrt(old'*S*old) mu*new/sqrt(new'*S*new)]

%[d-sum(old.^2) d-sum(new.^2)]