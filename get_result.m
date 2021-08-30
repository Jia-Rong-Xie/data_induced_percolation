clear;
clc;

field_id = 4;
switch field_id
    case 1
        name = 'chaos';
    case 2
        name = 'pt';
    case 3
        name = 'EPLDS';
    case 4
        name = 'nano';
end
data = dlmread(['collaboration_',name,'.txt']);

%% set parameters
if field_id == 4
    year_0 = 2008; 
    year_1 = 2009 : 2013;
    year_2 = 2011 : 2013;
    year_3 = 2014 : 2016;
    thre_pub = 5;
    thre_focus = 5;
    thre_e = 2;
    flap_num_aut = 20;
    data = [data; dlmread('collaboration_nano_2.txt')];
else
    year_0 = 1998;  % the zeroth year of the data, i.e. the first year - 1
    year_1 = 1999 : 2003;
    year_2 = 2001 : 2003; 
    year_3 = 2004 : 2006;
    thre_pub = 2; % threshold of publication, the nodes are the scientists who have published at least # papers
    thre_focus = 2; % threshold of focused authors
    thre_e = 1; % threshold of collaboration, a link is constructed between any pair of scientists who have at least # joint publication
    flap_num_aut = 10;
end

%% collaboration network construction
dd = find( data(:,2) <= year_3(end) );
data = data(dd,:);
dd = find( data(:,2) <= year_1(end) );
data_2 = data(dd,:);
C = sparse(max(data_2(:,5)), max(data_2(:,5)));
up_flap = 0;
bot_flap = 0;
for p_id = 1:length(unique(data_2(:,1)))
    up_flap = bot_flap + 1;
    bot_flap = up_flap + data_2(up_flap,4)-1;
    auts = data_2(up_flap:bot_flap,5);
    for aut_1 = auts'
        for aut_2 = auts'
            C(aut_1,aut_2) = C(aut_1,aut_2)+1;
        end
    end
end
A = sparse(size(C,1), size(C,1));
A(find(C>=thre_e)) = 1; 
A = A - diag(diag(A));

C = sparse(data_2(:,5), data_2(:,2)-year_0, 1, max(data_2(:,5)), 5);
auts = find( sum(C,2) >= thre_pub);
A = A(auts, auts);
dd = find(sum(A)>0);
auts = auts(dd);
A = A(dd, dd); % the final collaboration network

%% count the publication of old and new fields for each authors
data = data(find( ismember(data(:,5), auts)), :); 
C = sparse(data(:,5), data(:,2)-year_0, 1, max(data(:,5)), 8);
C = C(auts,:);
pub_sum = full(sum(C,2));
pub_sum_year2 = full(sum(C(:,3:5),2));
pub_sum_year3 = full(sum(C(:,6:8),2));

dd = find(data(:,3) == 1);
C = sparse(data(dd,5), data(dd,2)-year_0, 1, max(data(:,5)), 8);
C = C(auts,:);
pub_old = full(sum(C,2));
pub_old_year2 = full(sum(C(:,3:5),2));
pub_old_year3 = full(sum(C(:,6:8),2));

dd = find(data(:,3) == 2);
C = sparse(data(dd,5), data(dd,2)-year_0, 1, max(data(:,5)), 8);
C = C(auts,:);
pub_new = full(sum(C,2));
pub_new_year2 = full(sum(C(:,3:5),2));
pub_new_year3 = full(sum(C(:,6:8),2));

%% define the focused authors, and count the induced index
focus_ornot = (pub_new_year2 == 0) & (pub_old_year2 >= thre_focus);
focus_auts = find(focus_ornot);

m_focus = full(sum(A(focus_auts,focus_auts))); 
temp_mat = A(:,focus_auts);
temp_mat = temp_mat.*repmat(m_focus, length(auts), 1);
m_list = full(max(temp_mat,[],2)); % the induced index of each node
m_list(focus_auts) = m_list(focus_auts) - 1; % minus 1: exclude the node itself
m_list(focus_auts) = - 1; % the focused authors are not considered
m_list(find(sum(A(:,focus_auts),2) <= 0)) = -2; % for the nodes without focused neighbors, the induced index are not defined (marked by -2)

%% observation of the induced percolation phenomenon
dm = 3;
m_max = max(m_list);
mm_flaps = [(0:dm:m_max)',(0:dm:m_max)'+dm-1]; % induced index range
res_m = []; % the final result
for mm_index = 1:size(mm_flaps,1)
    thought_aut = find( (m_list >= mm_flaps(mm_index,1)) & (m_list <= mm_flaps(mm_index,2)) & ( pub_sum_year3 > 0) );
    if length(thought_aut) >= flap_num_aut
        res_ratio = mean(pub_old_year3(thought_aut)./pub_sum_year3(thought_aut));
        res_m = [res_m, res_ratio];
    else
        break
    end
end
