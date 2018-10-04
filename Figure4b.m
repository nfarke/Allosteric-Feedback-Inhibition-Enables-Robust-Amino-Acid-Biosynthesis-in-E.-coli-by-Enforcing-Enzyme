function Figure4b
load PAR.mat
%e1
e1_wt = log2(par_wt(end-1,:))-log2(median(par_wt(end-1,:),2));
e1_af = log2(par_af(end-1,:))-log2(median(par_wt(end-1,:),2));

%e2
e2_wt = log2(par_wt(end,:))-log2(median(par_wt(end,:),2));
e2_af = log2(par_af(end,:))-log2(median(par_wt(end,:),2));

%met1
met1_wt = log2(par_wt(end-3,:))-log2(median(par_wt(end-3,:),2));
met1_af = log2(par_af(end-3,:))-log2(median(par_wt(end-3,:),2));

%met2
met2_wt = log2(par_wt(end-2,:))-log2(median(par_wt(end-2,:),2));
met2_af = log2(par_af(end-2,:))-log2(median(par_wt(end-2,:),2));

%define colors
black = [0 0 0];
blue = [100/255 190/255 255/255];
col = [black;blue;black;blue;black;blue;black;blue];


E1 = [e1_wt;e1_af;e2_wt;e2_af;met1_wt;met1_af;met2_wt;met2_af];
boxplot(E1','labels',{'e1_complete_model','e1_only_tf_model','e2_complete_model','e2_only_tf_model','m1_complete_model','m1_only_tf_model','m2_complete_model','m2_only_tf_model',},'widths',0.4,'color',col);
set(gca,'xticklabelrotation',45)
set(gca,'fontsize',12);
ylabel('relative concentration (log2)');
h = findobj(gca,'tag','Outliers');
delete(h)
h  =  findobj(gca,'tag','Median');
set(h,'Color','k')
end