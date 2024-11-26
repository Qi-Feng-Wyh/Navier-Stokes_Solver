x=[6,7,8,9,10,11];
y=[1.495089e-03,3.736385e-04,9.340787e-05,2.335888e-05,5.848799e-06,1.477849e-06];
y=log2(y);
scatter(x,y,'markerfacecolor',[0,0,1]);
xlabel('log_2 N');
ylabel('log_2 e_N');
p=polyfit(x,y,1);
fprintf('k=%e\n', p(1));