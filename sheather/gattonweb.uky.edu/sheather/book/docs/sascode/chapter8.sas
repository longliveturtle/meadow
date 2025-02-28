data mich;
 infile 'data/MichelinFood.txt' firstobs=2 
     expandtabs;
 input rating in notin mi prop;
run;
quit;

*table 8.1;
proc print data=mich;
run;
quit;

*fig. 8.1;
goptions reset = all;
symbol1 v=circle c=black h=1;
 axis1 label = (h=2 font=times angle=90 
     "Sample proportion") value=(font=times 
     h=1) order=(0 to 1 by .2) offset=(2,2);
 axis2 label = (h=2 font=times 
     'Zagat Food Rating') value =(font=times 
     h=1) offset=(2,2) ;
proc gplot data = mich;
 plot prop*rating /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;
quit;

*p. 267, null deviance;
proc logistic data = mich;
 model in/mi= /scale=D;
run;
quit;

*p. 267, residual deviance, p. 274;
proc logistic data = mich;
 model in/mi=rating /scale=D;
run;
quit;


data makex;
 do x = 15 to 28 by 0.1;
   output;
 end;
 run;
 quit;

data makex;
 set makex;
 eqn = -10.8415 + 0.5012*x;
 phat = 1/(1+exp(-eqn));
 drop eqn;
run;
quit;

data join;
 set makex mich;
run;
quit;

*fig. 8.2;
goptions reset = all;
symbol1 v=circle c=black h=1;
symbol2 c=black i=join;
 axis1 label = (h=2 font=times angle=90 
     "Probability of inclusion in Michelin Guide")
     value=(font=times h=1) order=(0 to 1 by .2)
     offset=(2,2);
 axis2 label = (h=2 font=times
     'Zagat Food Rating') value =(font=times 
     h=1) offset=(2,2) ;
proc gplot data = join;
 plot prop*rating=1 phat*x=2 /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2 overlay;
run;
quit;

proc logistic data = mich noprint;
 model in/mi=rating;
 output out=probs p=predprob;
run;
quit;

data table8p2;
 set probs;
 odds = predprob/(1-predprob);
run;

*table 8.2;
proc print data=table8p2;
var rating predprob odds;
run;
quit;

proc logistic data = mich noprint;
model in/mi=rating;
output out=resdat p=predprob reschi=pearson resdev=devres;
run;
quit;
 
data resdata;
set resdat;
respres = prop-predprob;
run;
quit;

*table 8.3;
proc print data=resdata;
var rating prop predprob respres pearson devres;
run;
quit;

data resdata;
set resdata;
pearson2 = pearson**2;
devres2 = devres**2;
run;
quit;

proc means data= resdata noprint;
var pearson pearson2 devres devres2;
output out =resdatsum SUM=;
run;
quit; 

proc print data=resdatsum;
run;
quit;

data resdatsum;
set resdatsum;
vardev = (devres2 - devres**2/14)/(14-1);
varpear = (pearson2 - pearson**2/14)/(14-1);
stddev = sqrt(vardev);
stdpear = sqrt(varpear);
key = 1;
keep key stddev stdpear;
run;
quit;

proc sort data=resdatsum;
by key;
run;
quit;

data fig83;
set resdat;
key = 1;
run;
quit;

proc sort data=resdatsum;
by key;
run;
quit;

data fig83;
merge fig83 resdatsum;
by key;
run;
quit;

proc print data=fig83;
run;
quit;

data fig83;
set fig83;
stdzpearson=pearson/stdpear;
stdzdev=devres/stddev;
run;
quit;

*fig. 8.3;
goptions reset = all;
symbol1 v=circle c=black;
 axis1 label = (h=2 font=times angle=90 
     "Standardized Deviance Residuals")
     value=(font=times h=1) offset=(2,2);
 axis2 label = (h=2 font=times 'Food Rating') 
     value =(font=times h=1) offset=(2,2) ;
proc gplot data = fig83;
 plot stdzdev*rating /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2 overlay;
run;
quit;

goptions reset = all;
symbol1 v=circle c=black;
 axis1 label = (h=2 font=times angle=90 
     "Standardized Pearson Residuals")
     value=(font=times h=1) offset=(2,2);
 axis2 label = (h=2 font=times 'Food Rating') 
     value =(font=times h=1) offset=(2,2) ;
proc gplot data = fig83;
 plot stdzpearson*rating /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2 overlay;
run;
quit;

proc import datafile='data/MichelinNY.csv' out=michny;
 getnames=yes;
 run;
 quit;

data jitter;
 set michny;
 if inmichelin = 0 then jin = InMichelin+uniform(0)*.03;
 else jin = InMichelin-uniform(0)*.03;
 jfood = food+(uniform(0)*.3-.15);
run;
quit;

goptions reset = all;
symbol1 v=circle c=black h=1;
 axis1 label = (h=2 font=times angle=90 
     "In Michelin Guide? (0=No, 1=Yes)") 
     value=(font=times h=1) order=(0 to 1 by .2)
     offset=(2,2);
 axis2 label = (h=2 font=times 
     'Food Rating') value =(font=times 
     h=1) offset=(2,2);
proc gplot data = jitter;
 plot jin*jfood /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;
quit;

%macro boxplot(var =,data =,yvar =, xvarlab=,yvarlab=);
goptions reset = all;
proc sort data = &data;
 by &var;
run;
quit;
proc gplot data = &data;
 axis1 label=(f=times h=2 "&xvarlab") minor=none 
    order=(-1 0.2 0.8 2) value=(f=times h=1 
    t=1 ' ' t=2 '0' t=3 '1' t=4 ' ');
 axis2 label=(f=times h=2 angle=90 "&yvarlab") 
    value=(f=times h=1);
 symbol1 value = circle interpol=boxt bwidth=36;
 plot &yvar*&var/ haxis=axis1 vaxis=axis2 vminor=0;
run;
quit;
%mend;                                            
%boxplot(var=inmichelin,data=michny,yvar=food,yvarlab=Food Rating,xvarlab=In Michelin Guide? (0 = No, 1 = Yes));


*p. 279;
proc logistic data = michny;
 model inmichelin(event='1')=food;
run;
quit;


proc logistic data = michny noprint;
 model inmichelin(event='1')=food;
 output out=output p=preds reschi=pearson 
     resdev=devres h=hat;
run;
quit;

data z_scores;
 set output;
   stdzdd=devres/sqrt(1-hat);
   stdzdp=pearson/sqrt(1-hat);
 keep stdzdd stdzdp food;
run;
quit;

*fig. 8.6;
goptions reset = all;
symbol1 v=circle h=1;
axis1 value=(f=times h=1) label=(h=2 f=times 
     angle=90 "Standardized Deviance Residuals");
axis2 value=(f=times h=1)
     label=(f=times h=2 "Food Rating") offset=(2,2);
proc gplot data = z_scores;
 plot stdzdd*food/vaxis=axis1 
     hminor=0 haxis=axis2 vminor=0;
run; 
quit;
axis1 value=(f=times h=1) label=(h=2 f=times 
     angle=90 "Standardized Pearson Residuals");
proc gplot data = z_scores;
 plot stdzdp*food/vaxis=axis1 
     hminor=0 haxis=axis2 vminor=0;
run; 
quit;

proc loess data = michny;
 model inmichelin=food/smooth=0.66666667;
 ods output OutputStatistics=loess;
run;

data join;
set loess jitter makex;
run;
quit;

proc sort data=join;
by food x;
run;
quit;

*fig. 8.7;
goptions reset = all;
symbol1 v=circle c=black h=1;
symbol2 c=blue i=join l=1;
symbol3 c=red i=join l=2;
 axis1 label = (h=2 font=times angle=90 
     "In Michelin Guide? (0=No, 1=Yes)")
     value=(font=times h=1) order=(0 to 1 by .2)
     offset=(2,2);
 axis2 label = (h=2 font=times 
     'Food Rating') value =(font=times h=1) 
     offset=(2,2) ;
proc gplot data = join;
 plot jin*jfood=1 phat*x=2 pred*food=3/hminor=0 
     vminor=0 vaxis=axis1 haxis=axis2 overlay;
run;
quit;

*fig. 8.8;
%boxplot(var=inmichelin,data=michny,yvar=food,yvarlab=Food Rating,xvarlab=In Michelin Guide? (0 = No, 1 = Yes));
%boxplot(var=inmichelin,data=michny,yvar=decor,yvarlab=Decor Rating,xvarlab=In Michelin Guide? (0 = No, 1 = Yes));
%boxplot(var=inmichelin,data=michny,yvar=service,yvarlab=Service Rating,xvarlab=In Michelin Guide? (0 = No, 1 = Yes));
%boxplot(var=inmichelin,data=michny,yvar=cost,yvarlab=Price Rating,xvarlab=In Michelin Guide? (0 = No, 1 = Yes));

data michny;
 set michny;
 lcost=log(cost);
run;
quit;

proc logistic data = michny noprint;
 model inmichelin(event='1') =food decor service cost lcost;
 output out=points p=predprob;
run;
quit;


proc loess data = points;
 model predprob=food/smooth=0.66666667;
 ods output OutputStatistics=loesspred;
run;
quit;

*fig. 8.9;
data join;
set loess michny;
run;
quit;

proc sort data=join;
by food;
run;
quit;

goptions reset = all;
symbol1 v=circle c=black h=1;
symbol2 c=black i=join l=1;
 axis1 label = (h=2 font=times angle=90 
     "Y, In Michelin Guide? (0=No, 1=Yes)")
     value=(font=times h=1) order=(0 to 1 by .2)
     offset=(2,2);
 axis2 label = (h=2 font=times 
     'Food Rating, x1') value =(font=times h=1) 
     offset=(2,2) ;
proc gplot data = join;
 plot inmichelin*food=1 pred*food=2/hminor=0 vminor=0
     vaxis=axis1 haxis=axis2 overlay;
run;
quit;

data join2;
set points loesspred;
run;
quit;

proc sort data=join2;
by food;
run;
quit;

goptions reset = all;
symbol1 v=circle c=black h=1;
symbol2 c=black i=join l=1;
 axis1 label = (h=2 font=times angle=90 
     "Y, In Michelin Guide? (0=No, 1=Yes)")
     value=(font=times h=1) order=(0 to 1 by .2)
     offset=(2,2);
 axis2 label = (h=2 font=times 
     'Food Rating, x1') value =(font=times h=1) 
     offset=(2,2) ;
*variable preds are y-hats, variable pred is loess
	 fit;
proc gplot data = join2;
 plot predprob*food=1 pred*food=2/hminor=0 vminor=0
     vaxis=axis1 haxis=axis2 overlay;
run;
quit;


%macro logitmmplot(dsn=, yvar=, x1=, allx=, xlab=, ylab=);
proc loess data = &dsn;
 model &yvar=&x1/smooth=0.6667;
 ods output OutputStatistics=loessout1;
run;
quit;

proc logistic data = &dsn noprint;
 model &yvar(ref='0')= &allx;
 output out = outreg p=yhat;
run;
quit;

proc loess data = outreg;
 model yhat=&x1/smooth=0.6667;
 ods output OutputStatistics=loessout;
run;
quit;

data loessout2;
 set loessout;
 Pred2 = pred;
 drop pred;
run;
quit;

data fit;
 set outreg;
 set loessout1;
 set loessout2;
run;
quit;

proc sort data = fit;
 by &x1;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times "&xlab") 
     value=(f=times h=1);
 axis2 label=(h=2 angle=90 f=times
     "&ylab") order=(0 to 1 by 1) value=(f=times h=1);
symbol1 v = circle c=black;
symbol2 i=join c=blue;
symbol3 i=join c=red l=2;
proc gplot data = fit;
 plot /*points:*/ &yvar*&x1=1 /*loess:*/
     Pred*&x1=2 Pred2*&x1=3/ haxis=axis1 
     vaxis=axis2 vminor=0 hminor=0 overlay;
run;
quit;
%mend;

*fig. 8.10;
%logitmmplot(dsn=michny, yvar=inmichelin, x1=food, allx=food 
     decor service cost lcost,xlab=Food, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=decor, allx=food 
     decor service cost lcost, xlab=Decor, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=service, allx=food 
     decor service cost lcost, xlab=Service, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=cost, allx=food 
     decor service cost lcost, xlab=Price, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=lcost, allx=food 
     decor service cost lcost, xlab=log(Price), ylab=y);

data points;
set points;
linpred= -log((1-predprob)/predprob);
run;
quit;
%logitmmplot(dsn=points, yvar=inmichelin, x1=linpred, allx=food 
     decor service cost lcost, xlab=Linear Predictor, ylab=y);


%macro respdifplot(dsn=, yvar=, xvar=, 
     ylab=, xlab=, byvar=);
 data zero;
  set &dsn;
  if inmichelin=0;
  xvarzero = &xvar;
 run;
 quit;

 data one;
  set &dsn;
  if inmichelin=1;
  xvarone = &xvar;
 run;
 quit;

 data plotit;
  set zero one;
 run;
 quit;

 goptions reset = all;
 axis1 label=(h=2 f=times angle=90 "&ylab") 
     value=(h=1 f=times) offset=(1,1);
 axis2 label=(h=2 f=times "&xlab") v=(h=1
     f=times) offset=(1,1);
 symbol1 v=circle i=r h=1 c=black l=1;
 symbol2 v=triangle h=1 i=r c=red l=2;
 legend1 across=1 frame /*offset=(12 pct, -1 pct)*/
     position=(bottom right inside) label=
	 (h=1 f=times position=top
     'In Michelin Guide?') value=(h=1 
     f=times 'No' h=1 f=times 'Yes');
 proc gplot data = plotit;
  plot &yvar*xvarzero=1 &yvar*xvarone=2/ 
     vaxis=axis1 haxis=axis2 overlay
     vminor=0 hminor=0 legend=legend1;
 run;
 quit;
 %mend;

 *fig. 8.11;
 %respdifplot(dsn=michny, yvar=service, xvar=decor,
     ylab=Service Rating, xlab=Decor Rating,
	 byvar=inmichelin);

data michny;
set michny;
decserv =decor*service;
run;
quit;

*fig. 8.12;
%logitmmplot(dsn=michny, yvar=inmichelin, x1=food, allx=food 
     decor service cost lcost decserv,xlab=Food, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=decor, allx=food 
     decor service cost lcost decserv, xlab=Decor, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=service, allx=food 
     decor service cost lcost decserv, xlab=Service, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=cost, allx=food 
     decor service cost lcost decserv, xlab=Price, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=lcost, allx=food 
     decor service cost lcost decserv, xlab=log(Price), ylab=y);

proc logistic data = michny noprint;
 model inmichelin(ref='0') =food decor service cost lcost decserv;
 output out=points p=predprob;
run;
quit;

data points;
set points;
linpred= -log((1-predprob)/predprob);
run;
quit;

%logitmmplot(dsn=points, yvar=inmichelin, x1=linpred, allx=food 
     decor service cost lcost decserv, xlab=Linear Predictor, ylab=y);

*p. 290, model 1;
proc logistic data = michny;
 model inmichelin(ref='0') =food decor service cost lcost;
run;
quit;

*p. 290, model 2, p.291-292;
proc logistic data = michny;
 model inmichelin(ref='0') =food decor service cost lcost decserv;
run;
quit;

*p. 290 p-value;
data devpval;
pval = 1-probchi(136.431-129.820,1);
run;
quit;
proc print data=devpval;
run;
quit;


proc logistic data = michny noprint;
 model inmichelin(ref='0') =food decor service cost lcost decserv;
 output out=points resdev=devres h=hat;
run;
quit;

data points;
set points;
  stdres = devres/sqrt(1-hat);
run;
quit;

*fig. 8.13;
goptions reset = all;
symbol1 v=circle c=black h=1;
 axis1 label = (h=2 font=times angle=90 
     "Standardized Deviance Residuals")
     value=(font=times h=1)
     offset=(2,2);
 axis2 label = (h=2 font=times 
     'Leverage Values') value =(font=times h=1) 
     offset=(2,2) ;
proc gplot data = points;
 plot stdres*hat=1 /hminor=0 vminor=0
     vaxis=axis1 haxis=axis2 href=.085;
run;
quit;

*p. 292-293;
proc logistic data = michny;
 model inmichelin(ref='0') =food decor service lcost decserv;
run;
quit;

*p. 292 p-value;
data devpval;
pval = 1-probchi(131.229-129.820,1);
run;
quit;
proc print data=devpval;
run;
quit;


*fig. 8.14;
%logitmmplot(dsn=michny, yvar=inmichelin, x1=food, allx=food 
     decor service lcost decserv,xlab=Food, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=decor, allx=food 
     decor service lcost decserv, xlab=Decor, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=service, allx=food 
     decor service lcost decserv, xlab=Service, ylab=y);
%logitmmplot(dsn=michny, yvar=inmichelin, x1=lcost, allx=food 
     decor service lcost decserv, xlab=log(Price), ylab=y);

proc logistic data = michny noprint;
 model inmichelin(ref='0') =food decor service lcost decserv;
 output out=points p=predprob resdev=devres h=hat;
run;
quit;

data points;
set points;
linpred= -log((1-predprob)/predprob);
run;
quit;

%logitmmplot(dsn=points, yvar=inmichelin, x1=linpred, allx=food 
     decor service lcost decserv, xlab=Linear Predictor, ylab=y);

data points;
set points;
  stdres = devres/sqrt(1-hat);
run;
quit;

*fig. 8.15;
goptions reset = all;
symbol1 v=circle c=black h=1;
 axis1 label = (h=2 font=times angle=90 
     "Standardized Deviance Residuals")
     value=(font=times h=1)
     offset=(2,2);
 axis2 label = (h=2 font=times 
     'Leverage Values') value =(font=times h=1) 
     offset=(2,2) ;
proc gplot data = points;
 plot stdres*hat=1 /hminor=0 vminor=0
     vaxis=axis1 haxis=axis2 href=.073;
run;
quit;

data table8p5;
set points;
if abs(stdres) < 2 then delete;
drop devres hat lcost _LEVEL_ linpred decserv;
run;
quit;

*table 8.5;
proc print data= table8p5;
run;
quit;
 
