data anscombe;
 infile 'data/anscombe.txt' firstobs=2;
 input case x1-x4 y1-y4;
run;

*Fig. 3.1;
%macro fourplots;
%do i = 1 %to 4;
 goptions reset = all;
 proc gplot data = anscombe;
  axis1 label=(h=2 font=times angle=90 "y&i") 
     order=(2 to 14 by 2) value=(font=times h=1);
  axis2 label=(h=2 font=times "x&i") 
     order=(0 to 20 by 5) value=(font=times h=1);
  symbol1 value=circle interpol=r;
  plot y&i*x&i / vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
 run;
 quit;
%end;
%mend;
%fourplots;


%macro reg;
 proc reg data = anscombe;
  %do i = 1 %to 4;
   model y&i = x&i;
  %end;
 run;
%mend;
%reg;


/*Fig. 3.2:  First get the residuals from
proc reg, and then plot them against x.*/
%macro regout;
 proc reg data = anscombe noprint;
  %do i = 1 %to 4;
   model y&i = x&i;
   output out = resdata&i r=resids;
  %end;
 run;
%mend;
%regout;

%macro resplots;
%do i = 1 %to 4;
 goptions reset = all;
 proc gplot data = resdata&i;
  title1 height=2 "Data Set &i" ;
  axis1 label=(font=times h=2 angle=90 "Residuals") 
     order=(-2 to 2 by 1) value=(font=times h=1);
  axis2 label=(font=times h=2 "x&i") 
     order=(0 to 20 by 5) value=(font=times h=1);
  symbol1 value=circle;
  plot resids*x&i / vaxis=axis1 haxis=axis2 
     vminor=0 hminor=0;
 run;
 quit;
%end;
%mend;
%resplots;


*Fig. 3.3;
goptions reset = all;
proc gplot data = anscombe;
  axis1 label=(h=2 font=times angle=90 "y2") 
     order=(3 to 10 by 1) value=(font=times h=1);
  axis2 label=(h=2 font=times "x2") 
     order=(4 to 14 by 2) value=(font=times h=1) 
     offset=(2,2);
  symbol1 value=circle interpol=r;
  plot y2*x2 / vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;
quit;

goptions reset = all;
proc gplot data = resdata2;
  axis1 label=(font=times h=2 angle=90 
     "Residuals") order=(-2 to 2 by 1) 
     value=(font=times h=1);
  axis2 label=(font=times h=2 "x2") 
     order=(4 to 14 by 2) value=(font=times h=1) 
     offset=(2,2);
  symbol1 value=circle;
  plot resids*x2 / vaxis=axis1 haxis=axis2 
     vminor=0 hminor=0;
run;
quit;

data huber;
 infile 'data/huber.txt' firstobs=2 expandtabs;
 input x ybad ygood;
run;


proc reg data = huber;
 model ybad=x;
 output out= resbad r=badres;
 model ygood=x;
 output out = resgood r=goodres;
 run;
proc print data = resbad;
 var badres;
proc print data = resgood;
 var goodres;
run;


goptions reset = all;
proc gplot data = huber;
 symbol1 value=circle interpol=r;
 axis1 label=(h=2 font=times angle=90 "YBad") 
     order=(-12 to 3 by 3) value=(f=times h=1);
 axis2 label=(h=2 f=times "x") value=(f=times h=1) 
     order=(-4 to 10 by 2) offset=(2,2);
 plot ybad*x/ vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;
quit;

goptions reset = all;
proc gplot data = huber;
 symbol1 value=circle interpol=r;
 axis1 label=(f=times h=2 angle=90 "YGood") 
     order=(-12 to 3 by 3) value=(f=times h=1);
 axis2 label=(f=times h=2 "x") offset=(2,2)
     value=(f=times h=1) order=(-4 to 10 by 2);
 plot ygood*x/ vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;
quit;


proc reg data = huber noprint;
 model ybad=x;
 output out= levbad h=badlev;
 model ygood=x;
 output out = levgood h=goodlev;
run;
proc print data = levbad;
 var x badlev;
proc print data = levgood;
 var x goodlev;
run;


* 3.8;
goptions reset = all;
proc gplot data = huber;
 symbol1 interpol=rq value=circle color=black;
 axis1 label=(h=2 f=times angle=90 "YBad") 
     order=(-4 to 3 by 1) value=(f=times h=1);
 axis2 label=(f=times h=2 "x") offset=(2,2)
     value=(f=times h=1) order=(-4 to 10 by 2);
 plot ybad*x/vaxis=axis1 haxis=axis2 vminor=0 hminor=0;
run;
quit;

data final; 
 set huber;
 xq = x**2;
proc reg data = final;
 model ybad = x xq;
run;

data bonds;
 infile 'data/bonds.txt' firstobs=2 expandtabs;
 input case couponrate bidprice;
run;


proc reg data = bonds;
 model bidprice=couponrate/clb;
 output out=table3_4 h=leverage r=residuals 
     student=StdResids;
run;

*Table 3.4;
proc print data = table3_4;
 var couponrate bidprice leverage residuals stdresids;
run;

*Fig. 3.9;
goptions reset = all;
proc gplot data = bonds;
 symbol1 value=circle interpol=r;
 axis1 label=(h=2 font=times angle=90 
     "Bid Price ($)") order=(85 to 120 by 5) 
     value=(font=times h=1);
 axis2 label=(h=2 font=times "Coupon Rate (%)") 
     order=(2 to 14 by 2) value=(font=times h=1);
 plot bidprice*couponrate/ vaxis=axis1 haxis=axis2 
     vminor=0 hminor=0;
run;
quit;


*Fig. 3.10;
data flower;
 set table3_4;
 if stdresids < 1.8 && stdresids > -2 then delete;
 labely = stdresids;
 drop stdresids;
run;

data notflower;
 set table3_4;
 if stdresids > 1.8 then delete;
   else if stdresids < -2 then delete;
run;

data labels;
 set flower notflower;
run;

goptions reset = all;
proc gplot data = labels;
 symbol1 value=circle pointlabel= 
     ("#case" height=1 position=bottom) color=black;
 symbol2 value=circle color=black;
 axis1 label=(h=2 angle=90 font=times 
     "Standardized Residuals") order=(-3 to 3 by 1) 
     value=(font=times h=1);
 axis2 label=(h=2 font=times "Coupon Rate (%)") 
     order=(2 to 14 by 2) value=(font=times h=1);
 plot labely*couponrate stdresids*couponrate/overlay 
     vaxis=axis1 haxis=axis2 
            vminor=0 hminor=0 vref=-2 2;
run;


*Fig. 3.11;
data notflower;
 set bonds;
 if case = 4 then delete;
  else if case = 13 then delete;
  else if case = 35 then delete; 
run;

goptions reset = all;
proc gplot data = notflower;
 symbol1 value = circle interpol=r;
 title1 height=2 font=times "Regular Bonds";
 axis1 label=(h=2 font=times angle=90 
     "Bid Price ($)") order=(85 to 120 by 5) 
     value=(font=times h=1);
 axis2 label=(h=2 font=times
     "Coupon Rate (%)") order=(2 to 14 by 2) 
     value=(font=times h=1);
 plot bidprice*couponrate/ vaxis=axis1 
     haxis=axis2 vminor=0 hminor=0;
run;


*Regression output p. 100;
proc reg data = notflower;
 model bidprice=couponrate;
 output out = resids student=stdres;
run;

*Fig. 3.12;
goptions reset = all;
proc gplot data = resids;
 symbol1 value = circle;
 title1 height=2 "Regular Bonds";
 axis1 label=(h=2 font=times angle=90 
     "Standardized Residuals") 
     order=(-4 to 2 by 1) value=(font=times h=1);
 axis2 label=(h=2 font=times "Coupon Rate (%)") 
     order=(2 to 14 by 2) value=(font=times h=1);
 plot stdres*couponrate/ vaxis=axis1 haxis=axis2 
     vminor=0 hminor=0 vref=-2 2;
run;


*Fig. 3.13;
proc reg data= bonds noprint;
 model bidprice=couponrate;
 output out = cook cookd=cd;
run;

data badcook;
 set cook;
 if cd < 0.1212 then delete;
 labely = cd;
 drop cd;
data goodcook;
 set cook;
 if cd > 0.1212 then delete;
data labelcook;
 set badcook goodcook;

goptions reset = all;
proc gplot data = labelcook;
 symbol1 value = circle color=black 
     pointlabel=("#case" height=1 position=bottom);
 symbol2 value=circle color=black;
 axis1 label=(h=2 font=times angle=90 
     "Cook's Distance") order=(0 to 1.2 by 0.2) 
     value=(font=times h=1);
 axis2 label=(h=2 font=times "Coupon Rate (%)") 
     order=(2 to 14 by 2) value=(font=times h=1);
 plot labely*couponrate cd*couponrate/ overlay 
     vaxis=axis1 haxis=axis2 vminor=0 hminor=0 
     vref=0.1212;
run;

data prod;
 infile 'data/production.txt' firstobs=2 expandtabs;
 input case time size;
run;

*Fig. 3.14:  some default graphics are;
ods graphics on;
proc reg data = prod;
 model time=size;
 output out=output r=resids p=fitted student=stdres cookd=cd h=levg;
run;
ods graphics off;


%macro plotlm(regout =,);
proc loess data = &regout;
 model resids=fitted/smooth=0.6667;
 ods output OutputStatistics=loessout;
run;

data fit;
 set regout;
 set loessout;

proc sort data = fit;
 by fitted;
run;

goptions reset = all;
symbol1 v=circle c=black;
symbol2 i=join c=black;
 title1 height=2 font=times  "Residuals vs Fitted";
 axis1 label = (font=times h=2 angle=90 'Residuals') 
      value=(font=times h=1);
 axis2 label = (font=times h=2 'Fitted values') 
      value =(font=times h=1);
proc gplot data = fit;
 plot /*points:*/ resids*fitted=1 /*loess:*/ 
     Pred*fitted=2/ overlay hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2 vref=0;
run;

goptions reset = all htext=1.5;
 title1 height=2 font=times "Normal Q-Q";
 symbol1 value=circle color=black;

proc univariate data = &regout noprint;
 qqplot stdres/normal(mu=0 sigma=1 l=1 color=black) 
     font=times vminor=0 hminor=0 
     vaxislabel= "Standardized Residuals";
run;

data plot3;
 set &regout;
 sqrtres = sqrt(abs(stdres));
run;
proc loess data = plot3;
 model sqrtres=fitted/smooth=0.6667;
 ods output OutputStatistics=loessout;
run;
data fit;
 set plot3;
 set loessout;
proc sort data = fit;
 by fitted;
run;
goptions reset = all;
symbol1 v=circle c=black;
symbol2 i=join c=black;
 title1 height=2 font=times  "Scale-Location";
 axis1 label = (font=times h=2 angle=90 
     'Sqrt(Abs(Res))')  
     value=(font=times h=1);
 axis2 label = (font=times h=2 'Fitted values') 
     value =(font=times h=1);
proc gplot data = fit;
 plot /*points:*/ sqrtres*fitted=1 /*loess:*/ 
     Pred*fitted=2/ overlay hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;

*Fourth plot;
proc sort data = &regout;
 by levg;
 run;
proc loess data = &regout;
 model stdres=levg/smooth=0.67777;
 ods output OutputStatistics=loessout;
run;
data fit;
 set &regout;
 set loessout;
proc sort data = fit;
 by levg;
run;
goptions reset = all;
symbol1 v=circle c=black;
symbol2 i=join c=black;
 title1 height=2 font=times "Residuals vs Leverage";
 axis1 label = (h=2 font=times angle=90 
     "Standardized Residuals")
     value=(font=times h=1);
 axis2 label = (h=2 font=times 'Leverage') 
      value =(font=times h=1) ;
proc gplot data = fit;
 plot /*points:*/ stdres*levg=1 /*loess:*/ Pred*levg=2/
     overlay hminor=0 vminor=0 vaxis=axis1 haxis=axis2 
     vref=0 href=0;
run;
quit;
%mend;

%plotlm(regout=output);


data clean;
 infile 'data/cleaning.txt' firstobs=2 expandtabs;
 input case crews rooms;
run;

*Fig. 3.15;
goptions reset = all;
proc gplot data = clean;
 symbol1 value=circle interpol=r;
 axis1 label=(h=2 angle=90 font=times 
     "Number of Rooms Cleaned") 
     order=(0 to 80 by 10) value=(font=times h=1);
 axis2 label=(h=2 font=times "Number of Crews") 
     order=(2 to 16 by 2) value=(font=times h=1);
 plot rooms*crews/ vaxis=axis1 haxis=axis2 
     vminor=0 hminor=0;
run;

proc reg data = clean; 
 model rooms=crews;
 output out=predintervals p=fit lcl=lwr ucl=upr;
run;

proc print data = predintervals;
 var fit lwr upr;
run;

proc reg data = clean noprint;
 model rooms=crews;
 output out=regout r=resids student=stdres cookd=cd 
     p=fitted h=levg;
run;

goptions reset = all;
proc gplot data = regout;
 symbol1 value=circle;
 axis1 label=(h=2 angle=90 font=times 
     "Standardized Residuals") order=(-2.5 to 2.5 by 1) 
     value=(font=times h=1);
 axis2 label=(h=2 font=times "Number of Crews") 
     offset=(3,3) order=(2 to 16 by 2) 
     value=(font=times h=1);
 plot stdres*crews/ vaxis=axis1 haxis=axis2 
     vminor=0 hminor=0;
run;


*Fig. 3.17;
data plot3;
 set regout;
 sqabs = sqrt(abs(stdres));
run;

goptions reset = all;
proc gplot data = plot3;
 symbol1 value=circle interpol=r;
 axis1 label=(h=2 angle=90 font=times 
     "Sqrt(Abs(Std'zd Resids))") value=(font=times h=1)
     order=(0.2 to 1.6 by 0.2);
 axis2 label=(h=2 font=times "Number of Crews") 
     offset=(3,3) order=(2 to 16 by 2) 
     value=(font=times h=1);
 plot sqabs*crews/ vaxis=axis1 haxis=axis2 
     vminor=0 hminor=0;
run;

%plotlm(regout=regout);

proc sort data = clean;
 by crews;
run;

PROC MEANS STD; 
VAR rooms; 
BY crews;
run;

data stdevs;
 input std @@;
 cards;
3.000000 4.966555 4.690416 6.642665 7.927123 7.28991 12.000463
;
data crews;
 input crews @@;
 cards;
2 4 6 8 10 12 16
;
data plotit;
 set stdevs;
 set crews;
run;

goptions reset = all;
proc gplot data = plotit;
symbol1 v=circle i=r c=black;
 axis1 label = (h=2 font=times angle=90 
     "Stdev(Rooms Cleaned)") 
     order=(1 to 13 by 2) value=(font=times h=1);
 axis2 label = (h=2 font=times 'Number of Crews') 
     order=(2 to 16 by 2) value =(font=times h=1) 
     offset=(2,2);
 plot std*crews=1/hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;

data sqrttrans;
 set clean;
 sqrtcrews = sqrt(crews);
 sqrtrooms = sqrt(rooms);
run;

proc reg data = sqrttrans; 
 model sqrtrooms=sqrtcrews;
 output out=predint lcl=lwr ucl=upr
   r=resids p=fitted student=stdres cookd=cd h=levg;
run;

proc print data = predint;
 var fitted lwr upr;
run;

goptions reset = all;
proc gplot data = sqrttrans;
symbol1 v=circle i=r c=black;
 axis1 label = (h=2 font=times angle=90 "Sqrt(Rooms Cleaned)") 
     order=(2 to 9 by 1) value=(font=times h=1);
 axis2 label = (h=2 font=times 'Sqrt(Number of Crews)') 
     order=(1 to 4 by 0.5) value =(font=times h=1) offset=(2,2);
 plot sqrtrooms*sqrtcrews=1/hminor=0 vminor=0 vaxis=axis1 haxis=axis2;
run;

goptions reset = all;
proc gplot data = predint;
symbol1 v=circle c=black;
 axis1 label = (h=2 font=times angle=90 "Standardized Residuals") 
     order=(-2 to 2.5 by 0.5) value=(font=times h=1);
 axis2 label = (h=2 font=times 'Sqrt(Number of Crews)') 
     order=(1 to 4 by 0.5) value =(font=times h=1) offset=(0,2);
 plot stdres*sqrtcrews=1/hminor=0 vminor=0 vaxis=axis1 haxis=axis2;
run;

%plotlm(regout=predint);

data food;
 infile 'data/confood1.txt' firstobs=2 
     expandtabs;
 input week sales price;
run;

* 3.22;
goptions reset = all;
symbol1 v=circle c=black i=r;
axis1 order=(0 to 7000 by 1000) label=(angle=90 
     f=times h=2 "Sales")
     offset=(0,2) value=(f=times h=1);
axis2 order=(0.55 to 0.85 by 0.05) value=(f=times h=1)
     label=(f=times h=2 "Price");
proc gplot data = food;
 plot sales*price/vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;

data food;
 set food;
 lsales=log(sales);
 lprice=log(price);
run;

* 3.23;
goptions reset = all;
symbol1 v=circle c=black i=r;
axis1 order=(5 to 9 by 1) label=(angle=90 
     f=times h=2 "log(Sales)")
     value=(f=times h=1);
axis2 order=(-0.6 to -0.1 by 0.1) value=(f=times h=1)
     label=(f=times h=2 "log(Price)");
proc gplot data = food;
 plot lsales*lprice/vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;

proc reg data = food;
 model lsales=lprice;
  output out=sres student=stdres;
run;

*3.24;
goptions reset = all;
symbol1 v=circle c=black;
axis1 order=(-2.5 to 3.5 by 1) label=(angle=90 
     f=times h=2 "Standardized Residuals")
     value=(f=times h=1);
axis2 order=(-0.6 to -0.1 by 0.1) value=(f=times h=1)
     label=(f=times h=2 "log(Price)");
proc gplot data = sres;
 plot stdres*lprice/vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;

data transf;
 infile 'data/responsetransformation.txt' 
     firstobs=2 expandtabs;
 input x y;
run;

*3.25;
goptions reset = all;
symbol1 v=circle c=black;
axis1 order=(0 to 100 by 20) value=(f=times h=1)
     label=(h=2 f=times angle=90 "y");
axis2 order=(0 to 5 by 1) value=(f=times h=1)
     label=(f=times h=2 "x");
proc gplot data = transf;
 plot y*x/vaxis=axis1 haxis=axis2 vminor=0 hminor=0;
run;


proc reg data = transf;
 model y=x;
 output out=res student=stdres;
run;

data res;
 set res;
 sqtabsy=sqrt(abs(stdres));
run;

*3.26;
goptions reset = all;
symbol1 v=circle c=black;
axis1 order=(-2 to 5 by 1) value=(f=times h=1)
     label=(h=2 f=times angle=90 
     "Standardized Residuals");
axis2 order=(0 to 5 by 1) value=(f=times h=1)
     label=(f=times h=2 "x");
proc gplot data = res;
 plot stdres*x/vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;

goptions reset = all;
symbol1 v=circle c=black;
axis1 order=(0 to 2.5 by 0.5) value=(f=times h=1)
     label=(h=2 f=times angle=90 
     "Sqrt(Abs(Stdzd Resids))");
axis2 order=(0 to 5 by 1) value=(f=times h=1)
     label=(f=times h=2 "x");
proc gplot data = res;
 plot sqtabsy*x/vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;

%macro kerneld(var =,dat =,) ;
proc kde data = &dat method=sjpi out=dens;
var &var;
run;
goptions reset = all;
symbol1 i=join width=1;
title height=2 font=times 
     "Gaussian kernel density estimate";     
axis1 value=(f=times h=1)
     label=(h=2 f=times angle=90 "Density");
axis2 value=(f=times h=1)
     label=(f=times h=2 "&var");
proc gplot data = dens;
 plot density*&var/vaxis=axis1 haxis=axis2 vminor=0
     hminor=0;
run; 
%mend;

%macro boxplot(dsn=, var=);
 data boxplot;
  set &dsn;
  group=1;
 goptions reset = all;
 symbol1 v=circle  i=boxt bwidth=56;
 axis1 label=(f=times h=2 angle=90 "&var") 
      value=(f=times h=1);
 axis2 order=(1) label=(' ')
      value=(f=times h=1 t=1 ' ');
 proc gplot data = boxplot;
  plot &var*group/vaxis=axis1 haxis=axis2 vminor=0
      hminor=0;
 run;
 quit;
%mend boxplot;

%macro qqplot(var=, dsn=);
 goptions reset = all htext=1.5;
  title1 height=2 font=times "Normal Q-Q";
  symbol1 value=circle color=black;
 proc univariate data = &dsn noprint;
  qqplot &var/normal(mu=est sigma=est l=1 color=black) 
      vminor=0 hminor=0 font=times
      vaxislabel= "&var";
 run;
%mend qqplot;

*3.27;
%kerneld(dat=transf,var=y);
%boxplot(dsn=transf,var=y);
%qqplot(dsn=transf, var=y);
%kerneld(dat=transf,var=x);
%boxplot(dsn=transf,var=x);
%qqplot(dsn=transf, var=x);


%macro irp(resp=, fit=,dat=);
proc iml;

start h_irp_call(lambda) global(resp,fit) ;
Y=log(resp);
if (abs(lambda) > .001) then Y = (1/lambda)#(resp##lambda-J(nrow(resp),1,1));
Y = J(nrow(resp),1,1) || Y ;
predit = Y*(INV(Y`*Y)*Y`*fit) ;
f = (fit - predit)`*(fit - predit) ;
f = -f;
return(f);
finish h_irp_call;

start irp_call(lambda) global(resp,fit)	;
f =  h_irp_call(lambda);
return(f);
finish irp_call;


start unirpopt;
lambda = J(1,1,1);
optn =  j(1,11,.);
optn[1] = 1;
optn[2] = 1;

call nlpqn(rc,lambdares,"irp_call",lambda,optn);
print lambdares;
call symput('lambda',char(lambdares)) ;
finish;

  use &dat;
  read all ;
  resp = &resp ;
  fit = &fit ;
run unirpopt;
quit;

data invresplot;
set &dat;
y=&resp	;
cbrty=y**(&lambda);
ly=log(y) ;
run;


%macro regouts(dsn=,yvar=,predname=);
 proc reg data = invresplot noprint;
  model fitted = &yvar;
  output out= &dsn p=&predname;
 run;
%mend;

%regouts(dsn=data1,yvar=y, predname=lam1hat);
%regouts(dsn=data2,yvar=cbrty,predname=cbrtyhat);
%regouts(dsn=data3,yvar=ly,predname=lyhat);

proc sort data = invresplot;
 by y;
run;

%macro sortit;
%do i = 1 %to 3;
 proc sort data = data&i;
  by y;
 run;
%end;
%mend sortit;
%sortit;

data full;
 merge invresplot data1 data2 data3;
 by y;
run;


goptions reset = all;
symbol1 v=circle c=black;
symbol2 i=join c=black l=5 w=2;
symbol3 i=join c=black l=1 w=2;
symbol4 i=join c=black l=2 w=2;
axis1 label=(h=2 angle=90 f=times "yhat") 
      value=(h=1 f=times);
axis2 label=(h=2 f=times "y")  value=(h=1 f=times) offset=(2,0);
legend1 label=(f=times h=1.5 j=c 'Lambda' position=top) 
     position=(bottom right inside) across=1 frame 
     value=(f=times h=1.5 j=c 'Yhat' j=c '1' 
     j=c "&lambda" j=c '0' j=c);
proc gplot data = full;
 plot &fit*y=1 lam1hat*y=2 cbrtyhat*y=3 lyhat*y=4/
     overlay vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0 legend=legend1;
run;
%mend irp;

proc reg data = transf noprint;
 model y=x;
 output out=regout p=fitted;
run;

%irp(resp=y, fit=fitted,dat=regout);

data powtransf;
set regout;
py1 = 1/y;
py2 = y**(-.5);
py3 = y**(-1/3);
py4 = log(y);
py5 = y**(1/3);
py6 = y**(1/2);
py7 = y;
run;

%macro regit;
%do i = 1 %to 7;
proc reg data=powtransf outest=py&i;
model fitted=py&i;
run;
%end;
%mend regit;
%regit;

data py1;
set py1;
power =-1;
keep _RMSE_ power;
run;
data py2;
set py2;
power =-1/2;
keep _RMSE_ power;
run;
data py3;
set py3;
power =-1/3;
keep _RMSE_ power;
run;
data py4;
set py4;
power =0;
keep _RMSE_ power;
run;
data py5;
set py5;
power =1/3;
keep _RMSE_ power;
run;
data py6;
set py6;
power =1/2;
keep _RMSE_ power;
run;
data py7;
set py7;
power =1;
keep _RMSE_ power;
run;

data rssdata;
set py1 py2 py3 py4 py5 py6 py7	;
rss = (_RMSE_**2)*(250-2);
run;

goptions reset = all;
symbol1 i=join c=black;
axis1 label=(h=2 f=times a=90 
     "Residual Sum of Squares") order=(0 to 
     50000 by 10000) value=(h=1 f=times);
axis2 label=(h=2 f=times "Lambda") order=(-1 to 1 by 0.5) 
     value=(h=1 f=times);
proc gplot data = rssdata;
 plot rss*power/vaxis=axis1 haxis=axis2 
     hminor=0 vminor=0;
run;



*Figure 3.30;
ods output boxcox=b details=d; 
   ods exclude boxcox; 
proc transreg details data = transf;
 model boxcox(y/ convenient 
    lambda=0.28 to 0.4 by .001)=identity(x);
 output out=trans;
run;

* Store values for reference lines; 
data _null_; 
 set d; 
 if description = 'CI Limit' 
     then call symput('vref',   formattedvalue); 
 if description = 'Lambda Used' 
     then call symput('lambda', formattedvalue); 
run; 

proc print data = b;
run;

data _null_; 
 set b end=eof; 
 where ci ne ' '; 
 if _n_ = 1 
    then call symput('href1', 
    compress(put(lambda, best12.))); 
 if ci  = '<' 
     then call symput('href2', 
     compress(put(lambda, best12.))); 
 if eof 
     then call symput('href3', 
     compress(put(lambda, best12.))); 
run;

*Plot 1:;
goptions reset = all;
axis2 order=(0.28 to 0.4 by 0.02) 
     label=(f=times h=2 "Lambda") value=(h=1 
     f=times);
axis1 order=(-20 to 30 by 10) label=(angle=90 
     f=times h=2 "log-Likelihood") value=(h=1
     f=times);
proc gplot data = b;
plot loglike * lambda / vref=&vref href=&href1 &href2 &href3 
     vminor=0 hminor=0 vaxis=axis1 haxis=axis2; 
 footnote height=1.5 font=times 
     "95% CI: &href1 - &href3, " 
     "Lambda = &lambda"; 
 symbol v=none i=spline c=black; 
run; 
 footnote; 


*Plot 2;
goptions reset = all;
axis2 order=(0.32 to 0.345 by 0.005) 
     label=(f=times h=2 "Lambda") value=(h=1 
     f=times);
axis1 order=(24.5 to 27 by .5) label=(angle=90 
     f=times h=2 "log-Likelihood") value=(h=1
     f=times);
proc gplot data = b;
plot loglike * lambda / vref=&vref href=&href2  
     vminor=0 hminor=0 vaxis=axis1 haxis=axis2; 
 symbol v=none i=spline c=black; 
run; 

data new;
 set transf;
 ty = y**(1/3);
 run;

proc reg data = new;
 model ty=x;
run;

%kerneld(dat=new,var=ty);
%boxplot(dsn=new,var=ty);
%qqplot(dsn=new, var=ty);
goptions reset = all;
 axis1 label=(font=times h=2 'x') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Y^(1/3)') value=(font=times h=1);
 symbol1 value = circle i=r;
proc gplot data = new;
 plot ty*x/  haxis=axis1 vaxis=axis2 vminor=0 hminor=0;
run;


data sal;
 infile 'data/salarygov.txt' firstobs=2;
 input id $ nw ne score maxsal;
run;

proc reg data = sal;
  model maxsal=score;
  output out=outreg student=stdres;
run;

data outreg;
set outreg;
sqtabsres = sqrt(abs(stdres));
run;
proc print data=outreg;
run;

*3.32;
goptions reset = all;
goptions reset = all;
 axis1 label=(font=times h=2 'Score') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Max Salary') value=(font=times h=1) 
     order=(1000 to 9000 by 1000);
 symbol1 value = circle i=r;
proc gplot data = outreg;
 plot maxsal*score/ haxis=axis1 vaxis=axis2 vminor=0 hminor=0;
run;

goptions reset = all;
goptions reset = all;
 axis1 label=(font=times h=2 'Score') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Standardized Residuals') value=(font=times h=1) 
     order=(-4 to 8 by 2);
 symbol1 value = circle i=r;
proc gplot data = outreg;
 plot stdres*score/ haxis=axis1 vaxis=axis2 vminor=0 hminor=0;
run;

goptions reset = all;
goptions reset = all;
 axis1 label=(font=times h=2 'Score') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Sqrt(Abs(Std Res))') value=(font=times h=1) 
     order=(0 to 3 by 1);
 symbol1 value = circle i=r;
proc gplot data = outreg;
plot sqtabsres*score/ haxis=axis1 vaxis=axis2 vminor=0    hminor=0;
run;


%macro mboxcox;
start geoMean(orMatrix,nuMatrix);
nuMatrix = log(orMatrix);
nuMatrix = nuMatrix[+,];
nuMatrix = exp(nuMatrix/nrow(orMatrix));
finish;

start h_bcmll(lambda) global(X,namesX);
GMX = J(1,ncol(X),0);
Y = X;
run geoMean(X,GMX);
i=1;
do while (i<=ncol(X));  
if lambda[1,i] = 0 then Y[,i]= GMX[1,i]*log(Y[,i]);
else Y[,i] = (GMX[1,i]##(1-lambda[1,i])) * (Y[,i]##lambda[1,i] - J(nrow(X),1,1)) / lambda[1,i];
i = i + 1;
end;
h1 = I(nrow(X)) - J(nrow(X),1,1) * J(nrow(X),1,1)`/nrow(X);
h2 = Y`*h1*Y/(nrow(X)-1);
h3 = det(h2);
if h3 <= 0 then print "Cannot estimate this transformation";
f = -nrow(X)*log(h3)/2; 
return(f);
finish h_bcmll;

start bcmll(lambda) global(X,namesX);
f = h_bcmll(lambda);
return(f);
finish bcmll;

start MultivariateBoxCox;
lambda = J(1,ncol(X),0);
optn =  j(1,11,.);
optn[1] = 1;
optn[2] = 2;
call nlpqn(rc,lambdares,"bcmll",lambda,optn);

call nlpfdd(crit, grad, hess, "bcmll",lambdares);
print grad;
variance = inv(-hess);
print variance;
print lambdares[c=namesX];

stderr = sqrt(vecdiag(variance)) ;
lrt0 = 2#(bcmll(lambdares)-bcmll(0#lambdares));
lrt1 = 2#(bcmll(lambdares)-bcmll(J(1,ncol(X),1)));
wald0= lambdares#(t(stderr##(-1)))	;
wald1= (lambdares-J(1,ncol(X),1))#(t(stderr##(-1)));
wald0 = t(wald0) ;
wald1 = t(wald1) ;
plrt0 = 1-CDF('CHISQUARE',lrt0,ncol(X));
plrt1 = 1-CDF('CHISQUARE',lrt1,ncol(X));

res = t(lambdares) || stderr || wald0 || wald1;

tcols = {'Power', 'Std.Error', 'Wald 0', 'Wald 1'};
resrow = t(namesX);
rescol= t(tcols);
print res[r=resrow c=rescol] ;

lrtert = J(2,3,0);
lrtert[1,1] = lrt0;
lrtert[2,1] = lrt1;
lrtert[1,2] = ncol(X);
lrtert[2,2] = ncol(X);
lrtert[1,3] = plrt0;
lrtert[2,3] = plrt1;

resrow = {'LRT all = 0', 'LRT all = 1'};
rescol = {'LRT','df','p-value'};
rescol = t(rescol);

print lrtert[r=resrow c=rescol];
finish;

run MultivariateBoxCox;


%mend mboxcox;

proc iml;
  use sal;
  read all ;
  namesX={"Salary" "Score"};
X  = maxsal || score ;
%mboxcox;
quit;


*3.33;
%kerneld(dat=sal,var=maxsal);
%boxplot(dsn=sal,var=maxsal);
%qqplot(dsn=sal, var=maxsal);
%kerneld(dat=sal,var=score);
%boxplot(dsn=sal,var=score);
%qqplot(dsn=sal, var=score);

*3.34;
data transf;
 set sal;
 logmax = log(maxsal);
 rtscore = sqrt(score);
run;
goptions reset = all;
symbol1 v=circle i=r c=black;
axis1  value=(f=times h=1)
     label=(h=2 f=times angle=90 "log(Max Salary)");
axis2  value=(f=times h=1) order=(5 to 35 by 10)
     label=(f=times h=2 "Sqrt(Score)");
proc gplot data = transf;
 plot logmax*rtscore/vaxis=axis1 haxis=axis2 vminor=0
     hminor=0;
run;

*3.35;
%kerneld(dat=transf,var=logmax);
%boxplot(dsn=transf,var=logmax);
%qqplot(dsn=transf, var=logmax);
%kerneld(dat=transf,var=rtscore);
%boxplot(dsn=transf,var=rtscore);
%qqplot(dsn=transf, var=rtscore);

*3.36;
proc reg data = transf;
 model logmax = rtscore;
 output out = outreg student=stdres;
run;

data plot2;
 set outreg;
 rtabsres = sqrt(abs(stdres));
run;

*Plot 1;
goptions reset = all;
symbol1 v=circle;
axis1 value=(f=times h=1) order=(-4 to 4 by 2)
     label=(h=2 f=times angle=90 "Standardized Residuals");
axis2 value=(f=times h=1) order=(5 to 35 by 10)
     label=(f=times h=2 "Sqrt(Score)");
proc gplot data = outreg;
 plot stdres*rtscore/vaxis=axis1 haxis=axis2 vminor=0
     hminor=0;
run;

*Plot 2;
goptions reset = all;
symbol1 v=circle interpol=r;
axis1 value=(f=times h=1) order=(0 to 2 by 0.5)
     label=(h=2 f=times angle=90 "Sqrt(Abs(Stdzd Res))");
axis2 value=(f=swiss h=1) order=(5 to 35 by 10) 
     label=(f=time h=2 "Sqrt(Score)");
proc gplot data = plot2;
 plot rtabsres*rtscore/vaxis=axis1 haxis=axis2 vminor=0
     hminor=0;
run;

proc iml;
  use sal;
  read all ;
  namesX={'Score'};
X  = score ;
%mboxcox;
quit;

proc reg data = transf noprint;
 model maxsal=rtscore;
 output out=regout p=fitted;
run;

data invresplot;
set regout;
y=maxsal	;
cbrty=y**(-.19);
ly=log(y) ;
run;


%macro regouts(dsn=,yvar=,predname=);
 proc reg data = invresplot noprint;
  model fitted = &yvar;
  output out= &dsn p=&predname;
 run;
%mend;

%regouts(dsn=data1,yvar=y, predname=lam1hat);
%regouts(dsn=data2,yvar=cbrty,predname=cbrtyhat);
%regouts(dsn=data3,yvar=ly,predname=lyhat);

proc sort data = invresplot;
 by y;
run;

%macro sortit;
%do i = 1 %to 3;
 proc sort data = data&i;
  by y;
 run;
%end;
%mend sortit;
%sortit;

data full;
 merge invresplot data1 data2 data3;
 by y;
run;


*3.37;
goptions reset = all;
symbol1 v=circle c=black;
symbol2 i=join c=black l=5 w=2;
symbol3 i=join c=black l=1 w=2;
symbol4 i=join c=black l=2 w=2;
axis1 label=(h=2 angle=90 f=times "yhat") 
     value=(h=1 f=times);
axis2 label=(h=2 f=times "y")
value=(h=1 f=times) offset=(2,0);
legend1 label=(f=times h=1.5 j=c 'Lambda' position=top) 
     position=(bottom right inside) across=1 frame 
     value=(f=times h=1.5 j=c 'Yhat' j=c '1' 
     j=c "-.19" j=c '0' j=c);
proc gplot data = full;
 plot fitted*y=1 lam1hat*y=2 cbrtyhat*y=3 lyhat*y=4/
     overlay vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0 legend=legend1;
run;

*3.38;
data transf;
 set transf;
 maxsaln25 =maxsal**(-.25);
run;
%kerneld(dat=transf,var=maxsaln25);
%boxplot(dsn=transf,var=maxsaln25);
%qqplot(dsn=transf, var=maxsaln25);


proc reg data = transf noprint;
 model maxsaln25=rtscore;
 output out=regout student=stdres;
run;

data regout;
set regout ;
rtabsres = sqrt(abs(stdres));
run;

goptions reset = all;
symbol1 v=circle i=r;
axis1  value=(f=times h=1)
     label=(h=2 f=times angle=90 "(Max Salary)^-.25");
axis2  value=(f=times h=1) order=(5 to 35 by 10)
     label=(f=times h=2 "Sqrt(Score)");
proc gplot data = regout;
 plot maxsaln25*rtscore/vaxis=axis1 haxis=axis2 vminor=0
     hminor=0;
run;
quit;

*Plot 2;

goptions reset = all;
symbol1 v=circle;
axis1  value=(f=times h=1)
     label=(h=2 f=times angle=90 "Standardized Residuals");
axis2  value=(f=times h=1) order=(5 to 35 by 10)
     label=(f=times h=2 "Sqrt(Score)");
proc gplot data = regout;
 plot stdres*rtscore/vaxis=axis1 haxis=axis2 vminor=0
     hminor=0;
run;


*Plot 3;
goptions reset = all;
symbol1 v=circle interpol=r;
axis1  value=(f=times h=1)
     label=(h=2 f=times angle=90 "Sqrt(Abs(Stdzd Res))");
axis2  value=(f=times h=1) order=(5 to 35 by 10)
     label=(f=times h=2 "Sqrt(Score)");
proc gplot data = regout;
 plot rtabsres*rtscore/vaxis=axis1 haxis=axis2 vminor=0
     hminor=0;
run;

