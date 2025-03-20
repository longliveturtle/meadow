data profsal; 
 infile 'data/profsalary.txt' firstobs=2 expandtabs;
 input case salary exper;
run;
quit;

*fig. 5.1;
goptions reset = all;
 axis1 label=(font=times h=2 'Years of Experience') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Salary') value=(font=times h=1) ;
 symbol1 value = circle;
proc gplot data = profsal;
 plot salary*exper/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;

proc reg data = profsal;
 model salary = exper;
 output out=outreg student=stdres;
run;
quit;

*fig. 5.2;
goptions reset = all;
 axis1 label=(font=times h=2 
     'Years of Experience') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Standardized Residuals') 
     value=(font=times h=1) ;
 symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*exper/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;

*fig. 5.3;
goptions reset = all;
 axis1 label=(font=times h=2 'Years of Experience') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Salary') value=(font=times h=1) ;
 symbol1 value = circle interpol=rq;
proc gplot data = profsal;
 plot salary*exper/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;

data quad;
 set profsal;
 expsq = exper**2;
run;
quit;

proc reg data = quad;
 model salary=exper expsq;
 output out=regout r=resids student=stdres cookd=cd 
     p=fitted h=levg;
run;
quit;

*fig. 5.4;
goptions reset = all;
 axis1 label=(font=times h=2 
     'Years of Experience') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Standardized Residuals') 
     value=(font=times h=1) ;
 symbol1 value = circle;
proc gplot data = regout;
 plot stdres*exper/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;

*fig. 5.5;
goptions reset = all;
 axis1 label=(font=times h=2 
     'Years of Experience') 
     value=(font=times h=1);
 axis2 label=(h=2 font=times angle=90 
     'Leverage') 
     value=(font=times h=1) ;
 symbol1 value = circle;
proc gplot data = regout;
 plot levg*exper/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0 vref=0.042;
run;
quit;

%macro plotlm(regout =,);
proc loess data = &regout;
 model resids=fitted/smooth=0.6667;
 ods output OutputStatistics=loessout;
run;
quit;

data fit;
 set &regout;
 set loessout;
run;
quit;

proc sort data = fit;
 by fitted;
run;
quit;

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
quit;

goptions reset = all htext=1.5;
 title1 height=2 font=times "Normal Q-Q";
 symbol1 value=circle color=black;

proc univariate data = &regout noprint;
 qqplot stdres/normal(mu=0 sigma=1 l=1 color=black) 
     font=times vminor=0 hminor=0 
     vaxislabel= "Standardized Residuals";
run;
quit;

data plot3;
 set &regout;
 sqrtres = sqrt(abs(stdres));
run;
quit;

proc loess data = plot3;
 model sqrtres=fitted/smooth=0.6667;
 ods output OutputStatistics=loessout1;
run;
quit;

data fit2;
 set plot3;
 set loessout1;
run;
quit;

proc sort data = fit2;
 by fitted;
run;
quit;

goptions reset = all;
symbol1 v=circle c=black;
symbol2 i=join c=black;
 title1 height=2 font=times  "Scale-Location";
 axis1 label = (font=times h=2 angle=90 
     'Sqrt(Abs(Res))')  
     value=(font=times h=1);
 axis2 label = (font=times h=2 'Fitted values') 
     value =(font=times h=1);
proc gplot data = fit2;
 plot /*points:*/ sqrtres*fitted=1 /*loess:*/ 
     Pred*fitted=2/ overlay hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;
quit;

proc sort data = &regout;
 by levg;
 run;
quit;

proc loess data = &regout;
 model stdres=levg/smooth=0.67777;
 ods output OutputStatistics=loessout2;
run;
quit;

data fit3;
 set &regout;
 set loessout2;
run;
quit;

proc sort data = fit3;
 by levg;
run;
quit;
goptions reset = all;
symbol1 v=circle c=black;
symbol2 i=join c=black;
 title1 height=2 font=times "Residuals vs Leverage";
 axis1 label = (h=2 font=times angle=90 
     "Standardized Residuals")
     value=(font=times h=1);
 axis2 label = (h=2 font=times 'Leverage') 
      value =(font=times h=1) ;
proc gplot data = fit3;
 plot /*points:*/ stdres*levg=1 /*loess:*/ Pred*levg=2/
     overlay hminor=0 vminor=0 vaxis=axis1 haxis=axis2 
     vref=0 href=0;
run;
quit;
%mend;


*fig. 5.6;
%plotlm(regout=regout);

*p. 129-130;
proc reg data = quad noprint;
 model salary = exper expsq;
 output out=predout p=fit lcl=lwr ucl=upr;
run;
quit;

data predout;
set predout;
keep exper fit lwr upr;
if exper in(10) then output;
run;
quit;

proc sort data=predout noduplicates;
  by exper;
run;
quit;

proc print data = predout noobs;
 var exper fit lwr upr;
run;
quit;

proc import datafile="data/nyc.csv" out=nyc replace;
 getnames=yes;
run;
quit;

*p. 138-139;
proc reg data = nyc;
 model price = food decor service east;
run;
quit;

*p. 139;
proc reg data = nyc;
 model price = food decor east;
run;
quit;

data travel;
 infile 'data/travel.txt' firstobs = 2 expandtabs;
 input amount age segment $ c;
run;
quit;

data tog; 
 set travel;
 if segment = 'A' then amta=amount; 
 if segment = 'C' then amtc=amount; 
run;
quit;

*fig. 5.7;
goptions reset = all;
 symbol1 font=times h=1 v='A' c=black;
 symbol2 font=times h=1 v='C' c=red;
 axis1 label = (h=2 font=times angle=90 
     "Amount Spent") 
     value=(font=times h=1);
 axis2 label = (h=2 font=times 'Age') 
      value =(font=times h=1);
proc gplot data = tog;
 plot amta*age=1 amtc*age=2/ 
     hminor=0 vminor=0 vaxis=axis1 haxis=axis2 
     overlay;
run;
quit;

data interact;
 set travel;
 inter = age*c;
run;
quit;

*p. 141;
proc reg data = interact;
 model amount = age c inter;
run;
quit;

*p. 143;
proc reg data = interact;
 model amount = age c;
run;
quit;

%macro ftest(dsn=, yvar=, fullm=, reducedm=,fulldf=,reddf=);
proc reg data = &dsn outest=est tableout;
 m1: model &yvar = &fullm/noprint;
 m2: model &yvar = &reducedm/noprint;
run;
quit;

data makef;
 set est;
 keep _RMSE_;
 if _TYPE_ ^= 'PARMS' then delete;
run;
quit;

proc transpose data = makef out=makef2(rename=(col1=rmse_full 
     col2=rmse_red));
run;
quit;

data make3;
 set makef2;
 keep rmse_full rmse_red;
run;
quit;

data ftest;
 set make3;
 rss_red=(rmse_red**2)*&reddf;
 rss_full=(rmse_full**2)*&fulldf;
 f=((rss_red - rss_full)/(&reddf-&fulldf))/(rss_full/&fulldf);
 pval=1-probf(f,(&reddf-&fulldf),&fulldf);
run;
quit;

proc print data = ftest;
run;
quit;
%mend ftest;

*p. 144;
%ftest(dsn=interact, yvar=amount, fullm= age c inter, 
     reducedm= age, fulldf=921, reddf=923);

data interac;
 set nyc;
 food_e = food*east;
 dec_e = decor*east;
 serv_e = service*east;
run;
quit;

*p. 145;
proc reg data = interac;
 model price = food decor service east food_e dec_e serv_e;
run;
quit;

*p. 146;
proc reg data = interac;
 model price = food decor east;
run;
quit;

*p. 146;
%ftest(dsn=interac, yvar=price, fullm= 
     food decor service east food_e dec_e serv_e, 
     reducedm= food decor east, fulldf=160, 
   reddf=164);

