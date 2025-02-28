data food1;
  infile 'data/confood2.txt' firstobs=2 expandtabs;
  input week sales price promotion saleslag1 $;
run;

data food;
 set food1;
 if saleslag1 = 'NA' then saleslag = .;
 else saleslag = saleslag1;
 drop saleslag1;
 lsales = log(sales);
 lprice = log(price);
 lsaleslag = log(saleslag);
run;

*Fig. 9.1;
data no;
 set food;
 if promotion=0;
 lsalesno = lsales;
run;
data yes;
 set food;
 if promotion=1;
 lsalesyes = lsales;
run;
data plotit;
 set no yes;
run;

goptions reset = all;
axis1 label=(h=2 f=times angle=90 'log(Sales)') 
    value=(h=1 f=times) offset=(1,1);
axis2 label=(h=2 f=times 'log(Price)') v=(h=1
    f=times)  offset=(1,1);
symbol1 v=triangle h=1 c=black;
symbol2 v=plus h=1 c=black;
legend1 across=1 frame 
    position=(top right inside) label=
 (h=2 f=times position=top
    'Promotion') value=(h=2 
    f=times 'No' h= f=times 'Yes');
proc gplot data = plotit;
 plot lsalesno*lprice=1 lsalesyes*lprice=2/ 
    vaxis=axis1 haxis=axis2 overlay
    vminor=0 hminor=0 legend=legend1;
run;

*Fig. 9.2;
proc sort data = plotit;
 by week;
run;

goptions reset = all;
axis1 label=(h=2 f=times angle=90 'log(Sales)') 
    value=(h=1 f=times) offset=(1,1);
axis2 label=(h=2 f=times 'Week, t') v=(h=1
    f=times)  offset=(1,1);
symbol1 i=join c=black l=1;
symbol2 v=triangle h=1 c=black;
symbol3 v=plus h=1 c=black;
legend1 across=1 frame 
    position=(bottom right inside) label=
 (h=2 f=times position=top
    'Promotion') value=(h=0 f=swiss ' ' h=2 
    f=times 'No' h=2 f=times 'Yes');
proc gplot data = plotit;
 plot lsales*week=1 lsalesno*week=2 lsalesyes*week=3/ 
    vaxis=axis1 haxis=axis2 overlay
    vminor=0 hminor=0 legend=legend1;
run;

*Fig. 9.3;
goptions reset = all;
axis1 label=(h=2 f=times angle=90 'log(Sales t)') 
    value=(h=1 f=times) offset=(2,1);
axis2 label=(h=2 f=times 'log(Sales t-1)') v=(h=1
    f=times)  offset=(2,1);
symbol1 v=circle c=black h=1;
proc gplot data = food;
 plot lsales*lsaleslag=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;


*Fig. 9.4;
ods html;
ods graphics on;
proc arima data = food;
 estimate plot;
 identify var=lsales;
run;
ods graphics off;
ods html close;

*Fig. 9.5;
proc reg data = food;
 model lsales = lprice week promotion;
 output out = resids_m2 student=stdres p=fitted;
run;

*Plot 1;
goptions reset = all;
axis1 value=(h=1 f=times) label=(h=2 f=times angle=90 
     'Standardized Residuals');
axis2 label=(h=2 f=times 'log(Price t)') v=(h=1
     f=times) order=(-.6 to -.1 by .1);
symbol1 v=circle c=black h=1;
proc gplot data = resids_m2;
 plot stdres*lprice=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 2;
goptions reset = all;
axis1 value=(h=1 f=times) label=(h=2 f=times angle=90 
     'Standardized Residuals');
axis2 label=(h=2 f=times 'log(Price t)') v=(h=1
     f=times) order=(-2 to 54 by 8);
symbol1 v=circle c=black h=1 i=join l=1;
proc gplot data = resids_m2;
 plot stdres*week=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 3;
goptions reset = all;
axis1 value=(h=1 f=times) label=(h=2 f=times angle=90 
     'Standardized Residuals');
axis2 label=(h=2 f=times 'log(Price t)') v=(h=1
     f=times) order=(0 to 1 by 1) offset=(15,15);
symbol1 v=circle c=black h=1;
proc gplot data = resids_m2;
 plot stdres*promotion=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 4;
goptions reset = all;
axis1 value=(h=1 f=times) label=(h=2 f=times angle=90 
     'Standardized Residuals');
axis2 label=(h=2 f=times 'log(Price t)') v=(h=1
     f=times) order=(5.5 to 8.5 by 1);
symbol1 v=circle c=black h=1;
proc gplot data = resids_m2;
 plot stdres*fitted=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;
*End Fig. 9.5;

*Fig. 9.6;
ods html;
ods graphics on;
title 'Series Standardized Residuals';
proc arima data = resids_m2;
 estimate plot;
 identify var=stdres;
run;
ods graphics off;
ods html close;


*output p. 313;
data food; set food; subjectnum=1; run;
proc mixed data=food method=ml;
  class subjectnum;
  model lsales=lprice promotion week/solution 
     residual;
  repeated /type=ar(1) subject=subjectnum; 
run;

proc mixed data=food method=ml;
  class subjectnum;
  model lsales=lprice promotion week/solution 
     residual outp=fig9p7;
  repeated /type=ar(1) subject=subjectnum; 
run;

*Fig. 9.7;
ods html;
ods graphics on;
title 'Series Standardized Residuals';
proc arima data = fig9p7;
 estimate plot;
 identify var=resid;
run;
ods graphics off;
ods html close;


*Output p. 318;
proc sort data = food;
 by week;
run;
%macro makecorr;
proc iml;
 use food;
 read all var{lprice promotion week} into x;
 read all var{week} into weeks;
 read all var{lsales} into y;
 rho = 0.5503593;
 int = J(nrow(x), 1, 1);
 design = int||x;
 rowsig = weeks
    %do i = 2 %to 52;
	   || weeks
	%end;
	;
 rowwks = t(weeks);
 colsig = rowwks
    %do i = 2 %to 52;
	  // rowwks
	%end;
	;
 sigma = rho##abs(rowsig-colsig);
 s = root(sigma)`; 
 max_diff=max(abs(s*s`-sigma));
 old_s=root(sigma);
 old_max_diff=max(abs(old_s*old_s`-sigma));
 ystar = solve(s,y);
 xstar = solve(s,design);
 data = ystar||xstar;
 data2 = data||x;
 names2={ystar intstar lprstar promstar wkstar 
     lprice promotion week};
 create foodstar from data2 [colname=names2];
 append from data2;
quit;
%mend;

%makecorr;

proc reg data = foodstar;
 model ystar = intstar lprstar promstar wkstar/noint;
 output out=outres student=stdres p=fitted r=resids 
     h=levg;
run;



*Fig. 9.8;
*Plot 1;
goptions reset = all;
axis1 value=(h=2 f=times) label=(h=2 f=times
     angle=90 'log(Sales)*');
axis2 label=(h=2 f=times 'Intercept*') v=(h=2
     f=times) order=(0.5 to 1.1 by 0.2);
symbol1 v=circle c=black h=2;
proc gplot data = foodstar;
 plot ystar*intstar=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 2;
goptions reset = all;
axis1 value=(h=2 f=times) label=(h=2 f=times
     angle=90 'log(Sales)*') offset=(1,0);
axis2 label=(h=2 f=times 'log(Price)*') v=(h=2
     f=times) order=(-0.6 to 0 by 0.2);
symbol1 v=circle c=black h=2;
proc gplot data = foodstar;
 plot ystar*lprstar=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 3;
goptions reset = all;
axis1 value=(h=2 f=times) label=(h=2 f=times 
     angle=90 'log(Sales)*') offset=(1,0);
axis2 label=(h=2 f=times 'Promotion*') v=(h=2
     f=times) order=(-1 to 1.5 by 0.5);
symbol1 v=circle c=black h=2;
proc gplot data = foodstar;
 plot ystar*promstar=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 4;
goptions reset = all;
axis1 value=(h=2 f=times) label=(h=2 f=times 
     angle=90 'log(Sales)*') offset=(1,0);
axis2 label=(h=2 f=times 'Week*') v=(h=2
     f=times) order=(0 to 30 by 10);
symbol1 v=circle c=black h=2;
proc gplot data = foodstar;
 plot ystar*wkstar=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;


*Fig. 9.9;
ods html;
ods graphics on;
title 'Series Standardized Residuals';
proc arima data = outres;
 estimate plot;
 identify var=stdres;
run;
ods graphics off;
ods html close;

*Fig. 9.10;
*Plot 1;
goptions reset = all;
axis1 value=(h=2 f=times) label=(h=2 f=times
     angle=90 'Standardized LS Residuals') 
     offset=(0,1);
axis2 label=(h=2 f=times 'log(Price)') v=(h=2
     f=times) order=(-.6 to -.1 by 0.1);
symbol1 v=circle c=black h=2;
proc gplot data = outres;
 plot stdres*lprice=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 2;
proc sort data = outres;
 by week;
run;
goptions reset = all;
axis1 value=(h=2 f=times) label=(h=2 f=times 
     angle=90 'Standardized LS Residuals');
axis2 label=(h=2 f=times 'Week') v=(h=2
     f=times) order=(0 to 52 by 13) offset=(1,2);
symbol1 v=circle c=black h=2 i=join;
proc gplot data = outres;
 plot stdres*week=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 3;
goptions reset = all;
axis1 value=(h=2 f=times) label=(h=2 f=times 
     angle=90 'Standardized LS Residuals');
axis2 label=(h=2 f=times 'Promotion') v=(h=2
     f=times) order=(0 to 1 by 1) offset=(15, 15);
symbol1 v=circle c=black h=2;
proc gplot data = outres;
 plot stdres*promotion=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;

*Plot 4;
goptions reset = all;
axis1 value=(h=2 f=times) label=(h=2 f=times 
     angle=90 'Standardized LS Residuals') 
     offset=(0,1);
axis2 label=(h=2 f=times 'Fitted Values') v=(h=2
     f=times) order=(2 to 7 by 1);
symbol1 v=circle c=black h=2;
proc gplot data = outres;
 plot stdres*fitted=1/vaxis=axis1 haxis=axis2 
    vminor=0 hminor=0;
run;
quit;
*End Fig. 9.10;



* output out=outres student=stdres p=fitted r=resids 
     h=levg;

*Fig 9.11;
%macro plotlm(regout =,);
proc loess data = &regout;
 model resids=fitted/smooth=0.6667;
 ods output OutputStatistics=loessout;
run;

data fit;
 set &regout;
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

%plotlm(regout=outres);
