data bridge;
 infile 'data/bridge.txt' expandtabs firstobs=2;
 input case time darea ccost dwgs length spans;
run;
quit;

data bridge;
 set bridge;
 log_Time = log(time);
 log_darea = log(darea);
 log_ccost = log(ccost);
 log_dwgs = log(dwgs);
 log_length = log(length);
 log_spans = log(spans);
run;
quit;

*p. 234;
proc reg data = bridge;
 model log_time = log_darea log_ccost log_dwgs 
     log_length log_spans;
run;
quit;


proc reg data = bridge outest=rsqadj edf;
 model log_time = log_darea log_ccost log_dwgs 
     log_length log_spans/ selection=adjrsq aic bic sse;
run;
quit;

data clean;
 set rsqadj;
 idarea = 0;
 if log_darea ^= . then idarea = 1;
 iccost = 0;
 if log_ccost ^= . then iccost = 1;
 idwgs = 0;
 if log_dwgs ^= . then idwgs=1;
 ilen = 0;
 if log_length ^= . then ilen=1;
 ispans=0;
 if log_spans ^= . then ispans=1;
 sst = _SSE_ / (1-_RSQ_);
 rsqadj = 1 - (_SSE_/(45-_IN_ -1))/(sst/(45-1));
 aicc = _AIC_ + (2*(_IN_+2)*(_IN_+3)/(45 - _IN_ - 1));
  bic = 45*log(_SSE_/45)+log(45)*(_IN_+1) ;
 drop _MODEL_ _TYPE_ _DEPVAR_ _RMSE_ intercept log_darea 
     log_ccost log_dwgs log_length log_spans log_time 
     _P_ _EDF_ _BIC_ sst _RSQ_;
run;
quit;

proc sort data = clean;
 by  _IN_ _SSE_;
run;
quit;

proc print data=clean;
run;
quit;

data table7p1;
 set clean;
 by _IN_ ;
 if (FIRST._IN_ EQ 1) then output;
run;
quit;

*table 7.1;
proc print data=table7p1 noobs;
 var _IN_ idwgs ispans iccost idarea ilen 
     rsqadj _AIC_ aicc bic;
run;
quit;


*fig. 7.1;
goptions reset = all;
symbol1 v=circle c=black;
 axis1 label = (h=2 font=times angle=90 
     "Adjusted R-squared") 
     value=(font=times h=1);
 axis2 label = (h=2 font=times 'Subset Size') 
      value =(font=times h=1) offset=(2,2)
      order=(1 to 5 by 1);
proc gplot data = table7p1;
 plot rsqadj*_IN_ /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;
quit;

*p. 235-236;
proc reg data=bridge;
model log_time=log_dwgs log_spans;
model log_time=log_dwgs log_spans log_ccost;
run;
quit;


data train;
 infile 'data/prostateTraining.txt' 
     expandtabs firstobs=2;
 input ocase lcavol lweight age lbph svi 
     lcp gleason pgg45 lpsa;
run;
quit;

*fig. 7.2;
ods graphics on;
proc corr data = train plots=matrix;
  var lpsa lcavol lweight age lbph ;
run;
quit;
ods graphics off;

ods graphics on;
proc corr data = train plots=matrix;
  with 	svi lcp gleason pgg45;
  var lpsa lcavol lweight age lbph ;
run;
quit;
ods graphics off;

ods graphics on;
proc corr data = train plots=matrix;
  with lpsa lcavol lweight age lbph;
var svi lcp gleason pgg45;  
run;
quit;
ods graphics off;

ods graphics on;
proc corr data = train plots=matrix;
  var 	svi lcp gleason pgg45;
run;
quit;
ods graphics off;

proc reg data = train noprint;
 model lpsa = lcavol lweight age lbph svi lcp 
     gleason pgg45;
output out=outreg r=resids p=fitted student=stdres cookd=cd h=levg;
run;
quit;

%macro stdresplot(xvar=,xlab=);
goptions reset = all;
axis2 label=(f=times h=2 "&xlab") value=(h=1 
     f=times);
axis1 label=(angle=90 f=times h=2 
     "Standardized Residuals") value=(h=1 f=times);
symbol v=circle; 
proc gplot data = outreg;
plot stdres*&xvar /vminor=0 hminor=0 vaxis=axis1
     haxis=axis2; 
run; 
quit;
%mend stdresplot;

*fig. 7.3;
%stdresplot(xvar=lcavol,xlab=lcavol);
%stdresplot(xvar=lweight,xlab=lweight);
%stdresplot(xvar=age,xlab=age);
%stdresplot(xvar=lbph,xlab=lbph);
%stdresplot(xvar=svi,xlab=svi);
%stdresplot(xvar=lcp,xlab=lcp);
%stdresplot(xvar=gleason,xlab=gleason);
%stdresplot(xvar=pgg45,xlab=pgg45);
%stdresplot(xvar=fitted,xlab=Fitted Values);

*fig. 7.4;
goptions reset = all;
 axis1 label=(h=2 f=times 'Fitted Values') 
     value=(f=times h=1)
     order=(0 to 5 by 1);
 axis2 label=(h=2 angle=90 f=times
     'lpsa') value=(f=times h=1);
symbol1 value = circle i=r;
proc gplot data = outreg;
 plot lpsa*fitted/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
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

*fig. 7.5;
%plotlm(regout=outreg);

*output p 242-244;
proc reg data = train;
 model lpsa = lcavol lweight age lbph svi lcp 
     gleason pgg45/vif;
output out=regout r=resids p=fitted student=stdres cookd=cd h=levg;
run;
quit;

%macro mmplot(dsn=, yvar=, x1=, allx=, xlab=, ylab=);
proc loess data = &dsn;
 model &yvar=&x1/smooth=0.6667;
 ods output OutputStatistics=loessout1;
run;
quit;

proc reg data = &dsn noprint;
 model &yvar = &allx;
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
     "&ylab") value=(f=times h=1);
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

*fig. 7.6;
%mmplot(dsn=train, yvar=lpsa, x1=lcavol, 
     allx= lcavol lweight age lbph svi lcp gleason 
     pgg45, xlab=lcavol, ylab=lpsa);
%mmplot(dsn=train, yvar=lpsa, x1=lweight, 
     allx= lcavol lweight age lbph svi lcp gleason 
     pgg45, xlab=lweight, ylab=lpsa);
%mmplot(dsn=train, yvar=lpsa, x1=age, 
     allx= lcavol lweight age lbph svi lcp gleason 
     pgg45, xlab=age, ylab=lpsa);
%mmplot(dsn=train, yvar=lpsa, x1=lbph, 
     allx= lcavol lweight age lbph svi lcp gleason 
     pgg45, xlab=lbph, ylab=lpsa);
%mmplot(dsn=train, yvar=lpsa, x1=lcp, 
     allx= lcavol lweight age lbph svi lcp gleason 
     pgg45, xlab=lcp, ylab=lpsa);
%mmplot(dsn=train, yvar=lpsa, x1=gleason, 
     allx= lcavol lweight age lbph svi lcp gleason 
     pgg45, xlab=gleason, ylab=lpsa);
%mmplot(dsn=train, yvar=lpsa, x1=pgg45, 
     allx= lcavol lweight age lbph svi lcp gleason 
     pgg45, xlab=pgg45, ylab=lpsa);
%mmplot(dsn=regout, yvar=lpsa, x1=fitted, 
     allx= lcavol lweight age lbph svi lcp gleason 
     pgg45, xlab=Fitted Values, ylab=lpsa);

	 
%macro addvarplot(dsn=, yvar=, ylab=, 
     var1=,othervar=);
 proc reg data = &dsn noprint;
  model &yvar = &othervar;
  output out=plot1a r=y1;
 run; 
quit;

 proc reg data = &dsn noprint;
  model &var1 = &othervar;
  output out = plot1b r=x1;
 run;
 quit;

 data plot1;
  merge plot1a plot1b;
 run;
 quit;

 goptions reset = all;
 axis1 label=(h=2 f=times "&var1|others") 
     value=(f=times h=1);
 axis2 label=(h=2 angle=90 f=times
      "&ylab|others") value=(f=times h=1);
 symbol1 value = circle i = r;
 title height=2 "Added-Variable Plot";
 proc gplot data = plot1;
  plot y1*x1/ haxis=axis1 vaxis=axis2 vminor=0
      hminor=0;
 run;
 quit;
%mend;

*fig. 7.7;
%addvarplot(dsn=train, yvar=lpsa, var1=lcavol, 
     ylab=lpsa,othervar=lweight age lbph svi lcp 
     gleason pgg45);
%addvarplot(dsn=train, yvar=lpsa, var1=lweight, 
ylab=lpsa,othervar=lcavol age lbph svi lcp 
     gleason pgg45);
%addvarplot(dsn=train, yvar=lpsa, var1=age, 
ylab=lpsa,othervar=lcavol lweight lbph svi lcp 
     gleason pgg45);
%addvarplot(dsn=train, yvar=lpsa, var1=lbph, 
ylab=lpsa,othervar=lcavol lweight age svi lcp 
     gleason pgg45);
%addvarplot(dsn=train, yvar=lpsa, var1=svi, 
ylab=lpsa,othervar=lcavol lweight age lbph lcp 
     gleason pgg45);
%addvarplot(dsn=train, yvar=lpsa, var1=lcp, 
ylab=lpsa,othervar=lcavol lweight age lbph svi 
     gleason pgg45);
%addvarplot(dsn=train, yvar=lpsa, var1=gleason, 
ylab=lpsa,othervar=lcavol lweight age lbph svi 
     lcp pgg45);
%addvarplot(dsn=train, yvar=lpsa, var1=pgg45, 
     ylab=lpsa,othervar=lcavol lweight age lbph svi 
     lcp gleason);


proc reg data = train outest=rsqadj edf;
 model lpsa = lcavol lweight age lbph svi lcp 
     gleason pgg45/ selection=adjrsq aic bic sse;
run;
quit;


data clean;
 set rsqadj;
 ilcavol = 0;
 if lcavol ^= . then ilcavol = 1;
 ilweight = 0;
 if lweight ^= . then ilweight = 1;
 iage = 0;
 if age ^= . then iage=1;
 ilbph = 0;
 if lbph ^= . then ilbph=1;
 isvi=0;
 if svi ^= . then isvi=1;
 ilcp=0;
 if lcp ^= . then ilcp=1;
 igleason=0;
 if gleason ^= . then igleason=1;
 ipgg45=0;
 if pgg45 ^= . then ipgg45=1;
 sst = _SSE_ / (1-_RSQ_);
 rsqadj = 1 - (_SSE_/(67-_IN_ -1))/(sst/(67-1));
 aicc = _AIC_ + (2*(_IN_+2)*(_IN_+3)/(67 - _IN_ - 1));
   bic = 67*log(_SSE_/67)+log(67)*(_IN_+1) ;
 drop _MODEL_ _TYPE_ _DEPVAR_ _RMSE_ intercept lcavol lweight age lbph svi lcp 
     gleason pgg45 _P_ _EDF_  _BIC_ sst _RSQ_;
run;
quit;

proc sort data = clean;
 by  _IN_ _SSE_;
run;
quit;


data table7p2;
 set clean;
 by _IN_ ;
 if (FIRST._IN_ EQ 1) then output;
run;
quit;

*table 7.2;
proc print data=table7p2 noobs;
 var _IN_ ilcavol ilweight iage ilbph isvi ilcp 
     igleason ipgg45 rsqadj _AIC_ aicc BIC;
run;
quit;


*fig. 7.8;
goptions reset = all;
symbol1 v=circle c=black;
 axis1 label = (h=2 font=times angle=90 
     "Adjusted R-squared") 
     value=(font=times h=1);
 axis2 label = (h=2 font=times 'Subset Size') 
      value =(font=times h=1);
proc gplot data = table7p2;
 plot rsqadj*_IN_ /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;
quit;

*p. 246;
proc reg data=train;
model lpsa=lcavol lweight;
model lpsa=lcavol lweight svi lbph;
model lpsa=lcavol lweight svi lbph pgg45 lcp age;
run;
quit;

data test;
 infile 'data/prostateTest.txt' 
     expandtabs firstobs=2;
 input ocase lcavol lweight age lbph svi 
     lcp gleason pgg45 lpsa;
run;

*p. 247;
proc reg data=test;
model lpsa=lcavol lweight;
model lpsa=lcavol lweight svi lbph;
model lpsa=lcavol lweight svi lbph pgg45 lcp age;
run;
quit;


data trainno45;
set train;
obs = _N_ ;
if obs EQ 45 then delete;
run;
quit;

proc print data=trainno45;
run;
quit;

proc reg data = trainno45 outest=rsqadj edf;
 model lpsa = lcavol lweight age lbph svi lcp 
     gleason pgg45/ selection=adjrsq aic bic sse;
run;
quit;


data clean;
 set rsqadj;
 ilcavol = 0;
 if lcavol ^= . then ilcavol = 1;
 ilweight = 0;
 if lweight ^= . then ilweight = 1;
 iage = 0;
 if age ^= . then iage=1;
 ilbph = 0;
 if lbph ^= . then ilbph=1;
 isvi=0;
 if svi ^= . then isvi=1;
 ilcp=0;
 if lcp ^= . then ilcp=1;
 igleason=0;
 if gleason ^= . then igleason=1;
 ipgg45=0;
 if pgg45 ^= . then ipgg45=1;
 sst = _SSE_ / (1-_RSQ_);
 rsqadj = 1 - (_SSE_/(66-_IN_ -1))/(sst/(66-1));
 aicc = _AIC_ + (2*(_IN_+2)*(_IN_+3)/(66 - _IN_ - 1));
   bic = 66*log(_SSE_/66)+log(66)*(_IN_+1) ;
 drop _MODEL_ _TYPE_ _DEPVAR_ _RMSE_ intercept lcavol lweight age lbph svi lcp 
     gleason pgg45 _P_ _EDF_  _BIC_ sst _RSQ_;
run;
quit;

proc sort data = clean;
 by  _IN_ _SSE_;
run;
quit;


data table7p2no45;
 set clean;
 by _IN_ ;
 if (FIRST._IN_ EQ 1) then output;
run;
quit;

proc print data=table7p2no45 noobs;
 var _IN_ ilcavol ilweight iage ilbph isvi ilcp 
     igleason ipgg45 rsqadj _AIC_ aicc BIC;
run;
quit;

*fig. 7.9;
goptions reset = all;
symbol1 v=circle c=black;
 axis1 label = (h=2 font=times angle=90 
     "Adjusted R-squared") 
     value=(font=times h=1);
 axis2 label = (h=2 font=times 'Subset Size') 
      value =(font=times h=1)
	        order=(1 to 5 by 1);
proc gplot data = table7p2;
 plot rsqadj*_IN_ /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;
quit;

goptions reset = all;
symbol1 v=circle c=black;
 axis1 label = (h=2 font=times angle=90 
     "Adjusted R-squared") 
     value=(font=times h=1);
 axis2 label = (h=2 font=times 'Subset Size') 
      value =(font=times h=1)
      order=(1 to 5 by 1);
proc gplot data = table7p2no45;
 plot rsqadj*_IN_ /hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;
quit;

data train1;
 set train;
 lpsa1 = lpsa;
 run;
 quit;

data test2;
 set test;
 lpsa2 = lpsa;
 run;
 quit;

data tog;
 set train1 test2;
run;
quit;

*fig. 7.10;
goptions reset = all;
symbol1 i=r v=plus c=red h=1;
symbol2 i=r v=triangle c=black h=1;
 axis1 label = (h=2 font=times angle=90 "lpsa") 
     value=(font=times h=1) order=(-1 to 6 by 1);
 axis2 label = (h=2 font=times 'lweight') 
     value =(font=times h=1);* order=(2.5 to 
     6.5 by 1);
 legend1 label=(f=times h=2 j=c 'Data Set' 
     position=top) position=(bottom right inside) 
     across=1 frame value=(f=times h=2 j=c 
     'Training' j=c 'Test');
proc gplot data = tog;
 plot lpsa1*lweight=1 lpsa2*lweight=2/
     overlay hminor=0 vminor=0 vaxis=axis1 
     haxis=axis2 legend=legend1;
run;
quit;

*fig. 7.11;
%addvarplot(dsn= test, yvar=lpsa, 
     ylab=lpsa, var1=lweight, 
     othervar=lcavol svi lbph pgg45 lcp age);

