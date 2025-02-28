proc import datafile="data/FieldGoals2003to2006.csv" 
    out=myfootb replace;
 getnames=yes;
run;
quit;

*p. 2;
proc corr data=myfootb pearson;
var FGt FGtM1; 
run;
quit;

*fig. 1.1;
goptions reset=all;
proc gplot data = myfootb;
 symbol v = circle;
 axis1 label = (h=2 font=times angle=90 
     justify=c 'Field Goal Percentage') 
     order=(60 to 100 by 10) value=(h=1 font=times) 
     offset=(0,2);
 axis2 label=(h=2 font=times 
     'Field Goal Percentage in Year t-1') order= 
     (65 to 100 by 5) value=(h=1 font=times) 
     offset=(2,2);
 title height=2 font=times 
     'Unadjusted Correlation = -0.139';
 plot FGt*FGtM1=1/hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2;
run;
quit;

*p. 3;
proc glm data=myfootb;
class Name;
model FGt=FGtm1 Name Name*FGtm1;
run;
quit;

*p. 3;
proc glm data=myfootb;
class Name;
model FGt = FGtM1 Name/solution;
run;
quit;

ods trace on;
ods RESULTS on;
proc glm data = myfootb;
 ods output ParameterEstimates=Parms; 
 class name;
 model FGt = FGtM1 name/solution;
run;
quit;
ods RESULTS off;
ods trace off;

data new;
 set parms;
 keep Estimate;
run;
quit;

proc transpose data = new out=params (rename=(_name_=
     name col1=int col2=slope col3-col20=x1-x18));
run;
quit;

%macro create;
 data create;
  set params;
  do x = 65 to 100 by 1;
    y0 = int + slope*x;
     %do i = 1 %to 18;
	   y&i = int + slope*x + x&i;
	 %end;
	output;
  end;
 run;
quit;
%mend create;
%create;

data final;
 set create myfootb;
run;
quit;

*fig. 1.2;
goptions reset = all;
%macro plotit;
 proc gplot data = final;
  symbol1 interpol=none value=circle color=black;
  %do i = 2 %to 20;
    symbol&i interpol=l line=&i color=black;
  %end;
 axis1 label =(h=1.5 font=times angle=90 'Field Goal 
     Percentage in Year t') order=(60 to 100 by 10) 
     value=(h=1.5 font=times);
 axis2 label = (h=1.5 font=times "Field Goal Percentage in Year 
     t-1") order= (65 to 100 by 5) value=(h=1.5 font=times) 
     offset=(2,2);
 title h=2 "Slope of each line = -0.504" font=times;
 plot FGt*FGtM1 y0*x y1*x y2*x y3*x y4*x y4*x y6*x y7*x 
     y8*x y9*x y10*x y11*x y12*x y13*x y14*x y15*x 
     y16*x y17*x y18*x /overlay haxis=axis2  
     hminor=0 vaxis=axis1 vminor=0;
run;
quit;
%mend;
%plotit;


proc import datafile="data/circulation.txt" 
    out=news replace;
run;
quit;

*fig. 1.3;
goptions reset = all;
proc gplot data = news;
 axis1 label = (f=times h = 2 angle=90 'Sunday Circulation') 
     order = (0 to 2000000 by 500000) value=(f=times h=1);
 axis2 label =(f=times h=2 'Weekday Circulation') 
     order=(100000 to 1300000 by 400000) value=(f=times h=1);
 symbol1 v = circle color = black;
 symbol2 v = triangle color = red;
 legend1 label = (f=times h=1 justify=c 'Tabloid Dummy Variable' position=top justify=center) 
     position = (top left inside) across=1 frame 
     value = (f=times h=1 '0' justify = c '1' justify = c) ;
plot sunday*weekday = tabloid/ legend=legend1 
     haxis = axis2 hminor=0 vaxis = axis1 vminor=0;
run;
quit;

data logs;
 set news;
 lsun = log(sunday);
 lweek = log(weekday);
run;
quit;

*fig. 1.4;
goptions reset = all;
proc gplot data = logs;
 axis1 label = (f=times h=2 angle=90 'log(Sunday Circulation)') 
    order = (11.5 to 14.5 by 0.5) value=(f=times h=1);
 axis2 label = (f=times h=2 'log(Weekly Circulation)') 
    order=(11.5 to 14 by 0.5) value = (f=times h=1);
 symbol1 v = circle color = black;
 symbol2 v = triangle color = red;
 legend1 label = (h=1 'Tabloid Dummy Variable' 
     position = top justify = center) position = (top 
     left inside) across=1 frame value = (f=times h=1
     '0' justify = c '1' justify = c);
  plot lsun*lweek = tabloid/ legend=legend1 haxis =axis2 
     hminor=0 vaxis = axis1 vminor=0;
run;
quit;

proc import datafile="data/nyc.csv" out=nyc replace;
 getnames=yes; 
run;
quit;

*fig. 1.5;
ods graphics on;
proc corr data = nyc plots=matrix;
  var price food decor service;
run;
quit;
ods graphics off;

*fig. 1.6;
%macro boxplot(var =,data =,yvar =);
goptions reset = all;
proc sort data = &data;
 by &var;
run;
quit;
proc gplot data = &data;
 axis1 label=(f=times h=2 "&var") minor=none 
    order=(-1 0.2 0.8 2) value=(f=times h=1 
    t=1 ' ' t=2 '0' t=3 '1' t=4 ' ');
 axis2 label=(f=times h=2 angle=90 "&yvar") 
    value=(f=times h=1);
 symbol1 value = circle interpol=boxt bwidth=36;
 plot &yvar*&var/ haxis=axis1 vaxis=axis2 vminor=0;
run;
quit;
%mend;
%boxplot(var=east,data=nyc,yvar=price);


proc import datafile="data/Bordeaux.csv" out=wine 
     replace;
 getnames=yes;
run;
quit;

*fig. 1.7;
ods graphics on;
proc corr data = wine plots=matrix;
  var price ParkerPoints CoatesPoints;
run;
quit;
ods graphics off;

*fig. 1.8;
%boxplot(var = P95andAbove,data=wine,yvar=price);
%boxplot(var = FirstGrowth,data=wine,yvar=price);
%boxplot(var = CultWine,data=wine,yvar=price);
%boxplot(var = Pomerol,data=wine,yvar=price);
%boxplot(var = VintageSuperstar,data=wine,yvar=price);

data logs;
 set wine;
 log_price = log(price);
 log_ParkerPoints= log(ParkerPoints);
 log_CoatesPoints = log(CoatesPoints);
run;
quit;

*fig. 1.9;
ods graphics on;
proc corr data = logs plots=matrix;
  var log_price log_ParkerPoints log_CoatesPoints;
run;
quit;
ods graphics off;

*fig. 1.10;
%boxplot(var = P95andAbove,data=logs,yvar=log_price);
%boxplot(var = FirstGrowth,data=logs,yvar=log_price);
%boxplot(var = CultWine,data=logs,yvar=log_price);
%boxplot(var = Pomerol,data=logs,yvar=log_price);
%boxplot(var = VintageSuperstar,data=logs,yvar=log_price);




