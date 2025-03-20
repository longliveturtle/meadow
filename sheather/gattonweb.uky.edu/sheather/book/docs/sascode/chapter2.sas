proc import datafile="data/production.txt" 
    out=prod replace;
run;
quit;

*fig. 2.1;
goptions reset = all;
proc gplot data = prod;
 symbol1 v=circle;
 axis1 label = (font=times h=2 angle=90 'Run Time') 
     order=(140 to 260 by 20) value=(font=times h=1);
 axis2 label = (font=times h=2 'Run Size') 
     order=(50 to 350 by 50) value =(font=times h=1);
 plot runtime*runsize/hminor=0 vminor=0 vaxis=axis1 
     haxis=axis2 ;
run;
quit;

*p. 19;
proc reg data = prod;
 model runtime=runsize;
run;
quit;

*fig. 2.3;
goptions reset = all;
proc gplot data = prod;
 symbol1 v=circle interpol=r;
 axis1 label = (h=2 font=times angle=90 'Run Time') 
     order=(140 to 260 by 20) value=(font=times h=1);
 axis2 label = (font=times h=2 'Run Size') 
     order=(50 to 350 by 50) value =(font=times h=1);
 plot runtime*runsize/hminor=0 vminor=0 vaxis=axis1 
     haxis=axis2 ;
run;
quit;

data t_val;
  tpvalue =2*cdf('t',-6.98,18,);
run;
quit;

*p. 22;
proc print data=t_val;
run;
quit;

*p. 24;
proc reg data = prod;
 model runtime = runsize/ clb;
run;
quit;

*p. 27;
data cipi;
 input runsize @@;
 cards;
 50 100 150 200 250 300 350
 ;
run;
quit;

data new;
 set prod cipi;
run;
quit;

proc reg data = new;
 model runtime = runsize/clm cli noprint;
 output out = messy p=fit lclm=lwr uclm=upr lcl=lwr_pi ucl=upr_pi;
run;
quit;

data first;
 set messy;
 if _N_ > 20;
 drop case runtime runsize lwr_pi upr_pi;
run;
quit;

proc print data = first;
run;
quit;

data second;
 set messy;
 if _N_ > 20;
 drop case runtime runsize lwr upr;
run;
quit;

proc print data = second;
run;
quit;

*p. 30;
proc reg data = prod;
 model runtime = runsize;
run;
quit;

data change;
 infile 'data/changeover_times.txt' firstobs=2 
     expandtabs;
 input method $ y x;
run;
quit;

*p. 31;
proc reg data = change;
 model y = x;
run;
quit;

*p. 32;
data t_val;
  tpvalue =2*cdf('t',-2.254,118);
run;
quit;

proc print data=t_val;
run;
quit;

*fig. 2.5;
goptions reset = all;
proc gplot data = change;
 symbol v=circle interpol=r;
 axis1 label=(h=2 angle=90 font=times 
     'Change Over Time') order=(5 to 40 by 5) 
     value=(font=times h=1) offset=(2,2);
 axis2 label=(h=2 font=times 
     'Dummy variable, New') 
     offset=(10,10) order=(0 to 1 by 0.2) 
     value=(font=times h=1);
 plot y*x/hminor=0 vminor=0 vaxis=axis1 haxis=axis2;
run;
quit;


goptions reset = all;
proc boxplot data = change;
  symbol1 value=none;
  label x= 'Dummy variable, New';
  label y= 'Change Over Time';
  plot y*x/ boxwidth=35 cboxes=black 
     boxstyle=schematic idsymbol=circle hoffset=4
     height=3 font=times;
run;
quit;

goptions reset = all;
proc gplot data = change;
 axis1 label=(font=times h=2 'Method') 
     order=(-1 0.2 0.8 2) value=(font=times h=2 t=1 
     ' ' t=2 'Existing' t=3 'New' t=4 ' ');
 axis2 label=(h=2 font=times angle=90 
     'Change Over Time') 
     order=(5 to 40 by 5) value=(font=times h=1);
 symbol1 value = circle interpol=boxt bwidth=36;
 plot y*x/  haxis=axis1 vaxis=axis2 vminor=0;
run;
quit;


