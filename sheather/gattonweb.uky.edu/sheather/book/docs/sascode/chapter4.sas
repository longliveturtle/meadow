data inputit;
 infile "data/cleaningwtd.txt" firstobs=2 expandtabs;
 input case crews rooms stdev;
run;
quit;

data wtd;
 set inputit;
 wt = 1/(stdev**2);
run;
quit;

*p. 117;
proc reg data = wtd;
 model rooms = crews;
 weight wt;
 output out=outreg p=fit lcl=lwr ucl=upr;
run;
quit;

data outreg;
set outreg;
keep crews fit lwr upr;
if crews in(4,16) then output;
run;
quit;

proc sort data=outreg noduplicates;
  by crews ;
run;
quit;

proc print data = outreg noobs;
 var crews fit lwr upr;
run;
quit;

data weightit;
 set wtd;
 ynew = sqrt(wt)*rooms;
 x1new = sqrt(wt);
 x2new = sqrt(wt)*crews;
run;
quit;

*p. 120;
proc reg data = weightit;
 model ynew = x1new x2new/noint;
 output out=outreg p=fit lcl=lwr ucl=upr;
run;
quit;

data outreg;
set outreg;
keep crews fit lwr upr;
if crews in(4,16) then output;
run;
quit;

proc sort data=outreg noduplicates;
  by crews ;
run;
quit;

proc print data = outreg noobs;
 var crews fit lwr upr;
run;
quit;
