proc import datafile="data/Orthodont.txt" 
    out=orthodont replace;
 getnames=yes;
run;
quit;

data fig101;
set orthodont;
if subject='F01' then af1 = age;
if subject='F01' then df1 = distance;
if subject='F02' then af2 = age;
if subject='F02' then df2 = distance;
if subject='F03' then af3 = age;
if subject='F03' then df3 = distance;
if subject='F04' then af4 = age;
if subject='F04' then df4 = distance;
if subject='F05' then af5 = age;
if subject='F05' then df5 = distance;
if subject='F06' then af6 = age;
if subject='F06' then df6 = distance;
if subject='F07' then af7 = age;
if subject='F07' then df7 = distance;
if subject='F08' then af8 = age;
if subject='F08' then df8 = distance;
if subject='F09' then af9 = age;
if subject='F09' then df9 = distance;
if subject='F10' then af10 = age;
if subject='F10' then df10 = distance;
if subject='F11' then af11 = age;
if subject='F11' then df11 = distance;
run;
quit;

*fig. 10.1;
goptions reset = all;
symbol1 v=circle c=blue;
symbol2 i=join c=blue;
axis1 label =(h=2 font=times angle=90 "Distance") 
value=(font=times h=1)
	  order=(16 to 30 by 2);
axis2 label =(h=2 font=times "Age") value =(font=times h=1) offset=(2,2)
      order=(8 to 14 by 2);
axis3 major=none minor=none label=NONE value=NONE order=(16 to 30 by 2);
axis4 major=none minor=none label=NONE value=NONE;
	  goptions vsize=6;
	  goptions hsize=2;

proc sort data=fig101;
by age;
run;
quit;


proc gplot data = fig101; 
title 'F02';
plot df2*af2=2 df2*af2=1/overlay vaxis=axis1 haxis=axis4;
run;
quit;

proc gplot data = fig101;
title 'F08';
plot df8*af8=2 df8*af8=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig101;
 title 'F03';
 plot df3*af3=2 df3*af3=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig101;
 title 'F04';
 plot df4*af4=2 df4*af4=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig101;
 title 'F11';
 plot df11*af11=2 df11*af11=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig101;
 title 'F10';
 plot df10*af10=2 df10*af10=1/overlay vaxis=axis1 haxis=axis2;
run;
quit;
proc gplot data = fig101;
title 'F09';
plot df9*af9=2 df9*af9=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig101;
 title 'F06';
plot df6*af6=2 df6*af6=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig101;
  title 'F01';
 plot df1*af1=2 df1*af1=1/overlay vaxis=axis3 haxis=axis2;
 run;
 quit;
proc gplot data = fig101;
   title 'F05';
 plot df5*af5=2 df5*af5=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig101;
   title 'F07';
 plot df7*af7=2 df7*af7=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;

proc import datafile="data/orthdist.txt" 
    out=orthodontwide replace;
 getnames=yes;
run;
quit;

data orthodontwidef;
set orthodontwide;
if substr(subject,1,1)='F' then output;
run;
quit;

data orthodontwidef;
set orthodontwidef;
DistFAge8=Age8;
DistFAge10=Age10;
DistFAge12=Age12;
DistFAge14=Age14;
run;
quit;

*p. 334, fig. 10.2;
ods graphics on;
proc corr data=orthodontwidef plots=matrix;
var distfage14 distfage12 distfage10 distfage8;
run;
quit;
ods graphics off;

data female;
set orthodont;
if substr(subject,1,1)='F' then output;
run;
quit;

*p. 337; 
proc mixed data=female;
class subject;
model distance=age/solution intercept outp=frF;
random subject;
run;
quit;

*p. 337;
data res;
subject = sqrt(4.2786);
residual = sqrt(.6085);
run;
quit;
proc print data=res;
run;
quit;

*table 10.2;
data femalerandint;
set frF;
randint = pred-0.479545*age;
run;
quit;

proc means data=femalerandint noprint;
var randint;
by subject;
output out =randintest MEAN=;
run;
quit; 

ods trace on;
ods RESULTS on;
proc glm data=female;
 ods output ParameterEstimates=fixintest; 
class subject;
model distance = age subject/solution;
run; 
quit;
ods RESULTS off;
ods trace off;

data fixintest;
set fixintest;
estimate = estimate+21.1;
if _N_ < 3 then delete;
run;
quit;

proc print data=fixintest;
run;quit;

data fixintest;
set fixintest;
subject = substr(Parameter,11,3);
fixint=estimate;
run;
quit;

proc sort data=fixintest;
by subject;
run;
quit;

proc sort data=randintest;
by subject;
run;
quit;

data intest;
merge randintest fixintest;
by subject;
run;
quit;

data intest;
set intest;
if subject= 'F11' then order =1;
if subject= 'F04' then order =2;
if subject= 'F03' then order =3;
if subject= 'F08' then order =4;
if subject= 'F07' then order =5;
if subject= 'F02' then order =6;
if subject= 'F05' then order =7;
if subject= 'F01' then order =8;
if subject= 'F09' then order =9;
if subject= 'F06' then order =10;
if subject= 'F10' then order =11;
diff = randint-fixint;
run;
quit;

proc sort data=intest;
by order;
run;
quit;

proc print data=intest;
var subject randint fixint diff;
run;
quit;


*fig 10.3;
data fig103;
set frF;
if subject='F01' then ff1 = pred;
if subject='F01' then df1 = distance;
if subject='F02' then ff2 = pred;
if subject='F02' then df2 = distance;
if subject='F03' then ff3 = pred;
if subject='F03' then df3 = distance;
if subject='F04' then ff4 = pred;
if subject='F04' then df4 = distance;
if subject='F05' then ff5 = pred;
if subject='F05' then df5 = distance;
if subject='F06' then ff6 = pred;
if subject='F06' then df6 = distance;
if subject='F07' then ff7 = pred;
if subject='F07' then df7 = distance;
if subject='F08' then ff8 = pred;
if subject='F08' then df8 = distance;
if subject='F09' then ff9 = pred;
if subject='F09' then df9 = distance;
if subject='F10' then ff10 = pred;
if subject='F10' then df10 = distance;
if subject='F11' then ff11 = pred;
if subject='F11' then df11 = distance;
run;
quit;


data anno;
retain xsys ysys '1';
function = 'move'; x=0; y=0; output;
function = 'draw'; x=100; y=100; output;
run; 
quit;

goptions reset = all;
symbol1 v=circle c=blue;
axis1 label =(h=2 font=times angle=90 "Distance") 
value=(font=times h=1)
	  order=(16 to 28 by 2);
axis2 label =(h=2 font=times "Fitted") value =(font=times h=1) offset=(2,2)
      order=(16 to 28 by 2);
axis3 major=none minor=none label=NONE value=NONE order=(16 to 28 by 2);
axis4 major=none minor=none label=NONE value=NONE  order=(16 to 28 by 2);
	  goptions vsize=3;
	  goptions hsize=3;

proc gplot data = fig103 annotate=anno; 
title 'F04';
plot df4*ff4=1 /vaxis=axis1 haxis=axis4;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F11';
plot df11*ff11=1 /vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F02';
plot df2*ff2=1 /vaxis=axis1 haxis=axis4;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F08';
plot df8*ff8=1 /vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F03';
plot df3*ff3=1 /vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F01';
plot df1*ff1=1 /vaxis=axis1 haxis=axis4;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F05';
plot df5*ff5=1 /vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F07';
plot df7*ff7=1 /vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F10';
plot df10*ff10=1 /vaxis=axis1 haxis=axis2;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F09';
plot df9*ff9=1 /vaxis=axis3 haxis=axis2;
run;
quit;

proc gplot data = fig103 annotate=anno; 
title 'F06';
plot df6*ff6=1 /vaxis=axis3 haxis=axis2;
run;
quit;

data male;
set orthodont;
if substr(subject,1,1)='M' then output;
run;
quit;


*fig. 10.4;
data fig104;
set male;
if subject='M01' then am1 = age;
if subject='M01' then dm1 = distance;
if subject='M02' then am2 = age;
if subject='M02' then dm2 = distance;
if subject='M03' then am3 = age;
if subject='M03' then dm3 = distance;
if subject='M04' then am4 = age;
if subject='M04' then dm4 = distance;
if subject='M05' then am5 = age;
if subject='M05' then dm5 = distance;
if subject='M06' then am6 = age;
if subject='M06' then dm6 = distance;
if subject='M07' then am7 = age;
if subject='M07' then dm7 = distance;
if subject='M08' then am8 = age;
if subject='M08' then dm8 = distance;
if subject='M09' then am9 = age;
if subject='M09' then dm9 = distance;
if subject='M10' then am10 = age;
if subject='M10' then dm10 = distance;
if subject='M11' then am11 = age;
if subject='M11' then dm11 = distance;
if subject='M12' then am12 = age;
if subject='M12' then dm12 = distance;
if subject='M13' then am13 = age;
if subject='M13' then dm13 = distance;
if subject='M14' then am14 = age;
if subject='M14' then dm14 = distance;
if subject='M15' then am15 = age;
if subject='M15' then dm15 = distance;
if subject='M16' then am16 = age;
if subject='M16' then dm16 = distance;
run;
quit;


goptions reset = all;
symbol1 v=circle c=blue;
symbol2 i=join c=blue;
axis1 label =(h=2 font=times angle=90 "Distance") 
value=(font=times h=1)
	  order=(16 to 32 by 2);
axis2 label =(h=2 font=times "Age") value =(font=times h=1)
      order=(8 to 14 by 2);
axis3 major=none minor=none label=NONE value=NONE order=(16 to 32 by 2);
axis4 major=none minor=none label=NONE value=NONE order=(8 to 14 by 2);
	  goptions vsize=6;
	  goptions hsize=2;

proc sort data=fig104;
by age;
run;
quit;


proc gplot data = fig104; 
title 'M13';
plot dm13*am13=2 dm13*am13=1/overlay vaxis=axis1 haxis=axis4;
run;
quit;

proc gplot data = fig104; 
title 'M14';
plot dm14*am14=2 dm14*am14=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig104; 
title 'M09';
plot dm9*am9=2 dm9*am9=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig104; 
title 'M15';
plot dm15*am15=2 dm15*am15=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig104; 
title 'M06';
plot dm6*am6=2 dm6*am6=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig104; 
title 'M04';
plot dm4*am4=2 dm4*am4=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig104; 
title 'M01';
plot dm1*am1=2 dm1*am1=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig104; 
title 'M10';
plot dm10*am10=2 dm10*am10=1/overlay vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig104; 
title 'M16';
plot dm16*am16=2 dm16*am16=1/overlay vaxis=axis1 haxis=axis2;
run;
quit;

proc gplot data = fig104; 
title 'M05';
plot dm5*am5=2 dm5*am5=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig104; 
title 'M02';
plot dm2*am2=2 dm2*am2=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig104; 
title 'M11';
plot dm11*am11=2 dm11*am11=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig104; 
title 'M07';
plot dm7*am7=2 dm7*am7=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig104; 
title 'M08';
plot dm8*am8=2 dm8*am8=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig104; 
title 'M03';
plot dm3*am3=2 dm3*am3=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig104; 
title 'M12';
plot dm12*am12=2 dm12*am12=1/overlay vaxis=axis3 haxis=axis2;
run;
quit;

data orthodontwidem;
set orthodontwide;
if substr(subject,1,1)='M' then output;
run;
quit;

data orthodontwidem;
set orthodontwidem;
DistMAge8=Age8;
DistMAge10=Age10;
DistMAge12=Age12;
DistMAge14=Age14;
run;
quit;

*p. 340-341, fig. 10.5;
ods graphics on;
proc corr data=orthodontwidem plots=matrix;
var distmage14 distmage12 distmage10 distmage8;
run;
quit;
ods graphics off;

proc mixed data=male;
class subject;
model distance=age/solution intercept outp=mrF;
random subject;
run;
quit;

data res;
subject= sqrt(2.6407);
residual = sqrt(2.8164);
run;
quit;
proc print data=res;
run;
quit;

*fig. 10.6;
data fig106;
set mrf;
if subject='M01' then fm1 = pred;
if subject='M01' then dm1 = distance;
if subject='M02' then fm2 = pred;
if subject='M02' then dm2 = distance;
if subject='M03' then fm3 = pred;
if subject='M03' then dm3 = distance;
if subject='M04' then fm4 = pred;
if subject='M04' then dm4 = distance;
if subject='M05' then fm5 = pred;
if subject='M05' then dm5 = distance;
if subject='M06' then fm6 = pred;
if subject='M06' then dm6 = distance;
if subject='M07' then fm7 = pred;
if subject='M07' then dm7 = distance;
if subject='M08' then fm8 = pred;
if subject='M08' then dm8 = distance;
if subject='M09' then fm9 = pred;
if subject='M09' then dm9 = distance;
if subject='M10' then fm10 = pred;
if subject='M10' then dm10 = distance;
if subject='M11' then fm11 = pred;
if subject='M11' then dm11 = distance;
if subject='M12' then fm12 = pred;
if subject='M12' then dm12 = distance;
if subject='M13' then fm13 = pred;
if subject='M13' then dm13 = distance;
if subject='M14' then fm14 = pred;
if subject='M14' then dm14 = distance;
if subject='M15' then fm15 = pred;
if subject='M15' then dm15 = distance;
if subject='M16' then fm16 = pred;
if subject='M16' then dm16 = distance;
run;
quit;


data anno;
retain xsys ysys '1';
function = 'move'; x=0; y=0; output;
function = 'draw'; x=100; y=100; output;
run; 
quit;

goptions reset = all;
symbol1 v=circle c=blue;
axis1 label =(h=2 font=times angle=90 "Distance") 
value=(font=times h=1)
	  order=(16 to 32 by 2);
axis2 label =(h=2 font=times "Fitted") value =(font=times h=1)
      order=(20 to 32 by 2);
axis3 major=none minor=none label=NONE value=NONE order=(16 to 32 by 2);
axis4 major=none minor=none label=NONE value=NONE  order=(20 to 32 by 2);
	  goptions vsize=3;
	  goptions hsize=3;

proc gplot data = fig106 annotate=anno; 
title 'M06';
plot dm6*fm6=1 /vaxis=axis1 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M04';
plot dm4*fm4=1 /vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M01';
plot dm1*fm1=1 /vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M10';
plot dm10*fm10=1 /vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig106 annotate=anno; 
title 'M13';
plot dm13*fm13=1 /vaxis=axis1 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M14';
plot dm14*fm14=1 /vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M09';
plot dm9*fm9=1 /vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M15';
plot dm15*fm15=1 /vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig106 annotate=anno; 
title 'M07';
plot dm7*fm7=1 /vaxis=axis1 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M08';
plot dm8*fm8=1 /vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M03';
plot dm3*fm3=1 /vaxis=axis3 haxis=axis4;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M12';
plot dm12*fm12=1 /vaxis=axis3 haxis=axis4;
run;
quit;

proc gplot data = fig106 annotate=anno; 
title 'M16';
plot dm16*fm16=1 /vaxis=axis1 haxis=axis2;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M05';
plot dm5*fm5=1 /vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M02';
plot dm2*fm2=1 /vaxis=axis3 haxis=axis2;
run;
quit;
proc gplot data = fig106 annotate=anno; 
title 'M11';
plot dm11*fm11=1 /vaxis=axis3 haxis=axis2;
run;
quit;

data orthodont;
set orthodont;
fem=0;
if substr(subject,1,1)='F' then fem=1;
run;
quit;
*p. 343-344;
proc mixed data=orthodont;
class subject;
model distance=age fem fem*age/solution intercept;
random subject;
REPEATED / GROUP=fem;
run;
quit;

proc mixed data=orthodont;
class subject;
model distance=age fem fem*age/solution intercept;
random subject;
run;
quit;

data res;
subject = sqrt(3.2986);
residual = sqrt(1.9221);
run;
quit;
proc print data=res;
run;
quit;


%macro mixcomp;
ods trace on;
proc mixed ic data=orthodont;
ods output Mixed.InfoCrit=m105;
class subject;
model distance=age fem fem*age/solution intercept;
random subject;
REPEATED / GROUP=fem;
run;
quit;
ods trace off;
ods listing;

ods trace on;
proc mixed ic data=orthodont;
ods output Mixed.InfoCrit=m106;
class subject;
model distance=age fem fem*age/solution intercept;
random subject;
run;
quit;
ods trace off;
ods listing;

data m105;
set m105;
m105nll =  Neg2LogLike;
order=1;
run;
quit;

data m106;
set m106;
m106nll =  Neg2LogLike;
order=1;
run;
quit;

proc sort data=m105;
by order;
run;
quit;
proc sort data=m106;
by order;
run;
quit;

data together;
merge m105 m106;
by order;
run;
quit;

data lrpval;
set together;
teststat =  m106nll-m105nll	;
pval = 1-probchi(teststat,1);
keep teststat pval;
run;
quit;

proc print data=lrpval;
run;
quit;

%mend mixcomp;
*p. 345 p-value;
%mixcomp;



proc mixed data=orthodont;
class subject;
model distance=age fem fem*age/solution intercept outp=fitrand outpm=fitfix;
random subject/solution;
REPEATED / GROUP=fem;
run;
quit;

data fitrand;
set fitrand;
predrand=pred;
drop pred;
run;
quit;

data fitfix;
set fitfix;
predfix=pred;
drop pred;
run;
quit;

proc sort data=fitrand;
by subject age;
run;
quit;

proc sort data=fitfix;
by subject age;
run;
quit;

data together;
merge fitrand fitfix;
by subject age;
run;
quit;

data together;
set together ;
bayesred=predrand-predfix;
run;
quit;

proc means data=together noprint;
var bayesred;
by subject;
output out=together4qq MEAN=;
run;
quit; 


%macro qqplot(var=, dsn=);
 goptions reset = all htext=1.5;
  title1 height=2 font=times "Normal Q-Q";
  symbol1 value=circle color=black;
 proc univariate data = &dsn noprint;
  qqplot &var/normal(mu=est sigma=est l=1 color=black) 
      vminor=0 hminor=0 font=times
      vaxislabel= "&var";
 run;
quit;
%mend qqplot;

%qqplot(var=bayesred,dsn=together4qq);

data resfitfix;
set fitfix;
marginalres = distance-predfix;
keep subject age marginalres ;
run;
quit;


proc transpose data=resfitfix out=transfitfix prefix=MRAge;
by subject;
id age;
run;
quit;
*p. 347, fig. 10.8;
ods graphics on;
proc corr data=transfitfix;
var MRAge14 MRAge12 MRAge10 MRAge8;
run;
quit;
ods graphics off;

data resfitrand;
set fitrand;
condres = distance-predrand;
keep subject age condres;
run;
quit;


proc transpose data=resfitrand out=transfitrand prefix=CRAge;
by subject;
id age;
run;
quit;

*p. 348, fig. 10.9;
ods graphics on;
proc corr data=transfitrand;
var CRAge14 CRAge12 CRAge10 CRAge8;
run;
quit;
ods graphics off;

proc mixed data=orthodont;
class subject;
model distance=age fem fem*age/solution intercept outp=mod106;
random subject;
run;
quit;

data mod106;
set mod106;
condres=distance-pred;
keep subject pred condres;
run;
quit;

data malemod106;
set mod106;
if substr(subject,1,1)='M' then output;
run;
quit;

data femalemod106;
set mod106;
if substr(subject,1,1)='F' then output;
run;
quit;

goptions reset = all;
symbol1 v=circle c=blue;
axis1 label =(h=2 font=times angle=90 "Residuals (mm)") 
value=(font=times h=1)
	  order=(-5 to 5 by 2);

axis2 label =(h=2 font=times "Fitted Values (mm)") value =(font=times h=1)
      order=(16 to 32 by 2);
axis3 major=none minor=none label=NONE value=NONE order=(-5 to 5 by 2);

*fig. 10.10;
proc gplot data = malemod106; 
title 'Male';
plot condres*pred=1 /vaxis=axis1 haxis=axis2;
run;
quit;
proc gplot data = femalemod106; 
title 'Female';
plot condres*pred=1 /vaxis=axis3 haxis=axis2;
run;
quit;

proc mixed data=orthodont;
class subject;
model distance=age fem fem*age/solution intercept outp=mod105;
random subject;
REPEATED / GROUP=fem;
run;
quit;

data mod105;
set mod105;
condres=distance-pred;
keep subject pred condres;
run;
quit;

data malemod105;
set mod105;
if substr(subject,1,1)='M' then output;
run;
quit;

data femalemod105;
set mod105;
if substr(subject,1,1)='F' then output;
run;
quit;

goptions reset = all;
symbol1 v=circle c=blue;
axis1 label =(h=2 font=times angle=90 "Residuals (mm)") 
value=(font=times h=1)
	  order=(-6 to 6 by 2);
axis2 label =(h=2 font=times "Fitted Values (mm)") value =(font=times h=1)
      order=(16 to 32 by 2);
axis3 major=none minor=none label=NONE value=NONE order=(-6 to 6 by 2);

*fig. 10.11;
proc gplot data = malemod105; 
title 'Male';
plot condres*pred=1 /vaxis=axis1 haxis=axis2;
run;
quit;
proc gplot data = femalemod105; 
title 'Female';
plot condres*pred=1 /vaxis=axis3 haxis=axis2;
run;
quit;

*fig. 10.12;
data orthodont;
set orthodont;
male = ~fem;
run;
quit;

proc mixed data=orthodont;
class subject;
model distance=age fem fem*age/solution intercept outpm=mod106 VCIRY;
random subject;
run;
quit;

proc mixed data=orthodont;
class subject;
model distance=age fem fem*age/solution intercept outpm=mod105 VCIRY;
random subject;
REPEATED / GROUP=fem;
run;
quit;

%macro cholresboxplot(data=,model =);
goptions reset = all;
proc sort data = &data;
 by &var;
run;
quit;
proc gplot data = &data;
 axis1 label=(f=times h=2 "Sex") minor=none 
    order=(-1 0.2 0.8 2) value=(f=times h=1 
    t=1 ' ' t=2 'Female' t=3 'Male' t=4 ' ');
 axis2 label=(f=times h=2 angle=90 "Cholesky Residuals from Model &model") 
    value=(f=times h=1);
 symbol1 value = circle interpol=boxt bwidth=36;
 plot ScaledResid*male/ haxis=axis1 vaxis=axis2 vminor=0;
run;
quit;
%mend;                                            
%cholresboxplot(data=mod106,model=10.6);
%cholresboxplot(data=mod105,model=10.5);

proc import datafile="data/pigweights.csv" 
    out=pigs replace;
 getnames=yes;
run;
quit;

proc sort data=pigs;
by pigid weeknumber;
run;
quit;

* fig. 10.13;
%macro symbolfit(); 
     %do i = 1 %to 48;
symbol&i interpol=join
        value=circle;
		%end;
%mend;
axis1 label =(h=2 font=times angle=90 "Weight") 
value=(font=times h=1)
	  order=(20 to 90 by 10);
axis2 label =(h=2 font=times "Time") value =(font=times h=1)
      order=(1 to 9 by 1);
%symbolfit;

proc gplot data = pigs; 
plot weight*weeknumber=pigid/vaxis=axis1 haxis=axis2;
run;
quit;

data pigwide ;
set pigs;
keep pigid weeknumber weight;
run;
quit;

proc transpose data=pigwide out=pigwide let;
by pigid;
id weeknumber; 
run;
quit;

data pigwide;
set pigwide;
T1 = _1;
T2 = _2;
T3 = _3;
T4 = _4;
T5 = _5;
T6 = _6;
T7 = _7;
T8 = _8;
T9 = _9;
run;
quit;


*pg. 354, fig. 10.14;
ods graphics on;
proc corr data=pigwide plots=matrix;
var T1 T2 T3 T4 T5; 
run;
quit;
ods graphics off;

ods graphics on;
proc corr data=pigwide plots=matrix;
with  T1 T2 T3 T4 T5;
var T6 T7 T8 T9;
run;
quit;
ods graphics off;

ods graphics on;
proc corr data=pigwide plots=matrix;
with T6 T7 T8 T9;
var T1 T2 T3 T4 T5;
run;
quit;
ods graphics off;

ods graphics on;
proc corr data=pigwide plots=matrix;
var  T6 T7 T8 T9;
run;
quit;
ods graphics off;


*p. 356 model 10.13;
proc mixed data=pigs;
class pigid;
model weight=weeknumber/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;

*figure 10.16 ;
proc mixed data=pigs;
class pigid;
model weight=weeknumber/solution intercept outpm=mod10p13 VCIRY;
repeated /subject=pigid type=UN;
run;
quit;

*ods trace on;
*ods RESULTS on;

ods trace on;
proc mixed data=pigs;
ods output Mixed.R=c;
class pigid;
model weight=weeknumber/solution intercept outpm=mod10p13 VCIRY;
repeated /subject=pigid type=UN R;
run;
quit;
ods trace off;
ods listing;

proc iml;
*reset log print;
use mod10p13;
read all;
y = ScaledResid;
x =  weeknumber;
x =  J(nrow(x),1,1) || x ;
use c;
read all;
sn =  Col1 || Col2 || Col3 || Col4 || Col5 || Col6 || Col7 || Col8 ||
Col9;
sn = ginv(root((sn))`);
print sn;
 xp = x;

 do piglet = 1 to 48  ;
 xp[((piglet-1)*9+1):((piglet-1)*9+9)`,1] = sn *
 xp[((piglet-1)*9+1):((piglet-1)*9+9)`,1];
 xp[((piglet-1)*9+1):((piglet-1)*9+9)`,2] = sn *
 xp[((piglet-1)*9+1):((piglet-1)*9+9)`,2];
 end;
 print xp;
 x = xp[,2];
 create fig10p15 var{y x};
 append;

 quit;


*goptions reset = all;
*proc gplot data = fig10p15;
*symbol1 v=circle;
*plot y*x;
*run;
*quit;
 
proc loess data = fig10p15;
 model y=x/smooth=0.1;
 ods output OutputStatistics=loessout;
run;

data fit;
 set fig10p15;
 set loessout;

proc sort data = fit;
 by x;
run;

goptions reset = all;
symbol1 v=circle c=black;
symbol2 i=join c=black;
 axis1 label = (font=times h=2 angle=90 'Cholesky Residuals') 
      value=(font=times h=1);
 axis2 label = (font=times h=2 'x*') 
      value =(font=times h=1);
proc gplot data = fit;
 plot /*points:*/ y*x=1 /*loess:*/ 
     Pred*x=2/ overlay hminor=0 vminor=0 
     vaxis=axis1 haxis=axis2 vref=0;
run;


*do figure 10.16;
*do 10.16;
*do 10.16;
*do 10.16;
*do 10.16;
*do 10.16;
*do 10.16;
*do 10.16;
*do 10.16;
*do 10.16;
*do 10.16;


data pigspline;
set pigs;
Time = weeknumber;
TimeM2Plus = (Time>2)*(Time-2);
TimeM3Plus = (Time>3)*(Time-3);
TimeM4Plus = (Time>4)*(Time-4);
TimeM5Plus = (Time>5)*(Time-5);
TimeM6Plus = (Time>6)*(Time-6);
TimeM7Plus = (Time>7)*(Time-7);
TimeM8Plus = (Time>8)*(Time-8);
run;
quit;


*p 360, Model 10.15;
proc mixed data=pigspline;
class pigid;
model weight=Time TimeM2Plus TimeM3Plus TimeM4Plus TimeM5Plus TimeM6Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;

%macro mixcomp;
ods trace on;
proc mixed ic data=pigspline;
ods output Mixed.InfoCrit=m1015;
class pigid;
model weight=Time TimeM2Plus TimeM3Plus TimeM4Plus TimeM5Plus TimeM6Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;
ods trace off;
ods listing;

ods trace on;
proc mixed ic data=pigs;
ods output Mixed.InfoCrit=m1013;
class pigid;
model weight=weeknumber/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;
ods trace off;
ods listing;

data m1013;
set m1013;
m1013nll =  Neg2LogLike;
order=1;
run;
quit;

data m1015;
set m1015;
m1015nll =  Neg2LogLike;
order=1;
run;
quit;

proc sort data=m1015;
by order;
run;
quit;
proc sort data=m1013;
by order;
run;
quit;

data together;
merge m1015 m1013;
by order;
run;
quit;

data lrpval;
set together ;
teststat =  m1013nll-m1015nll	;
pval = 1-probchi(teststat,7);
keep teststat pval;
run;
quit;

proc print data=lrpval;
run;
quit;

%mend mixcomp;
*p. 361 p-value;
%mixcomp;

*p 363, Model 10.16;
proc mixed data=pigspline;
class pigid;
model weight=Time TimeM3Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;

*p. 363, p-value;
%macro mixcomp;
ods trace on;
proc mixed ic data=pigspline;
ods output Mixed.InfoCrit=m1015;
class pigid;
model weight=Time TimeM2Plus TimeM3Plus TimeM4Plus TimeM5Plus TimeM6Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;
ods trace off;
ods listing;

ods trace on;
proc mixed ic data=pigspline;
ods output Mixed.InfoCrit=m1016;
class pigid;
model weight=Time TimeM3Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;
ods trace off;
ods listing;

data m1016;
set m1016;
m1016nll =  Neg2LogLike;
order=1;
run;
quit;

data m1015;
set m1015;
m1015nll =  Neg2LogLike;
order=1;
run;
quit;

proc sort data=m1015;
by order;
run;
quit;
proc sort data=m1016;
by order;
run;
quit;

data together;
merge m1015 m1016;
by order;
run;
quit;

data lrpval;
set together ;
teststat =  m1016nll-m1015nll	;
pval = 1-probchi(teststat,4);
keep teststat pval;
run;
quit;

proc print data=lrpval;
run;
quit;
%mend;
%mixcomp;




*p 366, Model 10.17;
proc mixed data=pigspline;
class pigid;
model weight=Time TimeM3Plus TimeM5Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;



*p 366, pvalue;
%macro mixcomp;
ods trace on;
proc mixed ic data=pigspline;
ods output Mixed.InfoCrit=m1017;
class pigid;
model weight=Time TimeM3Plus TimeM5Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;
ods trace off;
ods listing;

ods trace on;
proc mixed ic data=pigspline;
ods output Mixed.InfoCrit=m1016;
class pigid;
model weight=Time TimeM3Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;
ods trace off;
ods listing;

data m1016;
set m1016;
m1016nll =  Neg2LogLike;
order=1;
run;
quit;

data m1017;
set m1017;
m1017nll =  Neg2LogLike;
order=1;
run;
quit;

proc sort data=m1017;
by order;
run;
quit;
proc sort data=m1016;
by order;
run;
quit;

data together;
merge m1017 m1016;
by order;
run;
quit;

data lrpval;
set together ;
teststat =  m1016nll-m1017nll	;
pval = 1-probchi(teststat,1);
keep teststat pval;
run;
quit;

proc print data=lrpval;
run;
quit;
%mend;
%mixcomp;



*p 367, Model 10.16 AR;
proc mixed data=pigspline;
class pigid;
model weight=Time TimeM3Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=AR(1) rcorr=1;
run;
quit;


*p 367, pvalue;
%macro mixcomp;
ods trace on;
proc mixed ic data=pigspline;
ods output Mixed.InfoCrit=m1016AR;
class pigid;
model weight=Time TimeM3Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=AR(1) rcorr=1;
run;
quit;
ods trace off;
ods listing;

ods trace on;
proc mixed ic data=pigspline;
ods output Mixed.InfoCrit=m1016;
class pigid;
model weight=Time TimeM3Plus TimeM7Plus TimeM8Plus/solution intercept;
repeated /subject=pigid type=UN rcorr=1;
run;
quit;
ods trace off;
ods listing;

data m1016;
set m1016;
m1016nll =  Neg2LogLike;
order=1;
run;
quit;

data m1016AR;
set m1016AR;
m1016ARnll =  Neg2LogLike;
order=1;
run;
quit;

proc sort data=m1016AR;
by order;
run;
quit;
proc sort data=m1016;
by order;
run;
quit;

data together;
merge m1016AR m1016;
by order;
run;
quit;

proc print data=together;
run;
quit;

data lrpval;
set together ;
teststat =  m1016ARnll-m1016nll;
pval = 1-probchi(teststat,45);
keep teststat pval;
run;
quit;

proc print data=lrpval;
run;
quit;
%mend;
%mixcomp;
