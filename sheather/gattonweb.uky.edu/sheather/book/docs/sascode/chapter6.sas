proc import datafile="data/nyc.csv" out=nyc replace;
 getnames=yes;
run;
quit;

*fig. 6.1;
ods graphics on;
proc corr data = nyc plots=matrix;
  var food decor service;
run;
quit;
ods graphics off;


proc reg data = nyc;
 model price = food decor service east;
 output out=outreg student=stdres p=fitted;
run;
quit;

*fig. 6.2;
goptions reset = all;
 axis1 label=(h=2 f=times 'Food') 
     value=(f=times h=1) offset=(3,3);
 axis2 label=(h=2 angle=90 f=times
     'Standardized Residuals') value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*food/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Decor') 
     value=(f=times h=1) offset=(3,0);
 axis2 label=(h=2 angle=90 f=times
     'Standardized Residuals') value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*decor/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Service') order=
     (14 to 24 by 2) value=(f=times h=1)
     offset=(3,3.5);
 axis2 label=(h=2 angle=90 f=times
     'Standardized Residuals') value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*service/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'East') order=
     (0 to 1 by 1) value=(f=times h=1)
     offset=(30,30);
 axis2 label=(h=2 angle=90 f=times
     'Standardized Residuals') value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*east/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

*fig. 6.3;
goptions reset = all;
 axis1 label=(h=2 f=times 'Fitted Values') 
     value=(f=times h=1);
 axis2 label=(h=2 angle=90 f=times 'Price') 
     value=(f=times h=1);
symbol1 value = circle i=r;
proc gplot data = outreg;
 plot price*fitted/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;

data gen;
 infile 'data/caution.txt';
 input x1 x2 y;
run;
quit;

*fig. 6.4;
ods graphics on;
proc corr data = gen plots=matrix;
  var y x1 x2;
run;
quit;
ods graphics off;


proc reg data = gen;
 model y = x1 x2;
 output out=outreg student=stdres p=fitted;
run;
quit;


*fig. 6.5;
goptions reset = all;
 axis1 label=(h=2 f=times 'x1') 
     value=(f=times h=1) offset=(0,3);
 axis2 label=(h=2 angle=90 f=times 
     'Studentized Residuals') 
     value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*x1/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'x2') 
     value=(f=times h=1) offset=(0,3);
 axis2 label=(h=2 angle=90 f=times 
     'Studentized Residuals') 
     value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*x2/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Fitted Values') 
     value=(f=times h=1) offset=(0,3);
 axis2 label=(h=2 angle=90 f=times 
     'Studentized Residuals') 
     value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*fitted/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;

*fig. 6.6;
goptions reset = all;
 axis1 label=(h=2 f=times 'Fitted Values') 
     value=(f=times h=1) offset=(0,3);
 axis2 label=(h=2 angle=90 f=times
     'y') order=(0 to 0.9 by 0.3)
     value=(f=times h=1);
symbol1 value = circle interpol=r;
proc gplot data = outreg;
 plot y*fitted/ haxis=axis1 vaxis=axis2 
     vminor=0 hminor=0;
run;
quit;


proc import datafile="data/nonlinearx.txt" out= nonlinearx replace;
 getnames=yes;
run;
quit;

*fig. 6.7;
goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'y')
	value=(f=times h=1);
axis2 label=(h=2 f=times 'x1')
	value=(f=times h=1);
symbol1 value = circle;
proc gplot data=nonlinearx;
plot y*x1/haxis=axis2 vaxis=axis1 ;
run;
quit;

goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'y' angle=90)
	value=(f=times h=1);
axis2 label=(h=2 f=times 'x2')
	value=(f=times h=1);
symbol1 value = circle;
proc gplot data=nonlinearx;
plot y*x2/haxis=axis2 vaxis=axis1 ;
run;
quit;

goptions reset=all;
axis1 label=(h=2 f=times 'x1')
	value=(f=times h=1);
axis2 label=(h=2 angle=90 f=times 'x2')
	value=(f=times h=1);
symbol1 value = circle;
proc gplot data=nonlinearx;
plot x2*x1/haxis=axis1 vaxis=axis2 ;
run;
quit;

proc reg data = nonlinearx;
 model y = x1 x2;
 output out=outreg student=stdres p=fitted;
run;
quit;

*fig 6.8;
goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'Standardized Residuals')
	value=(f=times h=1);
axis2 label=(h=2 f=times 'x1')
	value=(f=times h=1);
symbol1 value = circle;
proc gplot data=outreg;
plot stdres*x1/haxis=axis2 vaxis=axis1 ;
run;
quit;

goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'Standardized Residuals')
	value=(f=times h=1);
axis2 label=(h=2 f=times 'x2')
	value=(f=times h=1);
symbol1 value = circle;
proc gplot data=outreg;
plot stdres*x2/haxis=axis2 vaxis=axis1 ;
run;
quit;

goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'Standardized Residuals')
	value=(f=times h=1);
axis2 label=(h=2 f=times 'Fitted Values')
	value=(f=times h=1);
symbol1 value = circle;
proc gplot data=outreg;
plot stdres*fitted/haxis=axis2 vaxis=axis1 ;
run;
quit;

proc import datafile="data/nyc.csv" out=nyc replace;
 getnames=yes;
run;
quit;

*fig. 6.9;
goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'Price')
	value=(f=times h=1);
axis2 label=(h=2 f=times 'Food')
	value=(f=times h=1);
symbol1 value = circle interpol=r;
proc gplot data=nyc;
plot price*food/haxis=axis2 vaxis=axis1 ;
run;
quit;

goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'Price')
	value=(f=times h=1);
axis2 label=(h=2 f=times 'Decor')
	value=(f=times h=1);
symbol1 value = circle interpol=r;
proc gplot data=nyc;
plot price*decor/haxis=axis2 vaxis=axis1 ;
run;
quit;

goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'Price')
	value=(f=times h=1);
axis2 label=(h=2 f=times 'Service')
	value=(f=times h=1);
symbol1 value = circle interpol=r;
proc gplot data=nyc;
plot price*service/haxis=axis2 vaxis=axis1 ;
run;
quit;

goptions reset=all;
axis1 label=(h=2 angle=90 f=times 'Price')
	value=(f=times h=1);
axis2 label=(h=2 f=times 'East')
	value=(f=times h=1);
symbol1 value = circle interpol=r;
proc gplot data=nyc;
plot price*east/haxis=axis2 vaxis=axis1 ;
run;
quit;

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

*fig. 6.10;
%addvarplot(dsn=nyc, yvar=price, var1=Food, 
     ylab=price,othervar=decor service east);
%addvarplot(dsn=nyc, yvar=price, var1=decor,
 	ylab=price,othervar=food service east);
%addvarplot(dsn=nyc, yvar=price, var1=Service, 
     ylab=price,othervar=food decor east);
%addvarplot(dsn=nyc, yvar=price, var1=East, 
     ylab=price,othervar=food decor service);

data defects;
 infile 'data/defects.txt' firstobs=2 expandtabs;
 input case temperature density rate defective;
run;
quit;

*fig. 6.11;
ods graphics on;
proc corr data = defects plots=matrix;
  var defective temperature density rate;
run;
quit;
ods graphics off;

proc reg data = defects noprint;
 model defective = temperature density rate;
 output out = outreg student=stdres p=fitted;
run;
quit;

*fig. 6.12;
goptions reset = all;
 axis1 label=(h=2 f=times 'Temperature') 
     value=(f=times h=1) offset=(3,3)
     order=(0.5 to 3.5 by 1);
 axis2 label=(h=2 angle=90 f=times
     'Standardized Residuals') value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*temperature/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Density') 
     value=(f=times h=1) offset=(3,3)
     order=(18 to 33 by 3);
 axis2 label=(h=2 angle=90 f=times
     'Standardized Residuals') value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*density/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Rate') 
     value=(f=times h=1) offset=(3,3)
     order=(175 to 300 by 25);
 axis2 label=(h=2 angle=90 f=times
     'Standardized Residuals') value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*rate/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Fitted Values') 
     value=(f=times h=1) offset=(3,3);
 axis2 label=(h=2 angle=90 f=times
     'Standardized Residuals') value=(f=times h=1);
symbol1 value = circle;
proc gplot data = outreg;
 plot stdres*fitted/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

data outreg;
 set outreg;
 fitsq = fitted **2;
run;
quit;

proc reg data = outreg noprint;
 model defective=fitted fitsq;
 output out=plotit p=newfit;
run;
quit;

proc sort data = plotit;
 by fitted;
run;
quit;

*fig. 6.13; 
goptions reset = all;
 axis1 label=(h=2 f=times 'Fitted Values') 
     value=(f=times h=1) offset=(3,3);
 axis2 label=(h=2 angle=90 f=times
     'Defective') value=(f=times h=1)
     order=(-5 to 65 by 10);
symbol1 value = circle i=r c=black l=4;
symbol2 i=join c=black l=1;
proc gplot data = plotit;
 plot defective*fitted newfit*fitted/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0 overlay;
run;
quit;

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


call nlpqn(rc,lambdares,"irp_call",lambda,optn,);
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
quit;

%macro regouts(dsn=,yvar=,predname=);
 proc reg data = invresplot noprint;
  model fitted = &yvar;
  output out= &dsn p=&predname;
 run;
 quit;
%mend;

%regouts(dsn=data1,yvar=y, predname=lam1hat);
%regouts(dsn=data2,yvar=cbrty,predname=cbrtyhat);
%regouts(dsn=data3,yvar=ly,predname=lyhat);

proc sort data = invresplot;
 by y;
run;
quit;

%macro sortit;
%do i = 1 %to 3;
 proc sort data = data&i;
  by y;
 run;
 quit;
%end;
%mend sortit;
%sortit;

data full;
 merge invresplot data1 data2 data3;
 by y;
run;
quit;


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
quit;
%mend irp;

*fig. 6.14;
%irp(resp=defective, fit=fitted,dat=outreg);

*fig. 6.15;
ods output boxcox=b details=d; 
   ods exclude boxcox; 
proc transreg details data = defects;
 model boxcox(defective/ convenient 
    lambda=0.3 to 0.65 by .001)=identity(density temperature rate);
 output out=trans;
run;
quit;

data _null_; 
 set d; 
 if description = 'CI Limit' 
     then call symput('vref',   formattedvalue); 
 if description = 'Lambda Used' 
     then call symput('lambda', formattedvalue); 
run; 
quit;

proc print data = b;
run;
quit;

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
quit;

goptions reset = all;
axis2 order=(0.3 to 0.65 by 0.05) 
     label=(f=times h=2 "Lambda") value=(h=1 
     f=times);
axis1  label=(angle=90 
     f=times h=2 "log-Likelihood") value=(h=1
     f=times);
proc gplot data = b;
 plot loglike * lambda / vref=&vref href=&href1 &href2 &href3 
     vminor=0 hminor=0 vaxis=axis1 haxis=axis2; 
 symbol v=none i=spline c=black; 
run; 
quit;

data transf;
 set defects;
 rty = sqrt(defective);
run;
quit;

%macro plotit(xvar=);
goptions reset = all;
axis2 label=(f=times h=2 "&xvar") value=(h=1 
     f=times);
axis1 order=(0 to 8 by 2) label=(angle=90 
     f=times h=2 "Sqrt(Defective)") value=(h=1
     f=times);
symbol v=circle; 
proc gplot data = transf;
plot rty*&xvar /vminor=0 hminor=0 vaxis=axis1
     haxis=axis2; 
run; 
quit;
%mend plotit;

*fig. 6.16;
%plotit(xvar=Temperature);
%plotit(xvar=Density);
%plotit(xvar=Rate);

*p. 175;
proc reg data = transf;
 model rty = temperature density rate;
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

*fig. 6.17;
%stdresplot(xvar=Temperature,xlab=Temperature);
%stdresplot(xvar=Density,xlab=Density);
%stdresplot(xvar=Rate,xlab=Rate);
%stdresplot(xvar=fitted,xlab=Fitted Values);

*fig. 6.18;
goptions reset = all;
axis2 label=(f=times h=2 "Fitted Values") value=(h=1 
     f=times) ;
axis1 label=(angle=90 f=times h=2 
     "Sqrt(Defective)") value=(h=1 f=times);
symbol v=circle i=r; 
proc gplot data = outreg;
plot rty*fitted/vminor=0 hminor=0 
     vaxis=axis1 haxis=axis2; 
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

*fig. 6.19;
%plotlm(regout=outreg);

*fig. 6.20;
%addvarplot(dsn=transf, yvar=rty, var1=Temperature, 
     ylab=sqrt(Defective),
     othervar=Density Rate);
%addvarplot(dsn=transf, yvar=rty, var1=Density,
     ylab=sqrt(Defective), 
	 othervar=Temperature Rate);
%addvarplot(dsn=transf, yvar=rty, var1=Rate,
     ylab=sqrt(Defective),
	 othervar=Temperature Density);

proc import datafile='data/magazines.csv' out=ads
     replace;
 getnames=yes;
run;
quit;

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

*p. 177;
proc iml;
  use ads;
  read all ;
  namesX={"AdPages" "SubRevenue" "NewsRevenue"};
X  = adpages || subrevenue || newsrevenue ;
%mboxcox;
run;
quit;

*fig. 6.21;
ods graphics on;
proc corr data = ads plots=matrix;
  var adrevenue adpages subrevenue newsrevenue;
run;
quit;
ods graphics off;

data ads;
set ads;
ladpages = log(adpages);
lsubrevenue = log(subrevenue);
lnewsrevenue=log(newsrevenue);
run;
quit;

proc reg data = ads noprint;
 model adrevenue = ladpages lsubrevenue lnewsrevenue;
 output out = outreg p=fitted student=stdres;
run;
quit;

data invresplot;
set outreg;
y=adrevenue	;
cbrty=y**(.23);
ly=log(y) ;
run;
quit;


%macro regouts(dsn=,yvar=,predname=);
 proc reg data = invresplot noprint;
  model fitted = &yvar;
  output out= &dsn p=&predname;
 run;
quit;
%mend;

%regouts(dsn=data1,yvar=y, predname=lam1hat);
%regouts(dsn=data2,yvar=cbrty,predname=cbrtyhat);
%regouts(dsn=data3,yvar=ly,predname=lyhat);

proc sort data = invresplot;
 by y;
run;
quit;

%macro sortit;
%do i = 1 %to 3;
 proc sort data = data&i;
  by y;
 run;
quit;
%end;
%mend sortit;
%sortit;

data full;
 merge invresplot data1 data2 data3;
 by y;
run;
quit;

*fig. 6.22;
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
     j=c ".23" j=c '0' j=c);
proc gplot data = full;
 plot fitted*y=1 lam1hat*y=2 cbrtyhat*y=3 lyhat*y=4/
     overlay vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0 legend=legend1;
run;
quit;

*p. 179;
proc iml;
  use ads;
  read all ;
  namesX={"AdRevenue" "AdPages" "SubRevenue" "NewsRevenue"};
X  = adrevenue || adpages || subrevenue || newsrevenue ;
%mboxcox;
run;
quit;

data ads;
set ads;
ladrevenue=log(adrevenue);
run;
quit;

*fig. 6.23;
ods graphics on;
proc corr data = ads plots=matrix;
  var ladrevenue ladpages lsubrevenue lnewsrevenue;
run;
quit;
ods graphics off;

*p. 183;
proc reg data = ads;
 model ladrevenue = ladpages lsubrevenue lnewsrevenue;
output out=outreg r=resids p=fitted student=stdres cookd=cd h=levg;
run;
quit;

*fig. 6.24;
%stdresplot(xvar=lAdPages,xlab=log(AdPages));
%stdresplot(xvar=lSubRevenue,xlab=log(SubRevenue));
%stdresplot(xvar=lNewsRevenue,xlab=log(NewsRevenue));
%stdresplot(xvar=fitted,xlab=Fitted Values);

*fig. 6.25;
goptions reset = all;
symbol1 v=circle i=r;
axis1 label=(h=2 angle=90 f=times 
    "log(AdRevenue)") value=(h=1 f=times);
axis2 label=(h=2 f=times "Fitted Values") 
    value=(h=1 f=times);
proc gplot data = outreg;
 plot ladrevenue*fitted/
     vaxis=axis1 haxis=axis2 vminor=0 
     hminor=0;
run;
quit;

*fig. 6.26;
%plotlm(regout=outreg);

*fig. 6.27;
%addvarplot(dsn=ads, yvar=ladrevenue, var1=ladpages, 
     ylab=log(AdRevenue),
     othervar=lsubrevenue lnewsrevenue);
%addvarplot(dsn=ads, yvar=ladrevenue, var1= lsubrevenue, 
     ylab=log(AdRevenue),
     othervar= ladpages lnewsrevenue);
%addvarplot(dsn=ads, yvar=ladrevenue, var1= lnewsrevenue, 
     ylab=log(AdRevenue),
     othervar= ladpages lsubrevenue);

proc import datafile="data/circulation.txt" 
    out=news replace;
run;
quit;

data logs;
 set news;
 lsun = log(sunday);
 lweek = log(weekday);
run;
quit;

*fig. 6.28;
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

proc reg data = logs;
 model lsun = lweek tabloid;
output out=outreg r=resids p=fitted student=stdres cookd=cd h=levg;
run;
quit;

*fig. 6.29;
%stdresplot(xvar=lweek,xlab=log(Weekday Circulation));
%stdresplot(xvar=tabloid,xlab=Tabloid.with.a.Serious.Competitor);
%stdresplot(xvar=fitted,xlab=Fitted Values);

*fig. 6.30;
goptions reset = all;
axis1 label=(f=times h=2 angle=90 
    'log(Sunday Circulation)') value=(f=times h=1);
axis2 value=(f=times h=1) order=(11.5 to 14.5 by 1)
    label =(f=times h=2 'Fitted Values');
symbol1 v = circle i=r;
proc gplot data = outreg;
 plot lsun*fitted/haxis = axis2 hminor=0 
     vaxis=axis1 vminor=0;
run;
quit;

*fig. 6.31;
%plotlm(regout=outreg);

data new;
 input weekday tabloid;
 cards;
 210000 1
 210000 0
;
run;
quit;

data new;
set new;
lweek = log(weekday);
run;
quit;

data newnews;
set new logs;
run;
quit;

*p. 188;
proc reg data = newnews noprint;
 model lsun = lweek tabloid;
 output out=outreg p=fit lcl= lwr ucl=upr;
run;
quit;

data outreg;
set outreg;
if weekday in(210000) then output;
run;
quit;

proc print data=outreg noobs;
var fit lwr upr;
run;
quit;

*fig. 6.32;
%addvarplot(dsn=logs, yvar=lsun, var1=lweek, 
     ylab=log(Sunday),
     othervar=tabloid);
%addvarplot(dsn=logs, yvar=lsun, var1=tabloid, 
     ylab=log(Sunday),
     othervar=lweek);

data profsal; 
 infile 'data/profsalary.txt' firstobs=2 expandtabs;
 input case salary exper;
run;
quit;

proc loess data = profsal;
 model salary=exper/smooth=0.6667;
 ods output OutputStatistics=loessout;
run;
quit;

data fit;
 set profsal;
 set loessout;
run;
quit;

proc sort data = fit;
 by exper;
run;
quit;

*fig. 6.33;
goptions reset = all;
symbol1 v=circle c=black i=r l=2;
symbol2 i=join c=black l=1;
 axis1 label=(h=2 font=times angle=90 
     "Salary") 
     value=(font=times h=1);
 axis2 value =(font=times h=1)
     label=(h=2 font=times 'Years of Experience'); 
proc gplot data = fit;
 plot /*points:*/ salary*exper=1 /*loess:*/ Pred*exper=2/
     overlay hminor=0 vminor=0 vaxis=axis1 haxis=axis2;
run;
quit;

*fig. 6.34;
goptions reset = all;
symbol1 v=circle c=red i=rq l=2;
symbol2 i=join c=blue l=1;
 axis1 label=(h=2 font=times angle=90 
     "Salary") 
     value=(font=times h=1);
 axis2 value =(font=times h=1)
     label=(h=2 font=times 'Years of Experience'); 
proc gplot data = fit;
 plot /*points:*/ salary*exper=1 /*loess:*/ Pred*exper=2/
     overlay hminor=0 vminor=0 vaxis=axis1 haxis=axis2;
run;
quit;

*fig. 6.35;
proc loess data = defects;
 model defective=temperature/smooth=0.6667;
 ods output OutputStatistics=loessout1;
run;
quit;

data fit;
 set defects;
 set loessout1;
run;
quit;

proc sort data = fit;
 by temperature;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Temperature, x1') 
     value=(f=times h=1) order=(.75 to 3.25 by .5);
 axis2 label=(h=2 angle=90 f=times
     'Defective, Y') value=(f=times h=1);
symbol1 v = circle c=black;
symbol2 i=join c=black;
proc gplot data = fit;
 plot /*Points:*/ defective*temperature=1 /*loess:*/
     Pred*temperature=2/ haxis=axis1 vaxis=axis2
     hminor=0 vminor=0 overlay;
run;
quit;

proc reg data = defects noprint;
 model defective = temperature density rate;
 output out = outreg p=yhat;
run;quit;


proc loess data = outreg;
 model yhat=temperature/smooth=0.6667;
 ods output OutputStatistics=loessout;
run;
quit;

data fit;
 set outreg;
 set loessout;
run;
quit;

proc sort data = fit;
 by temperature;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Temperature, x1') 
     value=(f=times h=1) order=(.75 to 3.25 by .5);
 axis2 label=(h=2 angle=90 f=times
     'Y-hat') value=(f=times h=1);
symbol1 v = circle c=black;
symbol2 i=join c=black;
proc gplot data = fit;
 plot /*points:*/ yhat*temperature=1 /*loess:*/
     Pred*temperature=2/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0 overlay;
run;
quit;

*fig. 6.36;

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
 by temperature;
run;
quit;

goptions reset = all;
 axis1 label=(h=2 f=times 'Temperature, x1') 
     value=(f=times h=1) order=(.75 to 3.25 by .5);
 axis2 label=(h=2 angle=90 f=times
     'Defective') value=(f=times h=1)order=
     (-5 to 65 by 10);
symbol1 v = circle c=black;
symbol2 i=join c=blue;
symbol3 i=join c=red l=2;
proc gplot data = fit;
 plot /*points:*/ defective*temperature=1 /*loess:*/
     Pred*temperature=2 Pred2*temperature=3/ haxis=axis1 
     vaxis=axis2 vminor=0 hminor=0 overlay;
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

*fig. 6.37;
%mmplot(dsn=defects, yvar=defective, x1=temperature, 
     allx= temperature density rate, xlab=Temperature,
     ylab=Defective);
%mmplot(dsn=defects, yvar=defective, x1=density, 
     allx= temperature density rate, xlab=Density,
     ylab=Defective);
%mmplot(dsn=defects, yvar=defective, x1=rate, 
     allx= temperature density rate, xlab=Rate,
     ylab=Defective);
%mmplot(dsn=outreg, yvar=defective, x1=yhat, 
     allx= temperature density rate, xlab=Fitted Values,
     ylab=Defective);

proc reg data = transf noprint;
 model rty = temperature density rate;
 output out=outreg p=fitted;
run;
quit;

*fig. 6.38;
%mmplot(dsn=outreg, yvar=rty, x1=temperature, 
     allx= temperature density rate, xlab=Temperature,
     ylab=Sqrt(defective));
%mmplot(dsn=outreg, yvar=rty, x1=density, 
     allx= temperature density rate, xlab=Density,
     ylab=Sqrt(defective));
%mmplot(dsn=outreg, yvar=rty, x1=rate, 
     allx= temperature density rate, xlab=Rate,
     ylab=Sqrt(defective));
%mmplot(dsn=outreg, yvar=rty, x1=fitted, 
     allx= temperature density rate, xlab=Fitted Values,
     ylab=Sqrt(defective));

data bridge;
 infile 'data/bridge.txt' expandtabs firstobs=2;
 input case time darea ccost dwgs length spans;
run;
quit;

*p. 196;
proc iml;
  use bridge;
  read all ;
  namesX={"Time" "Darea" "CCost" "Dwgs" "Length" "Spans" };
X  = time || darea || ccost || dwgs || length || spans ;
%mboxcox;
run;
quit;

*fig. 6.39;
ods graphics on;
proc corr data = bridge plots=matrix;
  var Time DArea CCost Dwgs Length Spans;
run;
quit;
ods graphics off;

ods graphics on;
proc corr data = bridge plots=matrix;
  with spans;
  var Time DArea CCost Dwgs Length;
run;
quit;
ods graphics off;


data tbridge;
 set bridge;
 log_Time = log(time);
 log_darea = log(darea);
 log_ccost = log(ccost);
 log_dwgs = log(dwgs);
 log_length = log(length);
 log_spans = log(spans);
run;
quit;

*fig. 6.40;
ods graphics on;
proc corr data = tbridge plots=matrix;
  var log_Time log_DArea log_CCost log_Dwgs 
     log_Length log_Spans;
run;
quit;
ods graphics off;
ods graphics on;
proc corr data = tbridge plots=matrix;
  with log_Spans;
  var log_Time log_DArea log_CCost log_Dwgs 
     log_Length ;
run;
quit;
ods graphics off;

proc reg data = tbridge noprint;
 model log_time = log_darea log_ccost log_dwgs 
     log_length log_spans;
output out=outreg r=resids p=fitted student=stdres cookd=cd h=levg;
run;
quit;

*fig. 6.41;
%stdresplot(xlab=log(DArea), xvar=log_darea);
%stdresplot(xlab=log(CCost), xvar = log_ccost);
%stdresplot(xlab=log(Dwgs), xvar = log_Dwgs);
%stdresplot(xlab=log(Length), xvar = log_length);
%stdresplot(xlab = log(Spans), xvar = log_length);
%stdresplot(xlab = Fitted Values, xvar = fitted);

*fig. 6.42;
goptions reset = all;
 axis1 label=(h=2 f=times "Fitted Values") 
     value=(f=times h=1) order=(3.5 to 6 by 0.5);
 axis2 label=(h=2 angle=90 f=times 'log(Time)')
     value=(f=times h=1) order=(3.5 to 6.5 by 1);
symbol1 value = circle i=r;
proc gplot data = outreg;
 plot log_time*fitted/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

*fig. 6.43;
%plotlm(regout=outreg);

*p. 200;
proc reg data = tbridge;
 model log_time = log_darea log_ccost log_dwgs 
     log_length log_spans;
run;
quit;

*fig. 6.44;
%mmplot(dsn=outreg, yvar=log_time, x1=log_darea, 
     allx= log_darea log_ccost log_dwgs log_length
     log_spans, xlab=log(Darea), ylab=log(Time)); 
%mmplot(dsn=outreg, yvar=log_time, x1=log_ccost, 
     allx= log_darea log_ccost log_dwgs log_length
     log_spans, xlab=log(CCost), ylab=log(Time)); 
%mmplot(dsn=outreg, yvar=log_time, x1=log_dwgs, 
     allx= log_darea log_ccost log_dwgs log_length
     log_spans, xlab=log(Dwgs), ylab=log(Time)); 
%mmplot(dsn=outreg, yvar=log_time, x1=log_length, 
     allx= log_darea log_ccost log_dwgs log_length
     log_spans, xlab=log(Length), ylab=log(Time)); 
%mmplot(dsn=outreg, yvar=log_time, x1=log_spans, 
     allx= log_darea log_ccost log_dwgs log_length
     log_spans, xlab=log(Spans), ylab=log(Time)); 
%mmplot(dsn=outreg, yvar=log_time, x1=fitted, 	
     allx= log_darea log_ccost log_dwgs log_length
     log_spans, xlab=Fitted Values, ylab=log(Time)); 

*p 202;
proc corr data = tbridge noprob nosimple;
  var log_ccost log_darea log_dwgs log_length
   log_spans;
run;
quit;

*fig. 6.45;
%addvarplot(dsn=tbridge, yvar=log_time, 
     var1=log_DArea, ylab=
	 log(Time), othervar=log_ccost log_dwgs 
	 log_length log_spans);
%addvarplot(dsn=tbridge, yvar=log_time, 
     var1=log_CCost, ylab=
	 log(Time), othervar=log_darea log_dwgs 
	 log_length log_spans);
%addvarplot(dsn=tbridge, yvar=log_time, 
     var1=log_Dwgs, ylab=
	 log(Time), othervar=log_ccost log_darea 
	 log_length log_spans);
%addvarplot(dsn=tbridge, yvar=log_time, 
     var1=log_Length, ylab=
	 log(Time), othervar=log_ccost log_darea 
	 log_dwgs log_spans);
%addvarplot(dsn=tbridge, yvar=log_time, 
     var1=log_Spans, ylab=
	 log(Time), othervar=log_ccost log_darea 
	 log_dwgs log_length);

*p. 203;
proc reg data = tbridge;
 model log_time = log_darea log_ccost log_dwgs 
   log_length log_spans/ vif;
run;
quit;

proc import datafile="data/Bordeaux.csv" out=wine 
     replace;
 getnames=yes;
run;
quit;	

data wine;
 set wine;
 log_price = log(price);
 log_ParkerPoints= log(ParkerPoints);
 log_CoatesPoints = log(CoatesPoints);
run;
quit;

*p. 206-207;
proc reg data= wine;
model log_price = log_ParkerPoints log_CoatesPoints P95andAbove FirstGrowth 
CultWine Pomerol VintageSuperstar/ vif;
output out=outreg r=resids p=fitted student=stdres cookd=cd h=levg;
run;	
quit;

%macro boxplot(xvar =,xlab=);
proc sort data = outreg;
 by &xvar;
run;
quit;
proc gplot data = outreg;
 axis1 label=(f=times h=2 "&xlab") minor=none 
    order=(-1 0.2 0.8 2) value=(f=times h=1 
    t=1 ' ' t=2 '0' t=3 '1' t=4 ' ');
 axis2 label=(f=times h=2 angle=90 "Standardized Residuals") 
    value=(f=times h=1);
 symbol1 value = circle interpol=boxt bwidth=36;
 plot stdres*&xvar/ haxis=axis1 vaxis=axis2 vminor=0;
run;
quit;
%mend;

*fig. 6.46;
%stdresplot(xvar=log_ParkerPoints,xlab=log(Parker Points));
%stdresplot(xvar=log_CoatesPoints,xlab=log(Coates Points));
%boxplot(xvar=P95andAbove,xlab=P95andAbove); 
%boxplot(xvar=FirstGrowth,xlab=FirstGrowth); 
%boxplot(xvar=CultWine,xlab=CultWine); 
%boxplot(xvar=Pomerol,xlab=Pomerol); 
%boxplot(xvar=VintageSuperstar,xlab=VintageSuperstar);
%stdresplot(xvar=fitted,xlab=Fitted Values);	

*fig. 6.47;
goptions reset = all;
 axis1 label=(h=2 f=times "Fitted Values") 
     value=(f=times h=1); 
 axis2 label=(h=2 angle=90 f=times 'log(Price)')
     value=(f=times h=1);
symbol1 value = circle i=r;
proc gplot data = outreg;
 plot log_price*fitted/ haxis=axis1 vaxis=axis2 vminor=0
     hminor=0;
run;
quit;

*fig. 6.48;
%plotlm(regout=outreg);

*fig. 6.49;
%mmplot(dsn=outreg, yvar=log_price, x1=log_parkerpoints, 
     allx= log_ParkerPoints log_CoatesPoints P95andAbove FirstGrowth CultWine Pomerol VintageSuperstar, xlab=log(ParkerPoints), ylab=log(Price)); 
%mmplot(dsn=outreg, yvar=log_price, x1=log_coatespoints, 
     allx= log_ParkerPoints log_CoatesPoints P95andAbove FirstGrowth CultWine Pomerol VintageSuperstar, xlab=log(CoatesPoints), ylab=log(Price)); 
%mmplot(dsn=outreg, yvar=log_price, x1=fitted, 
     allx= log_ParkerPoints log_CoatesPoints P95andAbove FirstGrowth CultWine Pomerol VintageSuperstar, xlab=Fitted Values, ylab=log(Price)); 

*fig. 6.50;
%addvarplot(dsn=wine, yvar=log_price, 
     var1=log_parkerpoints, ylab=
	 log(Price), othervar= log_CoatesPoints P95andAbove FirstGrowth CultWine Pomerol VintageSuperstar);
%addvarplot(dsn=wine, yvar=log_price, 
     var1=log_coatespoints, ylab=
	 log(Price), othervar= log_ParkerPoints P95andAbove FirstGrowth CultWine Pomerol VintageSuperstar);
%addvarplot(dsn=wine, yvar=log_price, 
     var1=P5andAbove, ylab=
	 log(Price), othervar= log_ParkerPoints log_CoatesPoints FirstGrowth CultWine Pomerol VintageSuperstar);
%addvarplot(dsn=wine, yvar=log_price, 
     var1=firstgrowth, ylab=
	 log(Price), othervar= log_ParkerPoints log_CoatesPoints P95andAbove CultWine Pomerol VintageSuperstar);
%addvarplot(dsn=wine, yvar=log_price, 
     var1=cultwine, ylab= log(Price), othervar= log_ParkerPoints log_CoatesPoints P95andAbove FirstGrowth Pomerol VintageSuperstar);
%addvarplot(dsn=wine, yvar=log_price, 
     var1=pomerol, ylab=
	 log(Price), othervar= log_ParkerPoints log_CoatesPoints P95andAbove FirstGrowth CultWine VintageSuperstar);
%addvarplot(dsn=wine, yvar=log_price, 
     var1=vintagesuperstar, ylab=
	 log(Price), othervar= log_ParkerPoints log_CoatesPoints P95andAbove FirstGrowth CultWine Pomerol);

*p. 208-209;
proc reg data= wine;
model log_price = log_ParkerPoints log_CoatesPoints FirstGrowth CultWine Pomerol VintageSuperstar;
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

*p. 209;
%ftest(dsn=wine, yvar=log_price, fullm= log_ParkerPoints log_CoatesPoints P95andAbove FirstGrowth CultWine Pomerol VintageSuperstar, reducedm= log_ParkerPoints log_CoatesPoints FirstGrowth CultWine Pomerol VintageSuperstar, fulldf=64, reddf=65);

*tables 6.4 and 6.5;
data table6p4;
set outreg;
if stdres > 2.4 then output;
keep wine stdres;
run;
quit;

data table6p5;
set outreg;
if stdres < -2.7 then output;
keep wine stdres;
run;
quit;

proc print data=table6p4;
run;
quit;

proc print data=table6p5;
run;
quit;


proc import datafile="data/storks.txt" out= storks replace;
 getnames=yes;
run;
quit;
	
*fig. 6.51;
goptions reset = all;
axis2 label=(f=times h=2 "Number of Storks") value=(h=1 
     f=times);
axis1 label=(angle=90 
     f=times h=2 "Number of Babies") value=(h=1
     f=times);
symbol v=circle i=r; 
proc gplot data = storks;
plot babies*storks /vaxis=axis1
     haxis=axis2; 
run;
quit;

*p. 212; 
proc reg data=storks;
model babies=storks;
run;
quit;

*fig. 6.52;
goptions reset = all;
axis2 label=(f=times h=2 "Number of Storks") value=(h=1 
     f=times);
axis1 label=(angle=90 
     f=times h=2 "Number of Babies") value=(h=1
     f=times);
symbol v=circle i=r; 
proc gplot data = storks;
plot babies*storks /vaxis=axis1
     haxis=axis2; 
run;
quit;

goptions reset = all;
axis2 label=(f=times h=2 "Number of Women") value=(h=1 
     f=times);
axis1 label=(angle=90 
     f=times h=2 "Number of Babies") value=(h=1
     f=times);
symbol v=circle i=r; 
proc gplot data = storks;
plot babies*women /vaxis=axis1
     haxis=axis2; 
run;
quit;

goptions reset = all;
axis2 label=(f=times h=2 "Number of Storks") value=(h=1 
     f=times);
axis1 label=(angle=90 
     f=times h=2 "Number of Women") value=(h=1
     f=times);
symbol v=circle i=r; 
proc gplot data = storks;
plot women*storks /vaxis=axis1
     haxis=axis2; 
run;
quit;

*p. 213;
proc reg data=storks;	
model babies=women storks;
run;
quit;
