* fishing.line.sas, fishing line experiment, Table 19.33, p760;
;
* Data as arranged in Table 19.33;
data fish;
  input reel brand @@;
  do splitplot=1 to 4;
    input B $ weight @@;
	if B='S' then stress=1; else stress=2;
	output; drop B;
  end;
  lines;
  1 1  N 6.70  S 6.40  S 7.20  N 7.00 
  2 2  S 8.10  S 8.90  N 8.00  N 6.10 
  3 2  S 8.00  S 8.00  N 8.75  N 8.50 
  4 1  N 8.50  S 9.50  N 9.70  S 9.40 
;
* Data by observation;
data fish;
  input reel brand stress weight;
  * stress: 1 = stressed (S), 2 = non-stressed (N);
  lines;
          1    1      2    6.70
          1    1      1    6.40
          1    1      1    7.20
          1    1      2    7.00
          2    2      1    8.10
          2    2      1    8.90
          2    2      2    8.00
          2    2      2    6.10
          3    2      1    8.00
          3    2      1    8.00
          3    2      2    8.75
          3    2      2    8.50
          4    1      2    8.50
          4    1      1    9.50
          4    1      2    9.70
          4    1      1    9.40
;
run;
