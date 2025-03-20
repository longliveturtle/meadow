readme.txt : Stata_Code_Data_Graphics.zip

[current versions of these files available at 
	http://www.stat.tamu.edu/~sheather/book/]

Contents of zip file

data directory
	contains 40 files
	adcosts.csv - travel.txt
graphics directory
	contains 170 files
	fAPPp1.eps - f10p14.eps

main directory
	
	contains 11 .do files
	appendix.do
	chapter1.do
	chapter2.do
	chapter3.do
	chapter4.do
	chapter5.do
	chapter6.do
	chapter7.do
	chapter8.do
	chapter9.do
	chapter10.do
	
	contains 12 .ado files
	analysis_deviance.ado
	draw_matrix.ado
	irp.ado
	leaps_and_bounds.ado
	lowess_ties_optim.ado
	mboxcox.ado
	mmp.ado
	plot_bc.ado
	plot_lm.ado
	rlcipi.ado
	step_crit.ado
	vselect.ado
	
	contains 12 .sthlp files
	analysis_deviance.sthlp
	draw_matrix.sthlp
	irp.sthlp
	leaps_and_bounds.sthlp
	lowess_ties_optim.sthlp
	mboxcox.sthlp
	mmp.sthlp
	plot_bc.sthlp
	plot_lm.sthlp
	rlcipi.sthlp
	step_crit.sthlp
	vselect.sthlp


After unzipping the files, open Stata and set the current working 
directory to the main directory of the files 
(where the .do files are stored).

The analysis from the book can be completed by typing at the terminal window

do appendix.do
do chapter"X".do

where "X" = 1,...,10.

Selected parts of the analysis may be completed by opening the .do files using 
Stata's .do file edittor and highlighting commands, then pressing the "do" button 
within the edittor.

The same result may be obtained by copying the commands and pasting them into the
Stata terminal window.  This can be done using any text edittor.

The .ado files contain original programs that the book analysis uses.  The user
can learn more about these programs by typing at the terminal window

help "program name"

where "program name" is one of irp, mmp, etc.

These programs use the moremata package.

To install this package in Stata, type

"ssc install moremata"

Charles Lindsey
3/17/2009  