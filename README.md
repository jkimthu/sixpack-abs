# sixpack-abs
Feeling sore? Keep going. Correlations between the many many parameters of growth.

last update: 2019 Feb 6th


* * used in manuscript one * *

# figure 2
50. volume vs time, selected individual tracks from all conditions			% 2A
7.  growth rate vs ti													% 2B, supp
23. single period: growth rate vs nutrient phase							% 2C

# figure 3
21. monod: growth rate vs nutrient conc										% 3A
26. pdf of growth rates normalized by pop average 							% 3B

# figure 4
22. growth rate vs timescale: measured vs expected						   	% 4		

# figure 5
51. plot response to upshifts												% 5A
		i. all replicates
	   ii. mean of replicates
	  iii. mean and standard deviation between replicate means

52. plot response to downshifts 											% 5A
		i. all replicates
	   ii. mean of replicates
	  iii. mean and standard deviation between replicate means

54. quantify time to stabilization and stabilized gr value for upshifts		% 5B
55. quanitfy time to stabilization and stabilized gr value for downshifts   % 5B





TABLE of CONTENTS:

1. division size vs birth size  					population
2. inter-division time vs birth size 				population
3. growth rate vs birth size						population

4. division size vs birth size, heated scatter		single				%

5. inter-division time vs birth size				single
6. division size vs birth size						single

7. growth rate vs time 								population
8. birth size vs growth rate						population
9. inter-division time vs growth rate				population
10. birth size vs growth rate						single
11. inter-division time vs growth rate  			single

12. volume vs nutrient phase						population
13. fraction births vs nutrient phase				single
14. individual volume vs time 						single
15. size (including SA:V) vs nutrient phase 		population
16. growth rate vs interdivision time 				population

17. added volume vs time 							population
18. added volume vs birth volume					single
19. added volume vs time 							single


21. monod: growth rate vs nutrient conc				population
22. growth rate: measured vs expected				population
23. single period: growth rate vs nutrient phase							% 2C


24. violin plots with dV/dt calculation performed BEFORE 3hr trim
25. subsampling growth rates in time bins

26. pdf of growth rates normalized by pop ave

27. timescale specific responses in growth rate to upshifts & downshifts

28. added volume vs birth volume, Taheri	

29. violin plots to compare distributions between fluc and stable
30. fraction negative growth rates vs time, shifts 
31. added size vs birth size (normalized by average birth size), trimmed x axis
32. like 31, but collapsing all curves

33. division size vs birth size - length & width
34. monod: growth rate vs nutrient conc - width 
35. monod: growth rate vs nutrient conc - length
36. growth rate: measured vs. expected in length and width
37. size at birth over time
38. added width vs. birth width, Taheri-like
39. interdivision time over time

40. pdf of interdivision times, raw

41. pdf of birth volumes times, raw

42. computing T-statistic, assuming noise from single cell variability 
43. computing T-statistic, assuming noise from day-to-day variability

46. length and volume vs time
47. width vs time
48. width, length and volume by cell cycle fraction
49. width vs time, selected individual tracks from all conditions
50. volume vs time, selected individual tracks from all conditions			% 2A

51. plot response to upshifts												% 5A
		i. all replicates
	   ii. mean of replicates
	  iii. mean and standard deviation between replicate means

52. plot response to downshifts 											% 5A
		i. all replicates
	   ii. mean of replicates
	  iii. mean and standard deviation between replicate means

53. plot mean response of all replicates to shifts								
54. quantify time to stabilization and stabilized gr value for upshifts			% 5B
55. quantify time to stabilization and stabilized gr value for downshifts 		% 5B

56. quantify features of fluctuating growth rate
	- time to stabilization
	- mean stabilized growth rate (per period)
	- mean growth rate in high & low nutrient (per period)
	- amplitude of fluctuations

	VERSION B. bins growth rates per hr instead of per period like in 56.m


57. quantify number of unique tracks per condition per experiment
	trackData.mat, has a cell per experiment with eight values
	each row represents a condition: fluc, low, ave, high
	each column represents total data (bubble trimmed) or post 3h data 


58. relationship between G_low, G_ave, G_high and G_fluc
	goal: determine strength of day-by-day correlations in growth rate
	outcome: non-negligible positive correlations. normalize by daily Gs!

	Supplementary Figure 7. in manuscript one (growth rate)


70. growth rate vs time for individual xys from a given condition, one expt
71. growth rate vs time for sensitivity analysis using 2018-02-01 				% supp

72. overlaid periods by increasing time - looking for adaptation				% 6C

73. Size distributions at single timepoints
    Of interest for Vicente