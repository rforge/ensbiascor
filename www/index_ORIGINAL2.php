
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<html>
<body>

<!-- R-Forge Logo -->


<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr>
<td>
<img src="logo.png" />
</td> 
</tr>
</table>

<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<h1>ensbiascoR: Ensemble Resampling Bias Correction</h1>
<em><p><strong>An R package to resample or constrain large model ensemble to preserve physical aspects of the model simulations.
</strong> </p></em>

<!-- menu -->
<hr>
<h2>Contents</h2>
<ul>
	<li><a href="index.html">ensbiascoR: Introduction and Installation</a></li> 
	<li>ensbiascor - applications and examples:</li>
	<ul>
		<li><a href="ensbiascoR_example1.html">Example 1: Ensemble resampling bias correction (HadRM3P and weather@home). Monthly temperature and precipitation</a></li>
		<li><a href="ensbiascoR_example2.html">Example 2: Ensemble resampling bias correction (HadRM3P and weather@home). Maxima of temperature and relative humidity</a></li>
		<li><a href="ensbiascoR_example3.html">Example 3: Land coupling constraints on multi-model ensembles (CMIP5)</a></li>
	</ul>
	<li><a href="http://r-forge.r-project.org/projects/ensbiascor/">Project summary page</a></li>
</ul>
<hr>
<!-- end of menu -->


<p> <h4>Introduction</h4>
Resampling initial condition or multi-model climate ensemble simulations offers a complementary methodology to traditional bias correction by retaining physical aspects of the original model simulations. The R-package "ensbiascor" offers tools and examples, including data, to do so.
</p>

<h4><p> Related papers: </h4>
Sippel, S., Otto, F. E. L., Forkel, M., Allen, M. R., Guillod, B. P., Heimann, M., Reichstein, M., Seneviratne, S. I., Thonicke, K. & Mahecha, M. D. (2016) A novel bias correction methodology for climate impact simulations. Earth System Dynamics, 7, 71-88. doi:10.5194/esd-7-71-2016.
</p>

<br>
<tr>
<td>
<!-- <embed src="illu_ensbiascor.pdf" height="500px" id="fig-1"/> -->
<embed src="illu_ensbiascor.png" id="fig-1"/>
</td>
</tr>
<tr>
<td>
<p style="font-family:arial;font-size:16px;"> <b>Figure 1.</b>  Illustration of ensemble-based resampling methodology. (a) Empirical cumulative density function of summer mean temperatures over Central Europe in ERA-Interim. The nonparametric fit to the cumulative density using a Gaussian kernel for observations and the model ensemble (HadRM3P, <a href="http://climateprediction.net/weatherathome/"> climateprediction.net/weatherathome</a>) are shown by the blue and red lines, respectively. (b) A transfer function between the observed and modelled distribution is derived using Cubic Hermite splines. (c) Quantile-quantile plot for the original and resampled ensemble for the JJA temperature constraint. (d) Fraction of original ensemble members in percentile bins of the observed distribution (blue line in (a)), i.e. "effective ensemble size" after resampling.</p>
</td>
<td>
<br>


<p> <h4> Installation and Usage </h4>
To install the most recent version of ensbiascor in R: <br><b>install.packages("ensbiascor", repos="http://R-Forge.R-project.org")</b><br>

<p>
If you use <em>ensbiascor</em> in scientific publications, please cite:
<br>
Sippel, S., Otto, F. E. L., Forkel, M., Allen, M. R., Guillod, B. P., Heimann, M., Reichstein, M., Seneviratne, S. I., Thonicke, K. & Mahecha, M. D. (2016) A novel bias correction methodology for climate impact simulations. Earth System Dynamics, 7, 71-88. doi:10.5194/esd-7-71-2016.</p>
The package has been developed at the Max Planck Institute for Biogeochemistry, Jena, Germany.

<p> <h4>Author details and further information:</h4>
 Sebastian Sippel (<a href="mailto:ssippel@bgc-jena.mpg.de">ssippel@bgc-jena.mpg.de</a>)
<br> For news and further information check my personal page: <a href="https://www.bgc-jena.mpg.de/bgi/index.php/People/SebastianSippel/">here</a>. </p> 
</body>
</html>
