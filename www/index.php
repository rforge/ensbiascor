<!DOCTYPE html>
<!-- saved from url=(0014)about:internet -->
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
<meta http-equiv="x-ua-compatible" content="IE=9" >

<title>ensbiascoR: Ensemble Resampling Bias Correction</title>

<style type="text/css">
body, td {
   font-family: sans-serif;
   background-color: white;
   font-size: 12px;
   margin: 8px;
}

tt, code, pre {
   font-family: 'DejaVu Sans Mono', 'Droid Sans Mono', 'Lucida Console', Consolas, Monaco, monospace;
}

h1 { 
   font-size:2.2em; 
}

h2 { 
   font-size:1.8em; 
}

h3 { 
   font-size:1.4em; 
}

h4 { 
   font-size:1.0em; 
}

h5 { 
   font-size:0.9em; 
}

h6 { 
   font-size:0.8em; 
}

a:visited {
   color: rgb(50%, 0%, 50%);
}

pre {	
   margin-top: 0;
   max-width: 95%;
   border: 1px solid #ccc;
   white-space: pre-wrap;
}

pre code {
   display: block; padding: 0.5em;
}

code.r, code.cpp {
   background-color: #F8F8F8;
}

table, td, th {
  border: none;
}

blockquote {
   color:#666666;
   margin:0;
   padding-left: 1em;
   border-left: 0.5em #EEE solid;
}

hr {
   height: 0px;
   border-bottom: none;
   border-top-width: thin;
   border-top-style: dotted;
   border-top-color: #999999;
}

@media print {
   * { 
      background: transparent !important; 
      color: black !important; 
      filter:none !important; 
      -ms-filter: none !important; 
   }

   body { 
      font-size:12pt; 
      max-width:100%; 
   }
       
   a, a:visited { 
      text-decoration: underline; 
   }

   hr { 
      visibility: hidden;
      page-break-before: always;
   }

   pre, blockquote { 
      padding-right: 1em; 
      page-break-inside: avoid; 
   }

   tr, img { 
      page-break-inside: avoid; 
   }

   img { 
      max-width: 100% !important; 
   }

   @page :left { 
      margin: 15mm 20mm 15mm 10mm; 
   }
     
   @page :right { 
      margin: 15mm 10mm 15mm 20mm; 
   }

   p, h2, h3 { 
      orphans: 3; widows: 3; 
   }

   h2, h3 { 
      page-break-after: avoid; 
   }
}

</style>





</head>

<body>
<!-- This is the project specific website template -->

<!-- It can be changed as liked or replaced by other content -->

<p><html>
<body></p>

<!-- R-Forge Logo -->

<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr>
<td>
<img src="data:image/png;base64,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" />
</td> 
</tr>
</table>

<!-- get project title  -->

<!-- own website starts here, the following may be changed as you like -->

<h1>ensbiascoR: Ensemble Resampling Bias Correction</h1>

<p><em><p><strong>An R package to resample or constrain large model ensemble to preserve physical aspects of the model simulations.
</strong> </p></em></p>

<!-- menu -->

<hr>

<h2>Contents</h2>

<ul>
  <li><a href="index.html">ensbiascoR: Introduction and Installation</a></li> 
    <li>ensbiascoR - applications and examples:</li>
    <ul>
        <li><a href="ensbiascoR_example1.html">Example 1: Ensemble resampling bias correction (HadRM3P and weather@home). Monthly temperature and precipitation</a></li>
        <li><a href="ensbiascoR_example2.html">Example 2: Ensemble resampling bias correction (HadRM3P and weather@home). Seasonal maxima of temperature and relative humidity</a></li>
        <li><a href="ensbiascoR_example3.html">Example 3: Land coupling constraints on multi-model ensembles (CMIP5)</a></li>
    </ul>
    <li><a href="http://r-forge.r-project.org/projects/ensbiascor/">Project summary page</a></li>
</ul>

<hr>

<!-- end of menu -->

<p> <h4>Introduction</h4>
Resampling initial condition or multi-model climate ensemble simulations offers a complementary methodology to traditional bias correction by retaining physical aspects of the original model simulations. The R-package &ldquo;ensbiascoR&rdquo; offers tools and examples, including data, to do so.
</p>

<h4><p> Related papers: </h4>

<p>Sippel, S., Otto, F. E. L., Forkel, M., Allen, M. R., Guillod, B. P., Heimann, M., Reichstein, M., Seneviratne, S. I., Thonicke, K. &amp; Mahecha, M. D. (2016) A novel bias correction methodology for climate impact simulations. Earth System Dynamics, 7, 71-88. doi:10.5194/esd-7-71-2016.
</p></p>

<p><img src="illu_ensbiascor.pdf" alt=""/></p>

<p> <h4> Installation and Usage </h4>
To install the most recent version of ensbiascor in R: <br><b>install.packages(&ldquo;ensbiascor&rdquo;, repos=&ldquo;http://R-Forge.R-project.org&rdquo;)</b><br>

<p>
If you use <em>ensbiascor</em> in scientific publications, please cite:
<br>
Sippel, S., Otto, F. E. L., Forkel, M., Allen, M. R., Guillod, B. P., Heimann, M., Reichstein, M., Seneviratne, S. I., Thonicke, K. & Mahecha, M. D. (2016) A novel bias correction methodology for climate impact simulations. Earth System Dynamics, 7, 71-88. doi:10.5194/esd-7-71-2016.</p>

<p>The package has been developed at the Max Planck Institute for Biogeochemistry, Jena, Germany.</p>

<p> <h4>Author details and further information:</h4>
 Sebastian Sippel (<a href="mailto:ssippel@bgc-jena.mpg.de">ssippel@bgc-jena.mpg.de</a>)
<br> For news and further information check my personal page: <a href="https://www.bgc-jena.mpg.de/bgi/index.php/People/SebastianSippel/">here</a>. </p> 

<p></body>
</html></p>

</body>

</html>

