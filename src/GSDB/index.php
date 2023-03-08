 <!DOCTYPE html>
<!-- Template by Quackit.com -->
<!-- Modified by oluwatosin oluwadare for Genome Structure Database -->
<?php
	include("php-wrapper/fusioncharts.php");
?>	
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <!-- The above 3 meta tags *must* come first in the head; any other head content must come *after* these tags -->

    <title>GSDB</title>

    <!-- Bootstrap 4 CSS. This is for the alpha 3 release of Bootstrap 4. This should be updated when Bootstrap 4 is officially released. -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-alpha.3/css/bootstrap.min.css" integrity="sha384-MIwDKRSSImVFAZCVLtU0LMDdON6KVCrZHyVQQj6e8wIEJkW4tvwqXrbMIya1vriY" crossorigin="anonymous">
           <link rel="icon" href="images/icon.PNG"> 
    <!-- Custom CSS: You can use this stylesheet to override any Bootstrap styles and/or apply your own styles -->
    <link href="css/custom.css" rel="stylesheet">
    
    <!-- For icons -->
    <link href="css/font-awesome-4.6.3/css/font-awesome.min.css" rel="stylesheet">

    <!---- CHARTS FROM DATABASE -->
	<script type="text/javascript" src="js/fusioncharts.js"></script>
	<script type="text/javascript" src="js/fusioncharts.charts.js"></script>
	<script type="text/javascript" src="js/themes/fusioncharts.theme.ocean.js"></script>
	
	    <!-- Navigation -->
    <nav id="topNav" class="navbar navbar-default navbar-static-top navbar-dark bg-inverse m-b-1">
        <button class="navbar-toggler hidden-md-up pull-right" type="button" data-toggle="collapse" data-target="#navbar">
            &#9776;
        </button>
        <a class="navbar-brand" href="index.php" >GSDB</a>
		<br>
        <div class="collapse navbar-toggleable-sm" id="navbar">
            <ul class="nav navbar-nav">
                <li class="nav-item">
                    <a class="nav-link" href="index.php">Home</a>
                </li>
               
			   <li class="nav-item">
                    <a class="nav-link" href="browse.php">Browse</a>
                </li>
				
				<li class="nav-item">
                    <a class="nav-link" href="softwaretools.html">Software Tools</a>
                </li>				
               
				<li class="nav-item">
                    <a class="nav-link" href="evaluate.php">Structure Evaluation</a>
                </li>
				
				<li class="nav-item">
                    <a class="nav-link" href="Tutorial.html">Tutorial</a>
                </li>
				 <li class="nav-item">
                    <a class="nav-link" href="contact.html">Contact Us</a>
                </li>
            </ul>
            <!-- Search 
            <form class="form-inline pull-xs-right">
                <input class="form-control" placeholder="GSDB ID" type="text">
                <button class="btn btn-secondary" type="submit">Search</button>
            </form>
			-->
            </ul>
        </div>
    </nav>

</head>

<body>
    

<div class="container-fluid">
		<!-- Center Column -->
	<div class="col-sm-12">
		<div class="row">
						
				<h1>GSDB : Genome Structure Database</h1>			
				<br/>
			
				<div style=" float: left;width: 400px;">
					  <div style="text-align: center;"><img src="images/genome_1.png" style='width:330px; ' alt="a genome visual" /></div><br/>	
				 </div >
					
				  <div style="float: none;overflow: hidden;padding: 20px"">
				        
						<p>GSDB is a database of three-dimensional (3D) chromosome and genome structures reconstructed from Hi-C data by multiple 3D modeling tools in the field.  In recent years, numerous Hi-C datasets have been generated and many chromosome/genome structure construction algorithms have been developed.
						 However, there is no public repository of 3D chromosome and genome structures available for the community to use.</p>
						<p>GSDB aims to fill this gap by providing a comprehensive repository of 3D chromosome and genome structures reconstructed from for Hi-C data. We hope this database will enable the exploration of the dynamic architectures of chromosomes and genomes for biomedical research.</p>
				
						<p> <i> Over 50,000 structures from 12 start-of-the-art Hi-C data structure prediction algorithms for 32 Hi-C datasets each containing varying resolutions. </i></p>
						<p><button  onclick="window.location='browse.php'" class="btn btn-outline-primary">Get Started</button></p>	
				  </div>
				  
				


							
		
		</div>
		<hr>
		<!-- Articles -->
		<div class="row">
		
			<?php include 'charts.php';?>
		  <div id="chart-2"class="col-sm-3" ><!-- Fusion Charts will render here--></div>
		  <div id="chart-3"class="col-sm-3" ><!-- Fusion Charts will render here--></div>
	      <div id="chart-4"class="col-sm-3" ><!-- Fusion Charts will render here--></div>
		  <div id="chart-1" class="col-sm-3"  ><!-- Fusion Charts will render here--></div>
		</div>	
	</div><!--/Center Column-->
		
</div><!--/container-fluid-->
	
	
<footer>

	<div class="small-print">
		<div >
		<center>
		<br/>
		<table>
		<tr>
		<td>
			<p>Copyright &copy; 2019 <a href="#"><a href="http://calla.rnet.missouri.edu/cheng/">BDM Lab</a> | <a href="mailto:chengji@missouri.edu">Contact</a></p>
			<img src="images/mu_resize.jpg" alt="University of Missouri">
			<img src="https://relay.hekademeia.org/track/image?src=hek-host-wc" alt="Hosting by Hekademeia" style="width: 160px;">
		</td>
		
		<td>
			</div><!-- BEGIN Free Counter Code! counter -->
			<script type="text/javascript" src="https://freecountercode.com/service/oyKAzHFID0hXdChGo3H5"></script>
			<!-- END Free Counter Code! counter -->
		</td>
		</tr>
		
		</table>
		</center>
		</div>
	</div>
</footer>

	
    <!-- Bootstrap core JavaScript
    ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->

    <!-- jQuery library -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.0.0/jquery.min.js" integrity="sha384-THPy051/pYDQGanwU6poAc/hOdQxjnOEXzbT+OuUAFqNqFjL+4IGLBgCJC3ZOShY" crossorigin="anonymous"></script>

    <!-- Tether -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/tether/1.2.0/js/tether.min.js" integrity="sha384-Plbmg8JY28KFelvJVai01l8WyZzrYWG825m+cZ0eDDS1f7d/js6ikvy1+X+guPIB" crossorigin="anonymous"></script>

    <!-- Bootstrap 4 JavaScript. This is for the alpha 3 release of Bootstrap 4. This should be updated when Bootstrap 4 is officially released. -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-alpha.3/js/bootstrap.min.js" integrity="sha384-ux8v3A6CPtOTqOzMKiuo3d/DomGaaClxFYdCu2HPMBEkf6x2xiDyJ7gkXU0MWwaD" crossorigin="anonymous"></script>

    <!-- Initialize Bootstrap functionality -->
    <script>
    // Initialize tooltip component
    $(function () {
      $('[data-toggle="tooltip"]').tooltip()
    })

    // Initialize popover component
    $(function () {
      $('[data-toggle="popover"]').popover()
    })
    </script> 
    
	<!-- Placeholder Images -->
	<!--<script src="js/holder.min.js"></script>-->
	
</body>

</html>
