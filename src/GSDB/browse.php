<?php
require_once(__DIR__."/connection.php");
?>
<!DOCTYPE html>
<!-- Template by Quackit.com -->
<!-- Modified by oluwatosin oluwadare for Genome Structure Database -->
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
    <link href="http://sysbio.rnet.missouri.edu/3dgenome/GSDB/css/custom.css" rel="stylesheet">
    
    <!-- For icons -->
    <link href="http://sysbio.rnet.missouri.edu/3dgenome/GSDB/css/font-awesome-4.6.3/css/font-awesome.min.css" rel="stylesheet">
	
	<!--- CDN DataTable-->
	 
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.css">
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/plug-ins/preview/searchPane/dataTables.searchPane.css">
	
	
    <!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
    <!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
    <!--[if lt IE 9]>
        <script src="https://oss.maxcdn.com/libs/html5shiv/3.7.0/html5shiv.js"></script>
        <script src="https://oss.maxcdn.com/libs/respond.js/1.4.2/respond.min.js"></script>
    <![endif]-->
    <!-- Navigation -->
    <nav id="topNav" class="navbar navbar-full navbar-static-top navbar-dark bg-inverse m-b-1">
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
               
               <!--
				<li class="nav-item">
                    <a class="nav-link" href="viewer.php">Viewer</a>
                </li>
				-->							
						
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
		
		<div id="table" class="col-sm-12">
		<h3 class="card-title">	<i class="fa fa-search" aria-hidden="true"></i> 
				Search Database 
				</h3>
				<hr>
			<table id="mytable1" class="display nowrap" cellspacing="0" frame="box" >
				<thead>
					<tr >						
						
						<th>Filename</th>
						<th>Hi-C Dataset Title</th>
						<th>3D Structure</th>
						<!--<th>Biosample Type</th> -->
						<th>Organism</th>
						<th>GSDB ID</th>											
						<th>Resolution</th>
						<th>Normalized Hi-C Data</th>
						<th>Project</th>
						<th>Project ID</th>
						<th>GEO Accession ID</th>
						
					</tr>
				</thead>
				
				<tbody>
				<?php
					// $query = "SELECT ID,Title,Biosample_Type,Organism,Project,Project_Id,GEO_Accesion_No  FROM general_info";
					$query = "SELECT ID,Title,Organism,Project,Project_Id,GEO_Accesion_No  FROM general_info";
					if ($stmt = $mysqli->prepare($query)) {
						/* execute query */
						$stmt->execute();
						/* store result */
						$stmt->store_result();				

					   /* Bind the result to variables */
					    // $stmt->bind_result($ID,$Title,$Biosample_Type,$Organism,$Project,$Project_Id,$GEO_Accesion_No);
						  $stmt->bind_result($ID,$Title,$Organism,$Project,$Project_Id,$GEO_Accesion_No);
						while ($stmt->fetch()) {
							$query2 = "SELECT Filename,Resolution FROM data_info WHERE GSDB_ID = ? ";
							if ($stmt2 = $mysqli->prepare($query2)){								
								$stmt2->bind_param("s", $ID);								
								/* execute query */
								$stmt2->execute();								
								/* Store the result (to get properties) */
								$stmt2->store_result();
							   /* Get the number of rows */
							   $num_of_rows2 = $stmt2->num_rows;
								/* Bind the result to variables */
							  $stmt2->bind_result($Filename,$Resolution);
							  while ($stmt2->fetch()) {
				?>
							<tr>	
									<td><a href="details.php?id=<?php echo $Filename; ?>"><?php echo $Filename; ?></a></td>										
									<td><?php echo $Title; ?></td>
									<?php 							
										$dir_download = 'http://calla.rnet.missouri.edu/genome3d/GSDB/Database/' . $ID . "/" . $Filename. "_3DStructures.tar.gz";										 
									?>
									<td><a href="details.php?id=<?php echo $Filename; ?>">View</a> | <a href="<?php echo " "; echo $dir_download; ?>">Download </a></td>			
									<!--<td><?php echo $Biosample_Type; ?></td>-->
									<td><?php echo $Organism; ?></td>
									<td><a href="details.php?id=<?php echo $ID; ?>"><?php echo $ID; ?></a></td>		
									<td><?php echo $Resolution; ?></td>
									
									<?php 							
										$data_download = 'http://calla.rnet.missouri.edu/genome3d/GSDB/Dataset/' . $ID . "/" . $Filename. ".tar.gz";		
										
									?>
									<td><center><a href="<?php echo " "; echo $data_download; ?>">Download </a></center></td>
									<td><?php echo $Project; ?></td>
									<td><?php echo $Project_Id; ?></td>
									<td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<?php echo $GEO_Accesion_No; ?>"><?php echo $GEO_Accesion_No; ?> </a></td>
							</tr>		
						
				<?php	
							  }										
							}		
										   
						}							
							/* free result */
							$stmt->free_result();
							/* close statement */
							$stmt->close();
					}

						/* close connection */
						$mysqli->close();
				?>
				</tbody>
			</table>		
			
		</div><!--/Center Column-->


	  

	</div><!--/container-fluid-->
	
	<footer>
		        
        <div class="small-print">        	
        		<p><p>Copyright &copy; 2019 <a href="#"><a href="http://calla.rnet.missouri.edu/cheng/">BDM Lab</a> | <a href="mailto:chengji@missouri.edu">Contact</a></p>
        		<img src="images/mu_resize.jpg" alt="">
			
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

	<!----CDN DataTable -->
	<script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.js"></script>
	<script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/plug-ins/preview/searchPane/dataTables.searchPane.js"></script>

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
   
   <!----CDN DataTable -->
	<script type="text/javascript" charset="utf-8">
			$(document).ready( function () {
			    $('#mytable1').DataTable( {
					searchPane: true,
					stateSave: false,
					 "order": [[ 4, "asc" ]]					
				} );
				
			} );
			

			
	</script>
	

	   
	<!-- Placeholder Images -->
	<script src="js/holder.min.js"></script>
	
</body>

</html>
