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
      
    <!-- Custom CSS: You can use this stylesheet to override any Bootstrap styles and/or apply your own styles -->
    <link href="css/custom.css" rel="stylesheet">
    
    <!-- For icons -->
    <link href="css/font-awesome-4.6.3/css/font-awesome.min.css" rel="stylesheet">
	
	<!--- CDN DataTable-->
	<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.css">
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/plug-ins/preview/searchPane/dataTables.searchPane.css">
	
    <!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
    <!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
    <!--[if lt IE 9]>
        <script src="https://oss.maxcdn.com/libs/html5shiv/3.7.0/html5shiv.js"></script>
        <script src="https://oss.maxcdn.com/libs/respond.js/1.4.2/respond.min.js"></script>
    <![endif]-->

</head>

<body>

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
                <li class="nav-item dropdown">
                    <a class="nav-link dropdown-toggle" href="#" id="dropdownMenuLink" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                    Data
                    </a>
                     <div class="dropdown-menu" aria-labelledby="dropdownMenuLink">
                        <a class="dropdown-item" href="search.php">Search</a>
                        <a class="dropdown-item" href="general.php">Summary</a>                        
                    </div>
                </li>               
                <li class="nav-item dropdown">
                    <a class="nav-link dropdown-toggle" href="#" id="dropdownMenuLink" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                    Materials and Method
                    </a>
                    <div class="dropdown-menu" aria-labelledby="dropdownMenuLink">
                        <a class="dropdown-item" href="publications.html">Publications</a>
                        <a class="dropdown-item" href="softwaretools.html">Software tools</a>                        
                    </div>
                </li>
				<li class="nav-item">
                    <a class="nav-link" href="viewer.php">Viewer</a>
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
                <input class="form-control" placeholder="GSDB ID" type="text" >
                <button class="btn btn-secondary" type="submit">Search</button>
            </form>
			-->
            </ul>
        </div>
    </nav>
    

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
						<th>No</th>
						<th>GSDB ID</th>
						<th>Title</th>
						<th>Biosample Type</th>
						<th>Organism</th>
						<th>Lab/Contributor(s)</th>					
						<th>Project</th>
						<th>Project ID</th>
						<th>GEO Accession ID</th>
					</tr>
				</thead>
				<tbody>
				<?php
					require_once('connection.php');
					$result=$conn->prepare("SELECT * FROM general_info ");
					$result->execute();
					for($i=0; $row = $result->fetch(); $i++){
				?>
					<tr>
						<td><?php echo $i+1; ?></td>
						<td><a href="details.php?id=<?php echo $row['ID']; ?>"><?php echo $row['ID']; ?></a></td>
						<td><?php echo $row['Title']; ?></td>
						<td><?php echo $row['Biosample_Type']; ?></td>
						<td><?php echo $row['Organism']; ?></td>
						<td><?php echo $row['Lab']; ?></td>						
						<td><?php echo $row['Project']; ?></td>
						<td><?php echo $row['Project_Id']; ?></td>
						<td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<?php echo $row['GEO_Accesion_No']; ?>"><?php echo $row['GEO_Accesion_No']; ?> </a></td>
						
					</tr>
					<?php } ?>
				</tbody>
			</table>		
			
		</div><!--/Center Column-->



	</div><!--/container-fluid-->
	
	<footer>
		        
        <div class="small-print">        	
        		<p><p>Copyright &copy; 2018 <a href="#"><a href="http://calla.rnet.missouri.edu/cheng/">BDM Lab</a> | <a href="mailto:chengji@missouri.edu">Contact</a></p>
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
					stateSave: true
				} );
				
			} );
					
	</script>
	

	   
	<!-- Placeholder Images -->
	<script src="js/holder.min.js"></script>
	
</body>

</html>
