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
	
	<!--- 3Dmol-->
	<script type="" src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script> 
	<script src="http://cdn.jsdelivr.net/3dmol.js/latest/3Dmol-min.js"></script>
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
               
			   <li class="nav-item">
                    <a class="nav-link" href="search.php">Search</a>
                </li>
				
				<li class="nav-item">
                    <a class="nav-link" href="softwaretools.html">Software Tools</a>
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
            <!-- Search -->
            <form class="form-inline pull-xs-right">
                <input class="form-control" placeholder="GSDB ID" type="text">
                <button class="btn btn-secondary" type="submit">Search</button>
            </form>

            </ul>
        </div>
    </nav>
    
    

<div class="container-fluid">	
	<!-- Left Column -->
		<div class="col-sm-4">			
		
			<div class="row">
				<article class="col-xs-12">
				<h3 class="card-title">	<i class="fa fa-tachometer" aria-hidden="true"></i> 
				Download Structures
				</h3>	
			<br>			
			<!-- Table 1 -->	
			<table id="mytable1" class="display  nowrap" cellspacing="0" frame="box" >
					<thead>
						<tr >						
							<th>Filename</th>
							<th>Resolution</th>
							<th>GSDB ID</th>		
							
						</tr>
					</thead>
					<tbody>
	<?php
				
				
					$mysqli = new mysqli("localhost", "root" , "genomeflow", "gsdb");

					/* check connection */
					if (mysqli_connect_errno()) {
						printf("Connect failed: %s\n", mysqli_connect_error());
						exit();
					}

					$query = "SELECT ID FROM general_info";
					if ($stmt = $mysqli->prepare($query)) {
						/* execute query */
						$stmt->execute();
						/* store result */
						$stmt->store_result();				

					   /* Bind the result to variables */
					   $stmt->bind_result($ID);
						while ($stmt->fetch()) {
							$query2 = "SELECT Filename,Resolution FROM data_info WHERE GSDB_ID = ?";
							if ($stmt2 = $mysqli->prepare($query2)){								
								$stmt2->bind_param("s", $ID);								
								/* execute query */
								$stmt2->execute();								
								/* Store the result (to get properties) */
								$stmt2->store_result();
							   /* Get the number of rows */
							   $num_of_rows2 = $stmt2->num_rows;
								/* Bind the result to variables */
							  $stmt2->bind_result($Filename,$Resoluiton);
							  while ($stmt2->fetch()) {
				?>
							<tr>
									<td><?php echo $Filename; ?></td>						
									<td><?php echo $Resoluiton; ?></td>
									<td><?php echo $ID; ?></td>	
									
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
				</article>
			</div>
			<hr>	
			<!-- Alert -->
			<div class="alert alert-success alert-dismissible" role="alert">
				<button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				 Click on the Filename to <strong>download</strong> constructed structures for multiple algorithms
			</div>
		</div>
		
	<!--/Left Column-->
	<!-- center Column -->
	  <div class="col-sm-1">
						

	  </div>
	  <!--/center Column -->	
		
		
	<!-- Right Column -->
	 <div class="col-sm-7">
	  
		<h3 class="card-title">	<i class="fa fa-eye" aria-hidden="true"></i> 
				Visualization								
		</h3>
		<p>Select the <b>Algorithm</b> the  <b>Dataset </b> and the  <b>chromosome</b> and Click on the <b>Display</b>  to view the sctructure on the viewer.<p> 
				
		<table cellpadding="10"> 
		  <tr> 
		    <td>
		
				<h6>Select Algorithm</h6>	
				<div class="dropdown dropdown-dark">
					<select id="alg" class="dropdown-select">				 
					  <option value="LorDG">LorDG</option>
					  <option value="3DMax">3DMax</option>
					  <option value="MOGEN">MOGEN</option>
					</select>
				 </div>
							 
				 
							 
				<h6>Select dataset</h6>
				<div class="dropdown dropdown-dark">
					<select id="filename" class="dropdown-select">
					<?php
						  require_once('connection.php');
						$result=$conn->prepare("SELECT * FROM general_info ");
						$result->execute();
						for($i=0; $row = $result->fetch(); $i++){						
							$ID=$row['ID'];
							$local_result=$conn->prepare("SELECT * FROM data_info WHERE GSDB_ID = '$ID'");
							$local_result->execute();
							for($j=0; $inner_row = $local_result->fetch(); $j++){	
								$filevalue=$inner_row['Filename'];
								$value= $ID."-".$filevalue;
								$Resvalue = $inner_row['Resolution'];
								$Value= $value."-".$Resvalue;
							?>						
								
								<option value= <?php echo $Value; ?>> <?php echo $value; ?> </option>;
						<?php } } ?>
						</select>				
				 </div>
				 <br>
				 
				 <h6>Select Resolution</h6>
				<div class="dropdown dropdown-dark">
					<select id="Resolution" class="dropdown-select">
					
					</select>				
				 </div>
				 <br>
				 
				 
					
				<h6>Select Structure file</h6>	
				<div class="dropdown dropdown-dark">
					<select id="chr" class="dropdown-select">				 
					  <option value="chr1">chr1</option>
					  <option value="chr2">chr2</option>
					  <option value="chr3">chr3</option>
					  <option value="chr4">chr4</option>
					  <option value="chr5">chr5</option>
					  <option value="chr6">chr6</option>
					  <option value="chr7">chr7</option>
					  <option value="chr8">chr8</option>
					  <option value="chr9">chr9</option>
					  <option value="chr10">chr10</option>
					  <option value="chr11">chr11</option>
					  <option value="chr12">chr12</option>
					  <option value="chr13">chr13</option>
					  <option value="chr14">chr14</option>
					  <option value="chr15">chr15</option>
					  <option value="chr16">chr16</option>
					  <option value="chr17">chr17</option>
					  <option value="chr18">chr18</option>
					  <option value="chr19">chr19</option>
					  <option value="chr20">chr20</option>
					  <option value="chr21">chr21</option>
					  <option value="chr22">chr22</option>
					  <option value="chrX">chrX</option>
					  <option value="chrY">chrY</option>
					   <option value="Genome">All</option>
					</select>
				</div>
				  
							
				<br><br>
				<p><input type="submit" class="btn btn-outline-primary"	 value="Display" name="submit" onclick="myfunc()" /></p> 		
			</td>
		    <td>
			<!-- Inserted the Visualization here -->
		     <iframe width=700, height=625, src="jsmol.html"></iframe> 		  
		    </td>
		</tr> 
		<!--- Comment for Evaluarion Result-->	
		<tr>
			 <td> <b>Evaluation Result :</b>	</td>				
			 <td>  <output id="out">	</td>
		 </tr>
		
			
		</table>
		
		
	</div><!--/Right Column -->	
		
	 
</div><!--/container-fluid-->
	
	<footer>
		        
        <div class="small-print">
        	<div class="container">
        		<p><p>Copyright &copy; 2018 <a href="#"><a href="http://calla.rnet.missouri.edu/cheng/">BDM Lab</a> | <a href="mailto:chengji@missouri.edu">Contact</a></p>
        		<img src="images/mu_resize.jpg" alt="">
									
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
				 $('#mytable1').DataTable({	
					"searching": false,
					
				} );
			} );
					
	</script>
	
	
	<!----For Visualization -->
	<script type="text/javascript" charset="utf-8">
		var prev_name="";
		function myfunc() {
				//algorithm
				var alg_e = document.getElementById ("alg");
				var algo = alg_e.options [alg_e.selectedIndex] .value;
				 //filename
				var file_e = document.getElementById ("filename");
				var fname = file_e.options [file_e.selectedIndex] .value;	
				//get current filename				
				if (fname!=prev_name){
					 document.getElementById("Resolution").innerHTML = "";
					 prev_name=fname;					 
				}
				 //chromosome
				var chr_e = document.getElementById ("chr");
				var chrpdb = chr_e.options [chr_e.selectedIndex] .value;	
				
				var base_link="http://calla.rnet.missouri.edu/genome3d/GSDB/Database/";
				
				//split filename to get ID
				var splitString = fname .split("-");
				var prefix = splitString[0];	// get id
				var gsdb_id = prefix;	
				var fname = splitString[1];		 //get filename	
				
				var Res = splitString[2];		//get Resolution
				var splitRes = Res.split(",");
				
				//Resolution
				 var ResList = document.getElementById('Resolution');	
			   // if empty no optiond create one
				if( !$('#Resolution').has('option').length > 0 ) {
					
					for (var i = 0; i < splitRes.length; i++) {
						ResList.options[i] = new Option(splitRes[i]);					
					}									
					var res = ResList.options[0] .value;						
				} else{
					// use the selected value
					
					var res = ResList.options [ResList.selectedIndex] .value;	
				}
			   														 
				var link = "";
				norm_VC ="VC";
				norm_YT ="Yaffe_Tanay";
				norm="";
				if (gsdb_id.trim()=="TE1402WS" || gsdb_id.trim()=="VH2561BL" ||gsdb_id.trim()=="BB8015WF" ||gsdb_id.trim()=="GG6098MH"){
					norm =norm_YT;
				}else{
					norm = norm_VC;
				}
				
			   if (splitRes.length == 1){
				   link = base_link + gsdb_id + "/" + fname +  "/" + norm + "/" + algo ;
			   }
			   else{
				     var res = res.toLowerCase();
				     link = base_link + gsdb_id + "/" + fname +  "/" + norm  +  "_" + res + "/" + algo ;
			   }
						
			alert(link);
			// document.write(link);	
				
				
		}
    </script>
	
	<!-- Placeholder Images -->
	<script src="http://sysbio.rnet.missouri.edu/3dgenome/GSDB/js/holder.min.js"></script>
	
</body>

</html>
