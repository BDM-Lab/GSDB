
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
    <link href="css/custom.css" rel="stylesheet">
    
    <!-- For icons -->
    <link href="css/font-awesome-4.6.3/css/font-awesome.min.css" rel="stylesheet">
	
	<!--- CDN DataTable-->
	<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.css">
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/plug-ins/preview/searchPane/dataTables.searchPane.css">
	
	<!--- 3Dmol-->
	<script type="" src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script> 
	

    <!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
    <!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
    <!--[if lt IE 9]>
        <script src="https://oss.maxcdn.com/libs/html5shiv/3.7.0/html5shiv.js"></script>
        <script src="https://oss.maxcdn.com/libs/respond.js/1.4.2/respond.min.js"></script>
    <![endif]-->
	
	<!-- include alertify script -->
	<script src="alertifyjs/alertify.js"></script>

	<!-- include alertify.css -->
	<link rel="stylesheet" href="alertifyjs/css/alertify.css">
	
	<!-- include semantic ui theme  -->
	<link rel="stylesheet" href="alertifyjs/css/themes/semantic.css">
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
		<h3 class="card-title">	<i class="fa fa-eye" aria-hidden="true"></i> 
				Side by Side Structure Visualization			
		</h3>
		<p>Select the <b>Algorithm</b>,<b> Dataset </b>, <b> Resolution </b>, <b> Chromosome</b>, and Click on the <b>Display</b> button to view the structure on the viewer.<p> 
		<hr>

	
	<h6><a href="details.php?id=<?=urlencode($_GET['id']); ?>"><font size="5"> &larr; </font>  <u> Back to single structure view page </u></a></h6>	
<table cellspacing="10" cellpadding="10">
<tr>
<td>
	<!-- Left Column -->	
			
		
				<h6>Select Algorithm</h6>	
				<div class="dropdown dropdown-dark">
					<select id="alg" class="dropdown-select">				 
					  <option value="LorDG">LorDG</option>
					  <option value="3DMax">3DMax</option>
					  <option value="MOGEN">MOGEN</option>
					  <option value="Pastis">Pastis</option>
					  <option value="Chromosome3D">Chromosome3D</option>
					  <option value="HSA">HSA</option>
					  <option value="miniMDS">miniMDS</option>
					  <option value="ShRec3D">ShRec3D</option>
					  <option value="GEM">GEM</option>
					  <option value="ChromSDE">ChromSDE</option>
					  <option value="SIMBA3D">SIMBA3D</option>					  
					 <option value="InfMod3DGen">InfMod3DGen</option> 
					 
				  
					</select>
				 </div>
				 
				<h6>Select dataset</h6>
				<div class="dropdown dropdown-dark">
					<select id="filename" class="dropdown-select">
					<?php
						  $var_value = $_GET['id'];	
							require_once('connection.php');
							$result=$mysqli->prepare("SELECT GSDB_ID,Filename,Resolution FROM data_info WHERE GSDB_ID = ? OR Filename = ?");
							$result->bind_param("ss",$var_value,$var_value);
							$result->execute();
							
							$result->bind_result($var_value,$filevalue,$Resvalue);
							for($i=0; $result->fetch(); $i++){
								$Value= $filevalue."-".$Resvalue;
								?>
								<option value= <?php echo $Value; ?>> <?php echo $filevalue; ?> </option>;
						<?php } ?>
						</select>				
				 </div>
				 <br>
				 
				 
				 <h6>Select Resolution</h6>
				<div class="dropdown dropdown-dark">
					<select id="Resolution" class="dropdown-select">
					
					</select>				
				 </div>
				 <br>
					
				<h6>Select Chromosome </h6>	
				<div class="dropdown dropdown-dark">
					<select id="chr" class="dropdown-select">				 
					  <option value="1">chr1</option>
					  <option value="2">chr2</option>
					  <option value="3">chr3</option>
					  <option value="4">chr4</option>
					  <option value="5">chr5</option>
					  <option value="6">chr6</option>
					  <option value="7">chr7</option>
					  <option value="8">chr8</option>
					  <option value="9">chr9</option>
					  <option value="10">chr10</option>
					  <option value="11">chr11</option>
					  <option value="12">chr12</option>
					  <option value="13">chr13</option>
					  <option value="14">chr14</option>
					  <option value="15">chr15</option>
					  <option value="16">chr16</option>
					  <option value="17">chr17</option>
					  <option value="18">chr18</option>
					  <option value="19">chr19</option>
					  <option value="20">chr20</option>
					  <option value="21">chr21</option>
					  <option value="22">chr22</option>
					  <option value="X">chrX</option>
					  <option value="Y">chrY</option>					
					</select>
				</div>
				  
							
				<br><br>
				<p><input type="submit" class="btn btn-outline-primary"	 value="Display" name="submit" onclick="view_struct1()" /></p> 		
				
			</td>
		    <td>
			<!-- Inserted the Visualization here -->
		       <iframe id="StructFrame" width=600, height=500, src="https://3dmol.csb.pitt.edu/viewer.html?url=https://gsdb.mu.hekademeia.org/structures/IFList_Chr_20_1mb_1450749845157.pdb&type=pdb&select=all&style=stick:color~crimson,radius~0.05" ></iframe> 		 
				<br>
				<b>Structure Evaluation: </b> 	<br/>	
				<iframe id="StructLog" width=600, height=105, src="" ></iframe> 
		
	
	<!--/Left Column-->
</td>

<td>
	<!-- center Column -->
	 
	  <!--/center Column -->	
</td>


<td>
	<!-- Right Column -->	 
		<table cellpadding="10"> 
		
		  <tr> 
		   <!--- Column 2 --->
		    <td>
		
				<h6>Select Algorithm</h6>	
				<div class="dropdown dropdown-dark">
					<select id="alg1" class="dropdown-select">				 
					  <option value="LorDG">LorDG</option>
					  <option value="3DMax">3DMax</option>
					  <option value="MOGEN">MOGEN</option>
					  <option value="Pastis">Pastis</option>
					  <option value="Chromosome3D">Chromosome3D</option>
					  <option value="HSA">HSA</option>
					  <option value="miniMDS">miniMDS</option>
					  <option value="ShRec3D">ShRec3D</option>
					  <option value="GEM">GEM</option>
					  <option value="ChromSDE">ChromSDE</option>
					  <option value="SIMBA3D">SIMBA3D</option>					  
					 <option value="InfMod3DGen">InfMod3DGen</option> 
					 
				  
					</select>
				 </div>
				 
				<h6>Select dataset</h6>
				<div class="dropdown dropdown-dark">
					<select id="filename1" class="dropdown-select">
					<?php
						  $var_value = $_GET['id'];	
							require_once('connection.php');
							$result=$mysqli->prepare("SELECT GSDB_ID,Filename,Resolution FROM data_info WHERE GSDB_ID = ? OR Filename = ?");
							$result->bind_param("ss",$var_value,$var_value);
							$result->execute();
							
							$result->bind_result($var_value,$filevalue,$Resvalue);
							for($i=0; $result->fetch(); $i++){
								$Value= $filevalue."-".$Resvalue;
								?>
								<option value= <?php echo $Value; ?>> <?php echo $filevalue; ?> </option>;
						<?php } ?>
						</select>				
				 </div>
				 <br>
				 
				 
				 <h6>Select Resolution</h6>
				<div class="dropdown dropdown-dark">
					<select id="Resolution1" class="dropdown-select">
					
					</select>				
				 </div>
				 <br>
					
				<h6>Select Chromosome </h6>	
				<div class="dropdown dropdown-dark">
					<select id="chr1" class="dropdown-select">				 
					  <option value="1">chr1</option>
					  <option value="2">chr2</option>
					  <option value="3">chr3</option>
					  <option value="4">chr4</option>
					  <option value="5">chr5</option>
					  <option value="6">chr6</option>
					  <option value="7">chr7</option>
					  <option value="8">chr8</option>
					  <option value="9">chr9</option>
					  <option value="10">chr10</option>
					  <option value="11">chr11</option>
					  <option value="12">chr12</option>
					  <option value="13">chr13</option>
					  <option value="14">chr14</option>
					  <option value="15">chr15</option>
					  <option value="16">chr16</option>
					  <option value="17">chr17</option>
					  <option value="18">chr18</option>
					  <option value="19">chr19</option>
					  <option value="20">chr20</option>
					  <option value="21">chr21</option>
					  <option value="22">chr22</option>
					  <option value="X">chrX</option>
					  <option value="Y">chrY</option>					
					</select>
				</div>
				  
							
				<br><br>
				<p><input type="submit" class="btn btn-outline-primary"	 value="Display" name="submit1" onclick="view_struct2();" /></p> 		
				
			</td>
		    <td>
			<!-- Inserted the Visualization here -->
		       <iframe id="StructFrame1" width=600, height=500, src="https://3dmol.csb.pitt.edu/viewer.html?url=https://gsdb.mu.hekademeia.org/structures/IFList_Chr_20_1mb_1450749845157.pdb&type=pdb&select=all&style=stick:color~crimson,radius~0.05" ></iframe> 		 
				<br>
				<b>Structure Evaluation: </b> 	<br/>	
				<iframe id="StructLog1" width=600, height=105, src="" ></iframe> 
			
				</td>
				
				
				
		</tr> 		
		</table>

</td>

</tr>
</table>	


	 
</div><!--/container-fluid-->
	
	<center>
		<br/>		        
        <div class="small-print">
        	<div >
        		<p><p>Copyright &copy; 2018 <a href="#"><a href="http://calla.rnet.missouri.edu/cheng/">BDM Lab</a> | <a href="mailto:chengji@missouri.edu">Contact</a></p>
        		<img src="images/mu_resize.jpg" alt="">
						<img src="https://relay.hekademeia.org/track/image?src=hek-host-wc" alt="Hosting by Hekademeia" style="width: 160px;">					
        	</div>
        </div>
	</center>


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
	
	
<script type="text/javascript">
			$(document).ready(function(){	
				var zip = new JSZip();
				zip.file("Hello.txt", "Hello world\n");

				jQuery("#data_uri").on("click", function () {
					zip.generateAsync({type:"base64"}).then(function (base64) {
						window.location = "data:application/zip;base64," + base64;
					}, function (err) {
						jQuery("#data_uri").text(err);
					});
				});
			})
	</script>
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
					"searching": false
				} );
			} );
			
			//Refresh page on load
			$( document ).ready(function() {
				view_struct1();
				view_struct2();
			});
		
	</script>
	
	
	<!----For Visualization -->
	<script type="text/javascript" charset="utf-8">
	
		function view_struct1(){
			//algorithm
			var alg_e = document.getElementById ("alg");
			
			//filename
			var file_e = document.getElementById ("filename");
			//chromosome
			var chr_e = document.getElementById ("chr");
			
			//Resolution
		
			var ResList = document.getElementById('Resolution');	
			var ResName = 'Resolution';
			//call function
			var result = new Array();
			result = myfunc(alg_e, file_e , chr_e, ResList,ResName );
			Struct_display = result[0];
			Struct_Log = result[1];
			// Display structure and log
			document.getElementById('StructFrame').src = Struct_display;	
			document.getElementById('StructLog').src = Struct_Log;				
		}	
		function view_struct2(){
			//algorithm
			var alg_e = document.getElementById ("alg1");
			
			//filename
			var file_e = document.getElementById ("filename1");
			//chromosome
			var chr_e = document.getElementById ("chr1");
			
			//Resolution
		
			var ResList = document.getElementById('Resolution1');	
			var ResName = 'Resolution1';
			//call function
			var result = new Array();
			result = myfunc(alg_e, file_e , chr_e, ResList,ResName );
			Struct_display = result[0];
			Struct_Log = result[1];
			// Display structure and log
			document.getElementById('StructFrame1').src = Struct_display;	
			document.getElementById('StructLog1').src = Struct_Log;				
		}	

	</script>
		
	<script type="text/javascript" charset="utf-8">
		var prev_name="";
		function myfunc(alg_e, file_e , chr_e, ResList,ResName ) {	
		        
				//algorithm
				var algo = alg_e.options [alg_e.selectedIndex] .value;
				 //filename			
				var fname = file_e.options [file_e.selectedIndex] .value;		
				
				//get current filename				
				if (fname!=prev_name){
					 document.getElementById(ResName).innerHTML = "";
					 prev_name=fname;					 
				}
				
				//chromosome			
				var chrp = chr_e.options [chr_e.selectedIndex] .value;
				var chrpdb = "chr" + chrp;
				
				
				var base_link="http://calla.rnet.missouri.edu/genome3d/GSDB/Database/";
				
				var gsdb_id = "<?php echo $var_value ?>"; 
				
				//split filename to get name and Resolution
				var splitString = fname .split("-");
				
				var fname = splitString[0];	// get id				
				var Res = splitString[1];	//get filename
				var splitRes = Res.split(",");		//Resolution
		
				// if empty no optiond create one
				var hname= '#' + ResName;
				
				if( !$(hname).has('option').length > 0 ) {
					
					for (var i = 0; i < splitRes.length; i++) {
						ResList.options[i] = new Option(splitRes[i]);				
					}		
					
					var len =splitRes.length-1;
					var res = ResList.options[len] .value;
					document.getElementById(ResName).value = res; //set selected vale
					
				} else{
					// use the selected value
					
					var res = ResList.options [ResList.selectedIndex] .value;	
				}
				
				
				
			 window.showAlert1 = function(){
					alertify.alert('3D structure not available', 'The selected Resolution is too high for Algorithm. ' + '<a href="resolution_message.html" target="_blank">Click for more details</a>');
			}
			
			 window.showAlert2 = function(){
					alertify.alert('3D structure not available', 'Due to time and memory constraint, only chromosome structure 10 to 23 could be constructed at this resolution. \t'  + '<a href="resolution_message.html" target="_blank">Click for more details</a>');
				}
				
			  window.showAlert3 = function(){
					alertify.alert('3D structure not available', 'Mouse cell( Mouse ES cell, Mouse Cortex) has Only Chromosome 1 - 20');
				}
			
			 window.showAlert4 = function(){
					alertify.alert('3D structure not available', 'Only Data for Chromosome 1 - 23 is available for the Hi-C data. Note: X = Chr. 23, Y = Chr. 24  ');
				}
			 window.showAlert5 = function(){
					alertify.alert('3D structure under construction', '3D Structure for Hi-C dataset will be uploaded once it is available.');
				}	
				
			 window.showAlert6 = function(){
					alertify.alert('3D structure not available', 'Genome 3D Structure for this Hi-C dataset is available for LorDG, MOGEN, ShRec3D, and miniMDS tools only.');
				}
			 
				//chromsome name for chromosome 23 = GM,K562,hESC,hIMR90
				if((gsdb_id.trim()=="BB8015WF" || gsdb_id.trim()=="GG6098MH" || gsdb_id.trim()=="OO7429SF"|| gsdb_id.trim()=="EA2504YQ"|| gsdb_id.trim()=="WT9059TG") && (chrpdb =="chrX")){				
					chrpdb ="chr23";					
				}
				if((gsdb_id.trim()=="BB8015WF" || gsdb_id.trim()=="GG6098MH" || gsdb_id.trim()=="OO7429SF"|| gsdb_id.trim()=="EA2504YQ"|| gsdb_id.trim()=="WT9059TG") && (chrpdb =="chrY")){				
					chrpdb ="chr24";					
				}
				if(fname.trim()=="GM06990_Genome" ){
					chrpdb="All";	
				}
				

				
			   // Highest Resolution for Pastis, chromosome3D, chromSDE is 250KB			   
			   
				if ((algo.trim()=="Pastis" || algo.trim()=="Chromosome3D" || algo.trim()=="ChromSDE"|| algo.trim()=="GEM"|| algo.trim()=="HSA"|| algo.trim()=="ChromSDE"|| algo.trim()=="ShRec3D" || algo.trim()=="InfMod3DGen" ) &&  (res.trim()=="50KB" || res.trim()=="40KB"  || res.trim()=="25KB" )){
					//alert ("The selected Resolution is too high for Algorithm. No structure available ");					
					 window.showAlert1();
				}
				else if ((algo.trim()=="SIMBA3D" ) &&  (res.trim()=="40KB"  || res.trim()=="25KB" || res.trim()=="10KB")){
					//alert ("The selected Resolution is too high for Algorithm. No structure available ");					
					 window.showAlert1();
				}
				else if ((algo.trim()=="Pastis" || algo.trim()=="Chromosome3D" || algo.trim()=="ChromSDE" || algo.trim()=="GEM" || algo.trim()=="HSA" || algo.trim()=="InfMod3DGen") &&  (res.trim()=="100KB" )){
					// alert ("The selected Resolution is too high for Algorithm. No structure available ");
					  window.showAlert1();
					// details: takes more than 72 hours 
				}
				else if ((algo.trim()=="ChromSDE" ) &&  (res.trim()=="500KB" ) && (chrpdb =="chr1" || chrpdb =="chr2" || chrpdb =="chr3" || chrpdb =="chr4" || chrpdb =="chr5" || chrpdb =="chr6" || chrpdb =="chr7" || chrpdb =="chr8" || chrpdb =="chr9")){
					//alert ("Structure not available for Chromosome. When Algorithm = ChromSDE, and Resolution = 500kb. Due to time and memory constraint, Only Chromosome Structure 10-23 is available");
					 window.showAlert2();
					// details: takes more than 72 hours 
				}
				else if ((algo.trim()=="HSA" ) &&  (res.trim()=="250KB" ) && (chrpdb =="chr1" || chrpdb =="chr2" || chrpdb =="chr3" || chrpdb =="chr4" || chrpdb =="chr5" || chrpdb =="chr6" || chrpdb =="chr7" || chrpdb =="chr8" || chrpdb =="chr9")){
					// alert ("Structure not available for Chromosome. When Algorithm = HSA, and Resolution = 250kb. Due to time and memory constraint, Only Chromosome Structure 10-23 is available");
					window.showAlert2();
					// details: takes more than 72 hours 
				}				
				else if ((algo.trim()=="ChromSDE" ) &&  (res.trim()=="250KB" ) ){
					//alert ("The selected Resolution is too high for Algorithm. No structure available");
					 window.showAlert1();
					// details: takes more than 72 hours 
				}
				else if((gsdb_id.trim()=="TE1402WS" || gsdb_id.trim()=="VH2561BL") && (chrpdb =="chr21" || chrpdb =="chr22" || chrpdb =="chr23" || chrpdb =="chr24" || chrpdb =="chrX" ||chrpdb =="chrY")){
					window.showAlert3();
				}
				else if((gsdb_id.trim()=="BB8015WF" || gsdb_id.trim()=="GG6098MH"|| gsdb_id.trim()=="OO7429SF") && (chrpdb =="chr24" || chrpdb =="chrY")){
					window.showAlert4();
				}
				else if((gsdb_id.trim()!="OO7429SF") && (algo.trim()=="InfMod3DGen")){
					window.showAlert5();
				}
				else if(fname.trim()=="GM06990_Genome" && algo.trim()!="LorDG" && algo.trim()!="MOGEN" && algo.trim()!="ShRec3D" && algo.trim()!="miniMDS" ){
					window.showAlert6();
				}

				else
				{
					var link = "";
					norm_VC ="VC";
					norm_YT ="Yaffe_Tanay";
					norm_KR ="KR"
					norm="";
					if (gsdb_id.trim()=="TE1402WS" || gsdb_id.trim()=="VH2561BL" ||gsdb_id.trim()=="BB8015WF" ||gsdb_id.trim()=="GG6098MH"){
						norm =norm_YT;
					}
					else if(gsdb_id.trim()=="OO7429SF"){
						norm = norm_KR;
					}else{
						norm = norm_VC;
					}
					
					
				   if (splitRes.length == 1){
					   link = base_link + gsdb_id + "/" + fname +  "/" + norm + "/" + algo + "/" + chrpdb  ;
					 
				   }
				   else{
						 var res = res.toLowerCase();						
						 link = base_link + gsdb_id + "/" + fname +  "/" + norm  +  "_" + res + "/" + algo + "/" + chrpdb  ;
				   }
				  
				 var pdblink  =  link + ".pdb" ;
				 var mol_link = "http://3dmol.csb.pitt.edu/viewer.html"; 
				 var Struct_display = '';
				 if (algo.trim()=="LorDG" || algo.trim()=="MOGEN" || algo.trim()=="miniMDS" || algo.trim()=="HSA"){
					var Struct_display =  mol_link + "?url=" + pdblink + "&type=pdb&select=all&style=stick:color~crimson,radius~0.045,";
				 }else{
					 var Struct_display =  mol_link + "?url=" + pdblink + "&type=pdb&select=all&style=line:color~crimson";
				 }
							
				 var Struct_Log = link + ".log";
				 
				
				// Display structure and log
				return[ Struct_display, Struct_Log ];				
				
				}
		}
   


   </script>
   
 
	
	<!-- Placeholder Images -->
	<script src="js/holder.min.js"></script>
	
</body>

</html>
