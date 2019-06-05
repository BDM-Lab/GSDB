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
    <link rel="icon" href="images/icon.PNG">
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

     <!-- Navigation -->
    <nav id="topNav" class="navbar navbar-full navbar-static-top navbar-dark bg-inverse m-b-1">
        <button class="navbar-toggler hidden-md-up pull-center" type="button" data-toggle="collapse" data-target="#navbar">
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
 
            </ul>
        </div>
    </nav>
    

</head>

<body>

<div class="container-fluid">	

		<h3 class="card-title">	<i class="fa fa-bar-chart" aria-hidden="true"></i> 
			Structure Evaluation 
		</h3>
		<hr> 	
	  <p><font size="3">Compare a chromosome or genome structure with an interaction frequency (IF)  matrix or another chromosome or genome structure.</font><p> 	
	<?php

	//echo ini_get("upload_max_filesize");
	 
	 if(isset($_POST['submit']))
    {
	    
	    $path = "/var/www/html/3dgenome/GSDB/evaluate/jobs/";
		
		//=========================================================================
		// Check if string is empty
		//=========================================================================
	    // function to check if string is empty ot not
	   function strictEmpty($var)
	   {		   
			// Delete this line if you want space(s) to count as not empty
			$var = trim($var);			
			if(isset($var) === true && $var === '') {		
				// It's empty
				return true;
			}
			else {			
				// It's not empty
				return false;
			}
		}
		
		function xrmdir($dir) {
			$items = scandir($dir);
			foreach ($items as $item) {
				if ($item === '.' || $item === '..') {
					continue;
				}
				$path = $dir.'/'.$item;
				if (is_dir($path)) {
					xrmdir($path);
				} else {
					unlink($path);
				}
			}
			rmdir($dir);
		}
		
		//=========================================================================
		//History document: contains list of access by users
		//=========================================================================
		// Function to get the client ip address
		function get_client_ip_env() {
			$ipaddress = '';
			if (getenv('HTTP_CLIENT_IP'))
				$ipaddress = getenv('HTTP_CLIENT_IP');
			else if(getenv('HTTP_X_FORWARDED_FOR'))
				$ipaddress = getenv('HTTP_X_FORWARDED_FOR');
			else if(getenv('HTTP_X_FORWARDED'))
				$ipaddress = getenv('HTTP_X_FORWARDED');
			else if(getenv('HTTP_FORWARDED_FOR'))
				$ipaddress = getenv('HTTP_FORWARDED_FOR');
			else if(getenv('HTTP_FORWARDED'))
				$ipaddress = getenv('HTTP_FORWARDED');
			else if(getenv('REMOTE_ADDR'))
				$ipaddress = getenv('REMOTE_ADDR');
			else
				$ipaddress = 'UNKNOWN';
		 
			return $ipaddress;
		}	
		


	   
	    /*make directory using remote address
	    $ip= get_client_ip_env() ;
		$path = $path.$ip;
		$old_mask = umask(0);
	    mkdir($path);
	  	umask($old_mask);	
		*/
		
		// make directory for job using Id inside Remote address
	     $jobid =   substr(md5(microtime()),rand(0,26),7); // $_POST["jobid"];		 
		 $path  =  $path.$jobid.'/';
		
		 if (file_exists($path )) {			
			xrmdir($path);
		  } 
		 
		 $old_mask = umask(0);
		 mkdir($path);
		 umask($old_mask);
		 
		
		 
		//make directory for input data 1
		  $input_data_1= $path.'input1/';
		  if (file_exists($input_data_1 )) {			
			xrmdir($input_data_1);
		  } 
		  $old_mask = umask(0);
		   mkdir( $input_data_1);  
		   umask($old_mask);
		  
		  
	    //make directory for input data 2
		  $input_data_2= $path.'input2/';
		  if (file_exists($input_data_2 )) {			
			xrmdir($input_data_2);
		  } 
		  $old_mask = umask(0);
		   mkdir( $input_data_2);  
		   umask($old_mask);  
		  
		 //==================================
		 //Uploading IF Matrix File from url 
		 //==================================
		  $matrix_url='';		  
		 if (!strictEmpty($_POST["ifmatrixlink"]))
		 {
			 $matrix_url=$_POST['ifmatrixlink'];
			 $data = file_get_contents($matrix_url);			
			 $matrix_url = $input_data_1.basename($matrix_url);
			 $upload =file_put_contents($matrix_url, $data);
			 if($upload) {
				//echo "matrix Upload Successful";
			 }else{
				//echo "matrix Upload Not successfull";
			 } 
		}
	
	   //==================================
	   // Upload pdb Structure from url 
	   //==================================
	     $pdb_url='';
		 if(  !strictEmpty($_POST["pdblink"]))
		 {
			 $pdb_url=$_POST['pdblink'];
			 $data = file_get_contents($pdb_url);			
			 $pdb_url = $input_data_2.basename($pdb_url);		 
			 $upload =file_put_contents($pdb_url, $data);
			 if($upload) {
				// echo "pdb Upload Successful";
			 }else{
				//echo "pdb upload Not successfull";
			 } 
		}		 
		 
		 
		 //==================================
		 //Uploading IF Matrix File
		 //==================================
		if(!empty($_FILES['ifmatrixfile']))
		  {
			//$path = "/var/www/html/3dgenome/GSDB/evaluate/jobs/";			
			$ifpath =  $input_data_1 . basename( $_FILES['ifmatrixfile']['name']);
			if(move_uploaded_file($_FILES['ifmatrixfile']['tmp_name'], $ifpath)) {		
			
			 // echo "The IF matrix file ".  basename( $_FILES['ifmatrixfile']['name']). " has been uploaded"; echo "<br>";
			 $matrix_file = $ifpath ;
			} else{
			 // echo "There was an error uploading the IF matrix file, please try again!"; echo "<br>";
			}
		  }	
		
		//==================================
		//Uploading the pdb structure File
	    //================================== 
		 if(!empty($_FILES['pdbfile']))
		  {
			//$path = "/var/www/html/3dgenome/GSDB/evaluate/jobs/";
			$pdbpath =  $input_data_2 . basename( $_FILES['pdbfile']['name']);
			if(move_uploaded_file($_FILES['pdbfile']['tmp_name'], $pdbpath)) {
			 // echo "The pdb structure file ".  basename( $_FILES['pdbfile']['name']).  " has been uploaded"; echo "<br>";
			  $pdb_file = $pdbpath;
			} else{
			  //echo "There was an error uploading the pdb structure file, please try again!"; echo "<br>";
			}
		  }	
		
		//Select either link or file that has been uploaded
		$matrix_file= strictEmpty($matrix_url)?$matrix_file:$matrix_url;
		$pdb_file= strictEmpty($pdb_url)?$pdb_file:$pdb_url;
		
		
		// echo $matrix_file; 	 echo $pdb_file;
		
		// Record Access in a history file			
		$history =   "/var/www/html/3dgenome/GSDB/evaluate/jobs/history.txt";
		$jd = fopen($history , 'a');
		$ipaddress = "IP: ". get_client_ip_env().PHP_EOL;		
	    fwrite($jd, $ipaddress);
		$Time= "Time: " . date("F j, Y, g:i a").PHP_EOL	;	
		fwrite($jd, $Time);
		$ifmatrix = "IF Matrix: ".$matrix_file.PHP_EOL;
		fwrite($jd , $ifmatrix);
		$id = "Job ID: ". $jobid .PHP_EOL;
		fwrite($jd , $id);
		$space= "==============================================================".PHP_EOL;		
		fwrite($jd , $space);		
		fclose($jd);
	  
	  
	    // create a shell script to execute java code 
		$job_log = 	  $path.$jobid.'_log.txt';		
	    //Decide based on the check box selected the shell script to call_user_func
		$type = $_POST['input_type'];  
		if ($type == "if") {          
			//echo 'IF_Structure Script to be called'; 
			shell_exec("/var/www/html/3dgenome/GSDB/evaluate/IF_Structure_evaluate.sh $matrix_file $pdb_file $path $jobid > $job_log &");		
		}
		else {
			//echo 'Structure_Structure Script to be called';
			shell_exec("/var/www/html/3dgenome/GSDB/evaluate/Structure_Structure_evaluate.sh $matrix_file $pdb_file $path $jobid > $job_log &");		
		}
	    
		// print output to iframe	
		$data = '';
		$time_elapse = 0;
		$result =   "/var/www/html/3dgenome/GSDB/evaluate/jobs/".$jobid."/".$jobid."_Result.log";
		while($data == '')
		{
			sleep(1);
			$time_elapse++;			
			$data = nl2br(file_get_contents($result));
			if($time_elapse>20)
			{
				$data = "Evaluation failed. Please check the input and try again!";
				
			}
			if($data != '' || $time_elapse>20)
			{
				break;
			}
		}
		

	    $url="evaluate_Result.php?output=$result";		
        header("Location:$url");
	

	}
		
	?>	
	<!-- Left Column-->
	<div class="col-sm-12">	

				<style type="text/css">
				  
				  td { font-size: 17px }
				</style>		
				
				<form id="form1" name="form1" method="POST" onsubmit="return validateForm()" action="" enctype="multipart/form-data" ">		
					<ul>
					 <li> <b>Type of input file 1</b> (required): <br>
						  <input type="radio" name="input_type" id="input_type_if" value="if" checked="checked"> Interaction frequency(IF) matrix <br> 
						  <input type="radio" name="input_type" id="input_type_pdb" value="pdb"> Chromosome structure in <a href="http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html"><u>PDB format</u></a><br>
					<br>				
					<li> <b>Input file 1</b> (required): <br>
						Please copy and paste the URL of input file 1 here. <a href="javascript:copy_link_1()"><u>IF Matrix Sample input</u> </a> &nbsp;  <a href="javascript:copy_link_2()"><u>Structure Sample input</u> </a><br>
						<textarea name="ifmatrixlink" id="ifmatrixlink" name="ifmatrixlink" rows="3" cols="100" ></textarea><br>
						Or upload the input file 1: <br><input type="file" id="ifmatrixfile" name="ifmatrixfile" size="70"> <br>
						
					<br>	
					<li> <b>Input chromosome structure file 2 in <a href="http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html"><u>PDB format</b> </u></a> (required): <br>	
						 Please copy and paste the URL of your structure file here. <a href="javascript:copy_link_3()"><u>Structure Sample input</u> </a><br>
						 <textarea id="pdblink" name="pdblink" rows="3" cols="100" ></textarea></td><br>
						 Or upload the input chromosome structure file 2:<br>
						 <input type="file" id="pdbfile" name="pdbfile" size="70"> <br>
						 
						 <br><br>
						 <li> <input name="submit" type="submit" style="font-size:16px" value="Compare" onclick="func()"  />  &nbsp; &nbsp;
						 <INPUT TYPE="reset" style="font-size:16px" VALUE="Clear form">
					</ul>
				</form>				
					<br/>
				 
	 
		 <br/>
			
	</div>	
	
	<!-- center Column-->
	
	
		
	
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
		
    function copy_link_1()
	{
		var xvalue=IO("evaluate/file1_example.txt");
		document.form1.ifmatrixlink.value = xvalue;
	}
	function copy_link_2()
	{
		var yvalue=IO("evaluate/file2_example.txt");
		document.form1.ifmatrixlink.value = yvalue;
	}
	
	function copy_link_3()
	{
		var zvalue=IO("evaluate/file3_example.txt");
		document.form1.pdblink.value = zvalue;
	}
	function IO(U, V) 
	{
		var X = !window.XMLHttpRequest ? new ActiveXObject('Microsoft.XMLHTTP') : new XMLHttpRequest();
		X.open(V ? 'PUT' : 'GET', U, false );
		X.setRequestHeader('Content-Type', 'text/html')
		X.send(V ? V : '');
		return X.responseText;
	}
	
	
	
   </script> 

   <!----CDN DataTable -->
	<script type="text/javascript" charset="utf-8">			
			$(document).ready( function () {
				 $('#mytable1').DataTable({	
					"searching": false
				} );
			} );
					
	</script>
	<script type="text/javascript" charset="utf-8">

		function getFileExtension(filename) {
		  return filename.slice((filename.lastIndexOf(".") - 1 >>> 0) + 2);
		}
		
		
		function validateForm(){
			//var jobid = document.getElementById("jobid").value;
			//var email = document.getElementById("email").value;
			var ifmatrixlink = document.getElementById("ifmatrixlink").value;
			var ifmatrixfile = document.getElementById("ifmatrixfile").value;
			var pdblink = document.getElementById("pdblink").value;
			var pdbfile = document.getElementById("pdbfile").value;
			
			/*if(email == null || email == ""){
				alert("Please supply email address!");
				return false;
			}
			
			if(jobid == null || jobid == ""){
				alert("Please supply job id!");
				return false;
			}*/
			if((ifmatrixlink == null || ifmatrixlink == "") && (ifmatrixfile == null || ifmatrixfile == "")){
				alert("Please provide a input file - IF matrix or pdb structure!");
				return false;
			}
			if(!(ifmatrixlink == null || ifmatrixlink == "") && !(ifmatrixfile == null || ifmatrixfile == "")){
				alert("Error for input file 1: \n Please supply one option: either paste a file link OR upload a file - not both!");
				document.getElementById("ifmatrixlink").value = '';
				return false;
			}
			if((pdblink == null || pdblink == "") && (pdbfile == null || pdbfile == "")){
				alert("Please supply a pdb structure file!");
				return false;
			}
			
			if(!(pdblink == null || pdblink == "") && !(pdbfile == null || pdbfile == "")){
				alert("Please supply one option: either paste pdb file link OR upload a pdb file - not both!");
				return false;
			}
			
			// validate the file extnetion for pdb file			 
			var allowedExtensions = /(\.pdb)$/i;
			
					
		   if (document.getElementById("input_type_pdb").checked == true && !allowedExtensions.exec(ifmatrixfile) && !(ifmatrixfile == null || ifmatrixfile == "")  ) {
				alert('Pdb structure option checked. \n Please upload a file having extension ".pdb" for the "Upload selected option file" field.');
				document.getElementById("ifmatrixfile").value = '';
				return false;
			}
			
			
			if (document.getElementById("input_type_pdb").checked == true && !allowedExtensions.exec(ifmatrixlink) && !(ifmatrixlink == null || ifmatrixlink == "")  ) {
				alert('Pdb structure option checked. \n Please upload a file having extension ".pdb" for the "Input type file link" field.');
				document.getElementById("ifmatrixlink").value = '';
				return false;
			}			
			
			
			if(!allowedExtensions.exec(pdbfile) && !(pdbfile == null || pdbfile == "") ){
				alert('Please upload file having extension .pdb  for the "Upload a pdb structure file" field.');
				document.getElementById("pdbfile").value = '';
				return false;
			}
			
		   if(!allowedExtensions.exec(pdblink) && !(pdblink == null || pdblink == "")  ){
				alert('Please upload file having extension .pdb for the "Pdb structure file" field.');
				document.getElementById("pdblink").value = '';
				return false;
			}

		}
		
	
	</script>

	
	
	
	
</body>

</html>
