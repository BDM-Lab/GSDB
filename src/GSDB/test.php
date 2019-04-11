<?php


$mysqli = new mysqli("localhost", "root" , "genomeflow", "gsdb");

/* check connection */
if (mysqli_connect_errno()) {
    printf("Connect failed: %s\n", mysqli_connect_error());
    exit();
}

$query = "SELECT ID,Title,Biosample_Type,Organism  FROM general_info";
if ($stmt = $mysqli->prepare($query)) {

    /* execute query */
    $stmt->execute();

    /* store result */
    $stmt->store_result();

    printf("Number of rows: %d.\n", $stmt->num_rows);

   /* Bind the result to variables */
   $stmt->bind_result($ID,$Title,$Biosample_Type,$Organism);
	while ($stmt->fetch()) {
        echo 'ID: '.$ID.'<br>';
        echo 'Title: '.$Title.'<br>';
        echo 'Biosample_Type: '.$Biosample_Type.'<br>';
        echo 'Organism: '.$Organism.'<br><br>';
		
		$query2 = "SELECT Filename,Resolution FROM data_info WHERE GSDB_ID = ?";
		if ($stmt2 = $mysqli->prepare($query2)){
			
			$stmt2->bind_param("s", $ID);
			
			/* execute query */
			$stmt2->execute();
			
			/* Store the result (to get properties) */
			$stmt2->store_result();

		   /* Get the number of rows */
		   $num_of_rows2 = $stmt2->num_rows;

		   printf("Number of Inner rows: %d.\n", $stmt2->num_rows);
			
			
			/* Bind the result to variables */
		  $stmt2->bind_result($Filename,$Resoluiton);
		  
		  while ($stmt2->fetch()) {
				echo 'Filename: '.$Filename.'<br>';
				echo 'Resoluiton: '.$Resoluiton.'<br>';
				
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