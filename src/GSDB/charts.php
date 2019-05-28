<?php			
			require_once('connection.php');
			
			/*======================================================================
			 LAB
			======================================================================*/
			$result=$conn->prepare("SELECT Lab, COUNT(*) AS freq FROM general_info GROUP BY Lab");
			$result->execute();
					
			 // The `$arrData` array holds the chart attributes and data
            $arrData_1 = array(
                "chart" => array(
                    "caption" => " Hi-C data Contributors",                    
					"plotgradientcolor"=> "",
					"bgcolor"=> "FFFFFF",
					"showalternatehgridcolor"=> "0",
					"showplotborder"=> "0",
					"divlinecolor"=> "CCCCCC",
					"showvalues"=> "1",
					"showcanvasborder"=> "0",
					"canvasbordercolor"=> "CCCCCC",
					"canvasborderthickness"=> "1",
					"yaxismaxvalue"=> "30000",
					"captionpadding"=> "30",
					"linethickness"=> "3",
					"sshowanchors"=> "0",
					"yaxisvaluespadding"=> "15",
					"showlegend"=> "1",
					"use3dlighting"=> "0",
					"showshadow"=> "0",
					"legendshadow"=> "0",
					"legendborderalpha"=> "0",
					"showborder"=> "0",
					"palettecolors"=> "#EED17F,#97CBE7,#074868,#B0D67A,#2C560A,#DD9D82",
					"bgAlpha"=> "0",
					"borderColor"=> "#666666",
					"borderThickness"=>  "4"
                  )
               );
			 $arrData_1 ["data"] = array();
			
			 // Push the data into the array
			for($i=0; $row = $result->fetch(); $i++){
				array_push($arrData_1 ["data"], array(
                "label" => $row[0],
                "value" => $row[1],               
                 )
               );	
			 
			}
			
			
			/*JSON Encode the data to retrieve the string containing the JSON representation of the data in the array. */
         
			$jsonEncodedData = json_encode($arrData_1 );				
			/*Create an object for the column chart using the FusionCharts PHP class constructor. Syntax for the constructor is ` FusionCharts("type of chart", "unique chart id", width of the chart, height of the chart, "div id to render the chart", "data format", "data source")`. Because we are using JSON data to render the chart, the data format will be `json`. The variable `$jsonEncodeData` holds all the JSON data for the chart, and will be passed as the value for the data source parameter of the constructor.*/
            $doughnutChart = new FusionCharts("doughnut2D", "Lab-Chart" , 350, 300, "chart-2", "json", $jsonEncodedData);
          	// Render the chart
            $doughnutChart->render();
		
			/*======================================================================
			//RESOLUTION
			======================================================================*/
			$result=$conn->prepare("SELECT Resolution, COUNT(*) AS freq FROM data_info GROUP BY Resolution");
			$result->execute();	
			 // The `$arrData` array holds the chart attributes and data
            $arrData_2 = array(
                "chart" => array(
                    "caption" => "Resolution",                    
					"showlabels"=> "0",
					"showlegend"=> "1",
					"enablemultislicing"=> "0",
					"slicingdistance"=> "15",
					"showpercentvalues"=> "1",
					"showpercentintooltip"=> "0",
					"bgAlpha"=> "0",
					
                  )
               );
			 $arrData_2 ["data"] = array();
			 // Push the data into the array
			for($i=0; $row = $result->fetch(); $i++){
				array_push($arrData_2 ["data"], array(
                "label" => $row[0],
                "value" => $row[1],               
                 )
               );	
			}
			/*JSON Encode the data to retrieve the string containing the JSON representation of the data in the array. */
            $jsonEncodedData = json_encode($arrData_2 );	
			/*Create an object for the column chart using the FusionCharts PHP class constructor. Syntax for the constructor is ` FusionCharts("type of chart", "unique chart id", width of the chart, height of the chart, "div id to render the chart", "data format", "data source")`. Because we are using JSON data to render the chart, the data format will be `json`. The variable `$jsonEncodeData` holds all the JSON data for the chart, and will be passed as the value for the data source parameter of the constructor.*/
            $doughnutChart = new FusionCharts("doughnut3D", "Resolution-Chart" , 370, 300, "chart-3", "json", $jsonEncodedData);
            // Render the chart
            $doughnutChart->render();		

			/*======================================================================
			//BIOSAMPLE
			======================================================================*/
			$cell = "cell line";
			$prim= "primary cell";
			$x=0;   //cell line
			$y=0;   //primary cell
			
			$mysqli = new mysqli("localhost", "root" , "genomeflow", "gsdb");
			/* check connection */
			if (mysqli_connect_errno()) {
				printf("Connect failed: %s\n", mysqli_connect_error());
				exit();
			}

			$query = "SELECT ID,Biosample_Type  FROM general_info";
			if ($stmt = $mysqli->prepare($query)) {
				/* execute query */
				$stmt->execute();
				/* store result */
				$stmt->store_result();
			   /* Bind the result to variables */
			   $stmt->bind_result($ID,$Biosample_Type);
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
							if (strcmp($cell,$Biosample_Type)==0){
								$x=$x+1;
							}
							if (strcmp($prim,$Biosample_Type)==0){
								$y=$y+1;
							}
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
			 // The `$arrData` array holds the chart attributes and data
            $arrData_3 = array(
                "chart" => array(
                    "caption" => "Biosample Type",                    
					"plotgradientcolor"=> "",
					"bgcolor"=> "FFFFFF",
					"showalternatehgridcolor"=> "0",
					"showplotborder"=> "0",
					"divlinecolor"=> "CCCCCC",
					"showvalues"=> "1",
					"showcanvasborder"=> "0",
					"canvasbordercolor"=> "CCCCCC",
					"canvasborderthickness"=> "1",
					"yaxismaxvalue"=> "30000",
					"captionpadding"=> "30",
					"linethickness"=> "3",
					"sshowanchors"=> "0",
					"yaxisvaluespadding"=> "15",
					"showlegend"=> "1",
					"use3dlighting"=> "0",
					"showshadow"=> "0",
					"legendshadow"=> "0",
					"legendborderalpha"=> "0",
					"showborder"=> "0",
					"palettecolors"=> "#EED17F,#97CBE7,#074868,#B0D67A,#2C560A,#DD9D82",
					"bgAlpha"=> "0",
					"borderColor"=> "#666666",
					"borderThickness"=>  "4"
                  )
               );
			
			 $arrData_3 ["data"] = array();
			 // Push the data into the array
			 array_push($arrData_3 ["data"], array("label" => $cell, "value" => $x, 	)             );
			 array_push($arrData_3 ["data"], array("label" => $prim, "value" => $y, 	)             );
			
			/*JSON Encode the data to retrieve the string containing the JSON representation of the data in the array. */
            $jsonEncodedData = json_encode($arrData_3 );	
			/*Create an object for the column chart using the FusionCharts PHP class constructor. Syntax for the constructor is ` FusionCharts("type of chart", "unique chart id", width of the chart, height of the chart, "div id to render the chart", "data format", "data source")`. Because we are using JSON data to render the chart, the data format will be `json`. The variable `$jsonEncodeData` holds all the JSON data for the chart, and will be passed as the value for the data source parameter of the constructor.*/
            $doughnutChart = new FusionCharts("doughnut2D", "Biosmaple-Chart" , 350, 300, "chart-4", "json", $jsonEncodedData);
            // Render the chart
            $doughnutChart->render();
			
						require_once('connection.php');
			/*======================================================================
			PROJECT
			======================================================================*/
			$homo = "ENCODE";
			$homoM ="GGR";
			$Mus= "Unknown";
			$x=0;   
			$y=0;   
			$z=0 ;    
			$mysqli = new mysqli("localhost", "root" , "genomeflow", "gsdb");
			/* check connection */
			if (mysqli_connect_errno()) {
				printf("Connect failed: %s\n", mysqli_connect_error());
				exit();
			}

			$query = "SELECT ID,Project  FROM general_info";
			if ($stmt = $mysqli->prepare($query)) {
				/* execute query */
				$stmt->execute();
				/* store result */
				$stmt->store_result();
			   /* Bind the result to variables */
			   $stmt->bind_result($ID,$Project);
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
							if (strcmp($homo,$Project)==0){
								$x=$x+1;
							}
							if (strcmp($homoM,$Project)==0){
								$y=$y+1;
							}
							if (strcmp($Mus,$Project)==0){
								$z=$z+1;
							}
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
			
			 // The `$arrData` array holds the chart attributes and data
            $arrData = array(
                "chart" => array(
                    "caption" => "Hi-C data Source Project ",
					"palette"=> "2",
					"animation"=> "1",
					"formatnumberscale"=> "1",					
					"numberprefix"=> "$",
					"pieslicedepth"=> "10",
					"startingangle"=> "35",
					"showborder"=> "0",
                    "showValues"=> "0",					
                    "theme"=> "zune",
					"bgAlpha"=> "0",
					
                  )
               );
			 $arrData["data"] = array();
			 // Push the data into the array
			// Push the data into the array
			 array_push($arrData["data"], array("label" => $homo, "value" => $x,  "issliced"=> "1"	)             );
			 array_push($arrData["data"], array("label" => $homoM, "value" => $y,  "issliced"=> "0"	)             );
			 array_push($arrData["data"], array("label" => $Mus, "value" => $z,  "issliced"=> "0"	)             );
			 
			/*JSON Encode the data to retrieve the string containing the JSON representation of the data in the array. */
            $jsonEncodedData = json_encode($arrData);	
			/*Create an object for the column chart using the FusionCharts PHP class constructor. Syntax for the constructor is ` FusionCharts("type of chart", "unique chart id", width of the chart, height of the chart, "div id to render the chart", "data format", "data source")`. Because we are using JSON data to render the chart, the data format will be `json`. The variable `$jsonEncodeData` holds all the JSON data for the chart, and will be passed as the value for the data source parameter of the constructor.*/
            $columnChart = new FusionCharts("doughnut3D", "Project-Chart" , 350, 300, "chart-1", "json", $jsonEncodedData);
            // Render the chart
            $columnChart->render();
			
	/**===============================================================================================
	 * Supplementary json_encode in case php version is < 5.2 (taken from http://gr.php.net/json_encode)
	 */

    function json_encode($a=false)
    {
        if (is_null($a)) return 'null';
        if ($a === false) return 'false';
        if ($a === true) return 'true';
        if (is_scalar($a))
        {
            if (is_float($a))
            {
                // Always use "." for floats.
                return floatval(str_replace(",", ".", strval($a)));
            }

            if (is_string($a))
            {
                static $jsonReplaces = array(array("\\", "/", "\n", "\t", "\r", "\b", "\f", '"'), array('\\\\', '\\/', '\\n', '\\t', '\\r', '\\b', '\\f', '\"'));
                return '"' . str_replace($jsonReplaces[0], $jsonReplaces[1], $a) . '"';
            }
            else
            return $a;
        }
        $isList = true;
        for ($i = 0, reset($a); $i < count($a); $i++, next($a))
        {
            if (key($a) !== $i)
            {
                $isList = false;
                break;
            }
        }
        $result = array();
        if ($isList)
        {
            foreach ($a as $v) $result[] = json_encode($v);
            return '[' . join(',', $result) . ']';
        }
        else
        {
            foreach ($a as $k => $v) $result[] = json_encode($k).':'.json_encode($v);
            return '{' . join(',', $result) . '}';
        }
    }

			
	?>