
<?php

	if (!defined('PDO::ATTR_DRIVER_NAME')) {
		echo 'PDO is unavailable<br/>';
	}
	elseif (defined('PDO::ATTR_DRIVER_NAME')) {
		echo 'PDO is available<br/>';
		print_r(PDO::getAvailableDrivers());
	}

  

	$servername = "localhost";
	$username = "root";
	$password = "genomeflow";
	$database = "gsdb";
	 

try {
    $conn = new PDO("mysql:host=$servername;dbname=$database", $username, $password);
    // set the PDO error mode to exception
    $conn->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);
    echo "Connected successfully !!!\n";
  
	echo " When you get this page it means connection to database was successfully   \n"; 
	
    }
catch(PDOException $e)
    {
    echo "Connection failed: " . $e->getMessage();
    }

?>