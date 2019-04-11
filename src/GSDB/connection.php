<?php

$db_server = "localhost";
$db_username = "root";
$db_password = "genomeflow";
$db_database = "gsdb";


try {
	$conn = new PDO("mysql:host=$db_server;dbname=$db_database", $db_username, $db_password);
	$conn->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);
	
   //echo "Connected successfully !!!\n";
    }
catch(PDOException $e)
    {
    echo "Connection failed: " . $e->getMessage();
  //	echo " When you get this page it means connection to database was Unsuccessfully   \n"; 
    }
?>
