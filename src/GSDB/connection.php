<?php

$db_server = getenv("MYSQL_HOSTNAME");
$db_username = getenv("MYSQL_USERNAME");
$db_password = getenv("MYSQL_ROOT_PASSWORD");
$db_database = getenv("MYSQL_DATABASE");

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
    
$mysqli = new mysqli($db_server,$db_username,$db_password,$db_database);

/* check connection */
if (mysqli_connect_errno()) {
  printf("Connect failed: %s\n", mysqli_connect_error());
  exit();
}
