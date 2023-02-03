<?php
sleep(10); //let mysql container come online
require_once(__DIR__."/../connection.php");
$dbname = getenv("MYSQL_DATABASE");
$mysqli->query("CREATE DATABASE {$dbname}");
$mysqli->select_db($dbname);
$mysqli->multi_query(file_get_contents(__DIR__.'/gsdb.sql'));