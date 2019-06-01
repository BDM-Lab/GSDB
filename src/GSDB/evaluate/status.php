<?php


	echo "";
	

if (isset($argc)) {
    echo "**********************************************************************************". "<br/>" ;
	echo " Assessing the similarity between Chromosome IF Matrix and Structure          ". "<br/>" ;
	echo " Reference: Oluwatosin Oluwadare, Max Highsmith, and Jianling Cheng        ". "<br/>" ;
	echo " For comments, Please email to: chengji@missouri.edu									". "<br/>" ;
	echo "**********************************************************************************". "<br/><br/>" ;
	
	$data = nl2br(file_get_contents($argv[1]));
	
	echo $data;
	echo "<br/><br/>";
	echo "**********************************************************************************". "<br/>" ;
	
}
else {
	echo "argc and argv disabled\n";
}

?>